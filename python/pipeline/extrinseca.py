# -*- coding: utf-8 -*-
"""Estimación de la extrínseca escáner <-> trayectoria (LAS-1.EP es propietario y no se decodifica).

HALLAZGO (b1): el Hesai XT32 del RS30 GIRA LENTAMENTE sobre su propio eje y (30 rpm, periodo ~2 s):
la normal del suelo vista en el marco del escáner describe un cono alrededor de y_s (mediana n.a = 0,42,
a = (0,01; 1,00; -0,02)). Por tanto la extrínseca depende del tiempo y no basta una transformación fija.
Modelo usado:

  X_mundo = Rw (X_ref - c) + c + dW ,   X_ref = p(t) + R(t) [ Er Rot(a_s, -theta(t)) (k r u(az, canal)) + Et ]
  theta(t) = omega_w (t - t_w) + phi_w     (una pareja omega, phi por ventana de tiempo)
  Er, Et : base de giro -> cuerpo de la trayectoria ; a_s: eje de giro en el escáner
  dW, Rw : corrección rígida del marco local->UTM (el ajuste Kabsch vs PosI tiene 0,18 m de rms)
  k      : factor de escala de distancia (diagnóstico 4 mm / 5 mm)

Pasos: (1) ángulo de giro por ventana a partir de la normal del suelo (RANSAC en el marco del escáner),
(2) inclinación de Er por gravedad y guiñada psi por búsqueda en 360° con correspondencias EXACTAS por
grupo temporal (el gps_time del LAS está truncado a 1/1024 s y cada sello agrupa los disparos de ese
intervalo), (3) refinamiento por mínimos cuadrados de todo el conjunto.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation as R

import pcap_hesai as ph

CUANTO_T = 1.0 / 1024.0
REL_T_GRUPO = 0.00015        # el sello k*CUANTO_T agrupa disparos en [k*q - 0,15 ms, k*q + 0,83 ms) (medido por correlación)


# --------------------------------------------------------------------------------------
# Modelo geométrico
# --------------------------------------------------------------------------------------
def tasa_giro(tray, t, dt=0.1):
    """Velocidad angular del cuerpo (grados/s) en los instantes t (diferencia de orientación en dt)."""
    r1 = tray.pose(t)[1]
    r2 = tray.pose(t + dt)[1]
    return np.degrees((r1.inv() * r2).magnitude()) / dt


def _widx(par, t):
    return np.clip(np.searchsorted(par["w_ini"], t, side="right") - 1, 0, len(par["w_ini"]) - 1)


def _theta(par, t):
    w = _widx(par, t)
    th = par["omega"][w] * (t - par["w_ini"][w]) + par["phi"][w]
    kt = par.get("kn_t")
    if kt is not None and len(kt):
        th = th + np.interp(t, kt, par["kn_d"])      # irregularidades de la fase del motor (nudos cada ~1 s)
    return th


def _Rw(par, t):
    """Rotación del mundo (alrededor de 'centro'): rotvec_w + dR*(t - t_ini) (deriva angular de la trayectoria hf
    respecto al LAS, por ventana)."""
    rv = np.tile(np.asarray(par["rotvec_w"], float), (len(t), 1))
    dr = par.get("dR")
    if dr is not None:
        w = _widx(par, t)
        rv = rv + np.asarray(dr)[w] * (t - par["w_ini"][w])[:, None]
    return R.from_rotvec(rv)


def rot_cuerpo(par, t, tray):
    """Rotación escáner -> mundo por instante: Rw * R(t) * Er * Rot(a_s, -theta)."""
    th = _theta(par, t)
    spin = R.from_rotvec(-th[:, None] * par["a_s"][None, :])
    _, rot = tray.pose(t)
    return _Rw(par, t) * rot * R.from_rotvec(par["rotvec_e"]) * spin


def _dW(par, t):
    """Desplazamiento del mundo por ventana: dW + dV*(t - t_ini) (deriva lineal: el LAS sale de una trayectoria SLAM
    optimizada que se separa de la hf unos dm a lo largo de un paso de ~1 min)."""
    d = np.asarray(par["dW"])
    if d.ndim != 2:
        return d
    w = _widx(par, t)
    dv = par.get("dV")
    out = d[w]
    if dv is not None:
        out = out + np.asarray(dv)[w] * (t - par["w_ini"][w])[:, None]
    kt = par.get("kdW_t")
    if kt is not None and len(kt):                    # corrección suave y no lineal (nudos cada ~5 s)
        kd = np.asarray(par["kdW"])
        out = out + np.stack([np.interp(t, kt, kd[:, j]) for j in range(3)], axis=1)
    return out


def origen_mundo(tray, t, par):
    p, rot = tray.pose(t)
    xr = rot.apply(np.tile(par["t_e"], (len(t), 1))) + p
    return _Rw(par, t).apply(xr - par["centro"]) + par["centro"] + _dW(par, t)


def a_mundo(tray, t, ps, par):
    """Puntos del escáner ps (N,3) con tiempo t (N,) -> mundo."""
    th = _theta(par, t)
    spin = R.from_rotvec(-th[:, None] * par["a_s"][None, :])
    p, rot = tray.pose(t)
    xb = R.from_rotvec(par["rotvec_e"]).apply(spin.apply(ps * par["k"])) + par["t_e"]
    xr = rot.apply(xb) + p
    return _Rw(par, t).apply(xr - par["centro"]) + par["centro"] + _dW(par, t)


def rayos_mundo(tray, t_blq, u, par):
    """Origen (n*32,3) y dirección (n*32,3) en el mundo de los disparos de n bloques (u: (n,32,3) en el escáner)."""
    Rt = rot_cuerpo(par, t_blq, tray)
    q = np.repeat(Rt.as_quat(), 32, axis=0)
    d = R.from_quat(q).apply(u.reshape(-1, 3))
    o = np.repeat(origen_mundo(tray, t_blq, par), 32, axis=0)
    return o, d


# --------------------------------------------------------------------------------------
# Paso 1: giro del cabezal a partir de la normal del suelo
# --------------------------------------------------------------------------------------
def _plano_ransac(P, rng, it=300, tol=0.03):
    mejor = (0, None, None)
    for _ in range(it):
        i = rng.choice(len(P), 3, replace=False)
        a, b, c = P[i]
        n = np.cross(b - a, c - a)
        nn = np.linalg.norm(n)
        if nn < 1e-6:
            continue
        n /= nn
        cnt = int((np.abs((P - a) @ n) < tol).sum())
        if cnt > mejor[0]:
            mejor = (cnt, n, a)
    n, a = mejor[1], mejor[2]
    if n @ a > 0:
        n = -n                       # normal hacia el sensor
    return n, mejor[0] / len(P)


def estimar_giro(cfg, tb, az, dist_m, w_ini, w_fin, rng, signo, log=print):
    """Eje a_s, y (omega, phi) por ventana, a partir de la normal del suelo por ventanas de 0,1 s."""
    u = ph.a_direcciones(az * signo)
    Ns, Ts, Ws = [], [], []
    for iw, (a, b) in enumerate(zip(w_ini, w_fin)):
        for t0 in np.arange(a + 0.5, b - 0.6, 0.1):
            m = (tb >= t0) & (tb < t0 + 0.1)
            d = dist_m[m]
            ok = (d > 0.5) & (d < 12)
            if ok.sum() < 2000:
                continue
            P = (u[m] * d[..., None])[ok]
            n, inl = _plano_ransac(P, rng)
            if inl > 0.35:
                Ns.append(n)
                Ts.append(t0 + 0.05)
                Ws.append(iw)
    Ns, Ts, Ws = np.array(Ns), np.array(Ts), np.array(Ws)
    c = Ns.mean(0)
    _, s, vt = np.linalg.svd(Ns - c)
    a_s = vt[2] if vt[2][1] > 0 else -vt[2]
    e1 = np.cross(a_s, [1.0, 0, 0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a_s, e1)
    ang = np.arctan2(Ns @ e2, Ns @ e1)
    omega, phi = np.zeros(len(w_ini)), np.zeros(len(w_ini))
    for iw in range(len(w_ini)):
        m = Ws == iw
        t, th = Ts[m], np.unwrap(ang[m])
        for _ in range(3):
            A = np.c_[t - w_ini[iw], np.ones(len(t))]
            co, *_ = np.linalg.lstsq(A, th, rcond=None)
            r = th - A @ co
            k = np.abs(r) < np.radians(20)
            if k.sum() < 10:
                break
            t, th = t[k], th[k]
        omega[iw], phi[iw] = co
        log("ventana %d: omega=%.4f rad/s (periodo %.3f s), resid sd %.1f deg, n=%d" % (
            iw, co[0], 2 * np.pi / abs(co[0]), np.degrees(r[k].std()) if k.sum() > 2 else -1, k.sum()))
    cono = float((Ns @ a_s).mean())
    th_all = omega[Ws] * (Ts - w_ini[Ws]) + phi[Ws]
    m_ref = np.mean([R.from_rotvec(-th * a_s).apply(n) for th, n in zip(th_all, Ns)], axis=0)
    m_ref /= np.linalg.norm(m_ref)
    return dict(a_s=a_s, omega=omega, phi=phi, m_ref=m_ref, cono_n_a=cono, n_normales=len(Ns))


# --------------------------------------------------------------------------------------
# Grupos temporales exactos
# --------------------------------------------------------------------------------------
def grupos_temporales(las_t, las_xyz, tb, az, dist_m, n_grupos, rng, min_pts=60, excluir=None):
    """Sortea n_grupos grupos temporales (sello k del gps_time del LAS). 'excluir': sellos k que no pueden salir
    (los usados en el ajuste, para que la validación hold-out use solo datos no usados en el ajuste)."""
    k_las = np.round(las_t / CUANTO_T).astype(np.int64)
    ks, cnt = np.unique(k_las, return_counts=True)
    ok = ks[cnt >= min_pts]
    if excluir is not None and len(excluir):
        ok = np.setdiff1d(ok, np.asarray(excluir, dtype=np.int64))
    if len(ok) == 0:
        return None
    ok = rng.choice(ok, size=min(n_grupos, len(ok)), replace=False)
    k_blq = np.floor((tb + REL_T_GRUPO) / CUANTO_T).astype(np.int64)
    o = np.argsort(k_las, kind="stable")
    ks_s = k_las[o]
    out = []
    for k in ok:
        ib = np.nonzero(k_blq == k)[0]
        if len(ib) == 0:
            continue
        bi, ci = np.nonzero(dist_m[ib] > 0.3)
        if len(bi) < 20:
            continue
        a, b = np.searchsorted(ks_s, k, "left"), np.searchsorted(ks_s, k, "right")
        out.append(dict(k=int(k), t=tb[ib][bi], az=az[ib][bi], ch=ci, r=dist_m[ib][bi, ci], P=las_xyz[o[a:b]]))
    if not out:
        return None
    t = np.concatenate([g["t"] for g in out])
    az_ = np.concatenate([g["az"] for g in out])
    ch = np.concatenate([g["ch"] for g in out])
    r = np.concatenate([g["r"] for g in out])
    lim = np.cumsum([0] + [len(g["t"]) for g in out])
    return dict(g=out, t=t, az=az_, ch=ch, r=r, lim=lim, k=np.array([g["k"] for g in out], dtype=np.int64))


def puntos_fuente(G, par):
    """Puntos en el marco del escáner: azimut (con signo y desfase), canal (con desfase de elevación) y distancia."""
    tabla = par.get("elev_deg")
    tabla = ph.ELEV_XT32_DEG if tabla is None else np.asarray(tabla, float)    # tabla de elevaciones (por defecto la estándar)
    el = np.radians(tabla[G["ch"]] + par.get("el_off", 0.0))
    az = np.radians(G["az"] * par["signo"] + par.get("az_off", 0.0))
    ce = np.cos(el)
    u = np.stack([ce * np.sin(az), ce * np.cos(az), np.sin(el)], axis=-1)
    return u * G["r"][:, None]


def residuos_grupos(G, tray, par, tope=0.5, por_grupo=False):
    """Para cada punto del LAS del grupo, distancia al eco gemelo del pcap (el LAS sólo conserva
    una parte de los ecos, así que el sentido LAS -> pcap es el que tiene correspondencia exacta)."""
    W = a_mundo(tray, G["t"] + par.get("dt", 0.0), puntos_fuente(G, par), par)
    d = []
    for g, a, b in zip(G["g"], G["lim"][:-1], G["lim"][1:]):
        dd, _ = cKDTree(W[a:b]).query(g["P"], distance_upper_bound=tope)
        d.append(np.minimum(dd, tope))
    return d if por_grupo else np.concatenate(d)


# --------------------------------------------------------------------------------------
# Calibración completa
# --------------------------------------------------------------------------------------
def _sub_grupos(G, sel):
    gs = [G["g"][i] for i in sel]
    cat = lambda k: np.concatenate([g[k] for g in gs])
    lim = np.cumsum([0] + [len(g["t"]) for g in gs])
    return dict(g=gs, t=cat("t"), az=cat("az"), ch=cat("ch"), r=cat("r"), lim=lim)


def refinar_fase(G, tray, par, w_ini, w_fin, paso_nudo=1.0, log=print, pasadas=2):
    """Fase del giro del cabezal: corrección delta(t) lineal a trozos con un nudo cada 'paso_nudo' s,
    ajustada nudo a nudo por búsqueda 1D del residuo medio (tope 0,3 m) con los grupos vecinos."""
    kn = np.concatenate([np.arange(a, b + 1e-6, paso_nudo) for a, b in zip(w_ini, w_fin)])
    par = dict(par)
    par["kn_t"], par["kn_d"] = kn, np.zeros(len(kn))
    gt = np.array([g["t"].mean() for g in G["g"]])
    for ps_ in range(pasadas):
        for i, tk in enumerate(kn):
            sel = np.nonzero(np.abs(gt - tk) < paso_nudo)[0]
            if len(sel) < 4:
                continue
            Gs = _sub_grupos(G, sel)
            base = par["kn_d"][i]
            mejor = (9, base)
            for dd in np.radians(np.arange(-6, 6.01, 0.5 if ps_ == 0 else 0.25)) + (base if ps_ == 0 else 0):
                par["kn_d"][i] = dd if ps_ == 0 else base + (dd - 0)
                sc = float(residuos_grupos(Gs, tray, par, 0.3).mean())
                if sc < mejor[0]:
                    mejor = (sc, par["kn_d"][i])
            par["kn_d"][i] = mejor[1]
        d = residuos_grupos(G, tray, par, 0.5)
        log("fase pasada %d: mediana %.3f m, p90 %.3f m" % (ps_ + 1, np.median(d), np.percentile(d, 90)))
    return par


def escanear_fase_ventanas(G, tray, par, log=print, etapas=((3.0, 90.0), (0.5, 3.0)), verbose=True):
    """Busca la fase phi de cada ventana (±90°, paso 3° y luego 0,5°) con el resto fijo."""
    par = dict(par)
    par["phi"] = np.array(par["phi"], float)
    wg = _widx(par, np.array([g["t"].mean() for g in G["g"]]))
    for i in range(len(par["phi"])):
        sel = np.nonzero(wg == i)[0]
        if len(sel) < 3:
            continue
        Gs = _sub_grupos(G, sel)
        p0 = par["phi"][i]
        for paso, ancho in etapas:
            mejor = (9, p0)
            for dphi in np.radians(np.arange(-ancho, ancho + 1e-9, paso)):
                par["phi"][i] = p0 + dphi
                sc = float(residuos_grupos(Gs, tray, par, 0.5).mean())
                if sc < mejor[0]:
                    mejor = (sc, p0 + dphi)
            p0 = mejor[1]
            par["phi"][i] = p0
        if verbose:
            log("ventana %d: fase ajustada a %.1f deg (score %.3f, %d grupos)" % (i, np.degrees(p0), mejor[0], len(sel)))
    return par


def refinar_dW_nudos(G, tray, par, w_ini, w_fin, paso=5.0, pasadas=2, log=print):
    """Desplazamiento del mundo lineal a trozos (nudos cada 'paso' s): absorbe la diferencia lenta entre la trayectoria hf y la
    trayectoria optimizada (cierres de bucle) con la que se generó el LAS. Cada nudo se ajusta con los grupos vecinos."""
    kn = np.concatenate([np.arange(a, b + 1e-6, paso) for a, b in zip(w_ini, w_fin)])
    par = dict(par)
    par["kdW_t"], par["kdW"] = kn, np.zeros((len(kn), 3))
    gt = np.array([g["t"].mean() for g in G["g"]])
    for ps_ in range(pasadas):
        for i, tk in enumerate(kn):
            sel = np.nonzero(np.abs(gt - tk) < paso)[0]
            if len(sel) < 6:
                continue
            Gs = _sub_grupos(G, sel)

            def f(q, i=i):
                par["kdW"][i] = q
                return np.r_[residuos_grupos(Gs, tray, par, 0.5), 0.3 * q]

            r = least_squares(f, par["kdW"][i].copy(), loss="soft_l1", f_scale=0.1, max_nfev=8, diff_step=1e-2)
            par["kdW"][i] = r.x
        d = residuos_grupos(G, tray, par, 0.5)
        log("nudos dW pasada %d: mediana %.3f m, p90 %.3f m" % (ps_ + 1, np.median(d), np.percentile(d, 90)))
    return par


def _rot_min(a, b):
    """Rotación mínima que lleva el vector a al b."""
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    if s < 1e-12:
        return R.identity()
    return R.from_rotvec(v / s * np.arctan2(s, a @ b))


def calibrar(cfg, tray, tb, az, dist_raw, las_t, las_xyz, w_ini, w_fin, log=print, signos=None, semillas=None):
    """Devuelve (par, informe). tb/az/dist_raw: disparos del pcap de todas las ventanas;
    w_ini/w_fin: límites de las ventanas (s)."""
    ce = cfg["pcap"]["extrinseca"]
    rng = np.random.default_rng(cfg["proyecto"]["semilla_aleatoria"])
    dist_m = dist_raw.astype(np.float64) * cfg["pcap"]["unidad_distancia_m"]
    w_ini, w_fin = np.asarray(w_ini, float), np.asarray(w_fin, float)
    om_max = ce.get("omega_max_deg_s", 0)
    if om_max and om_max > 0:
        # Solo tramos sin giro brusco del cuerpo (vueltas en cabecera): ahí la interpolación de pose y el desfase temporal
        # dominan el error. Se calibra y se valida con los puntos del LAS de los tramos estables; P usa los mismos tramos.
        estable = tasa_giro(tray, las_t) <= om_max
        log("tramos estables (giro <= %.0f deg/s): %.1f %% de los puntos del LAS de la ventana" % (om_max, 100 * estable.mean()))
        las_t, las_xyz = las_t[estable], las_xyz[estable]
    centro = las_xyz.mean(0)
    mejor = None
    usados = set()                       # sellos k de TODOS los grupos que intervienen en el ajuste (se excluyen del hold-out)

    def _anotar(G_):
        if G_ is None:
            raise ValueError("calibración: no hay grupos temporales con correspondencias LAS-pcap en la ventana")
        usados.update(int(k) for k in G_["k"])
        return G_
    for signo in (signos or (cfg["pcap"]["signo_azimut"], -cfg["pcap"]["signo_azimut"])):
        g0 = estimar_giro(cfg, tb, az, dist_m, w_ini, w_fin, rng, signo, log)
        G = _anotar(grupos_temporales(las_t, las_xyz, tb, az, dist_m, 25, rng))
        up_b = np.mean([tray.pose(np.array([0.5 * (a + b)]))[1].inv().apply([0, 0, 1.0])[0] for a, b in zip(w_ini, w_fin)], axis=0)
        up_b /= np.linalg.norm(up_b)
        Er0 = _rot_min(g0["m_ref"], up_b)
        best = (9, None)
        for psi in np.radians(np.arange(0, 360, 2.0)):
            Er = R.from_rotvec(up_b * psi) * Er0
            par = dict(a_s=g0["a_s"], omega=g0["omega"], phi=g0["phi"], w_ini=w_ini, rotvec_e=Er.as_rotvec(), t_e=np.zeros(3),
                       dW=np.zeros(3), rotvec_w=np.zeros(3), k=1.0, centro=centro, signo=signo, az_off=0.0, el_off=0.0, dt=0.0)
            sc = float(residuos_grupos(G, tray, par, 0.5).mean())
            if sc < best[0]:
                best = (sc, par)
        log("signo_az=%+d: mejor guinada -> score %.3f (cono n.a=%.3f)" % (signo, best[0], g0["cono_n_a"]))
        if mejor is None or best[0] < mejor["sc"]:
            mejor = dict(sc=best[0], par=best[1], signo=signo, g0=g0)
    signo, par = mejor["signo"], mejor["par"]
    par["omega"] = np.full(len(w_ini), float(np.mean(par["omega"])))
    if semillas:
        # Mismo equipo y mismo montaje: se prueban como arranque las extrínsecas (Er, omega) de las calles ya calibradas y
        # válidas. Con la guiñada de la búsqueda ciega, las ventanas largas de calles nuevas caían en mínimos locales.
        Gs = _anotar(grupos_temporales(las_t, las_xyz, tb, az, dist_m, 40, rng))
        cands = [("busqueda_ciega", par)]
        for k, s_ in enumerate(semillas):
            cands.append(("semilla_%d" % k, dict(par, rotvec_e=np.asarray(s_["rotvec_e"], float),
                                                  omega=np.full(len(w_ini), float(s_["omega"])))))
        mej_c = (9, None, None)
        for nom_c, c in cands:
            c = escanear_fase_ventanas(Gs, tray, dict(c), log, etapas=((4.0, 180.0), (1.0, 4.0)), verbose=False)
            sc_c = float(residuos_grupos(Gs, tray, c, 0.5).mean())
            if sc_c < mej_c[0]:
                mej_c = (sc_c, nom_c, c)
        log("arranque elegido: %s (score %.3f; %d candidatos)" % (mej_c[1], mej_c[0], len(cands)))
        par = mej_c[2]
    Gp = _anotar(grupos_temporales(las_t, las_xyz, tb, az, dist_m, 90, rng))
    par = escanear_fase_ventanas(Gp, tray, par, log)
    G = _anotar(grupos_temporales(las_t, las_xyz, tb, az, dist_m, max(20, ce["n_puntos_ajuste"] // 300), rng))
    nw = len(w_ini)

    # vector: rotvec_e(3) t_e(3) k(1) omega(1, común) phi(nw) az_off el_off dt_ms(3) dW(3*nw)
    def empaquetar(p):
        return np.r_[p["rotvec_e"], p["t_e"], p["k"], p["omega"].mean(), p["phi"], p["az_off"], p["el_off"], p["dt"] * 1e3,
                     np.broadcast_to(p["dW"], (nw, 3)).ravel(),
                     np.broadcast_to(p.get("dV", np.zeros((nw, 3))), (nw, 3)).ravel(),
                     np.broadcast_to(p.get("dR", np.zeros((nw, 3))), (nw, 3)).ravel()]

    def desempaquetar(x):
        q = dict(par)
        q.update(rotvec_e=x[0:3], t_e=x[3:6], k=x[6], omega=np.full(nw, x[7]), phi=x[8:8 + nw],
                 az_off=x[8 + nw], el_off=x[9 + nw], dt=x[10 + nw] * 1e-3, dW=x[11 + nw:11 + nw + 3 * nw].reshape(nw, 3),
                 dV=x[11 + 4 * nw:11 + 7 * nw].reshape(nw, 3), dR=x[11 + 7 * nw:11 + 10 * nw].reshape(nw, 3))
        return q

    x = empaquetar(par)
    libres = np.ones(len(x), bool)
    libres[3:6] = False                                      # palanca fija en 0: queda absorbida por dW de cada ventana (rumbo casi constante en cada paso)
    libres[6] = bool(ce["ajustar_escala_distancia"])
    libres[9 + nw] = bool(ce["ajustar_offset_elevacion"])   # el_off
    libres[10 + nw] = True                                     # dt (ms)
    libres[11 + nw:] = bool(ce["ajustar_traslacion_mundo"])
    libres[11 + 4 * nw:11 + 7 * nw] = bool(ce.get("ajustar_deriva_dW", True))
    libres[11 + 7 * nw:] = bool(ce.get("ajustar_deriva_rotacion", True))
    for tope in ce["dist_max_correspondencia_m"]:
        def fun(q, tope=tope):
            xx = x.copy()
            xx[libres] = q
            pp = desempaquetar(xx)
            return np.r_[residuos_grupos(G, tray, pp, tope), 0.05 * pp["t_e"], 0.02 * pp["dW"].ravel(), 0.5 * pp["dV"].ravel(), 5.0 * pp["dR"].ravel()]   # a priori débil: palanca y desplazamiento pequeños
        r = least_squares(fun, x[libres], loss="soft_l1", f_scale=tope / 3, max_nfev=ce["iteraciones_icp"], diff_step=1e-4)
        x[libres] = r.x
        log("tope %.2f m -> coste %.4f" % (tope, r.cost))
    par = desempaquetar(x)
    par["signo"] = signo
    log("tras ajuste global: mediana %.3f m" % np.median(residuos_grupos(G, tray, par, 2.0)))
    par = refinar_fase(G, tray, par, w_ini, w_fin, log=log)
    if ce.get("nudos_dW_s", 0):
        par = refinar_dW_nudos(G, tray, par, w_ini, w_fin, paso=ce["nudos_dW_s"], log=log)
    d = residuos_grupos(G, tray, par, 2.0)
    dg = residuos_grupos(G, tray, par, 2.0, por_grupo=True)
    wg = _widx(par, np.array([g["t"].mean() for g in G["g"]]))
    por_ventana = [dict(ventana=int(i), n=int(sum(len(dg[j]) for j in np.nonzero(wg == i)[0])),
                        mediana_m=float(np.median(np.concatenate([dg[j] for j in np.nonzero(wg == i)[0]]))) if (wg == i).any() else None)
                   for i in range(nw)]
    informe = dict(residuo_por_ventana=por_ventana, signo_azimut=int(signo), a_s=par["a_s"].tolist(), omega=par["omega"].tolist(), phi=par["phi"].tolist(),
                   w_ini=w_ini.tolist(), rotvec_e_deg=np.degrees(par["rotvec_e"]).tolist(), t_e_m=par["t_e"].tolist(),
                   dW_m=np.asarray(par["dW"]).tolist(), dV_m_s=np.asarray(par["dV"]).tolist(), kdW_t=np.asarray(par.get("kdW_t", [])).tolist(), kdW=np.asarray(par.get("kdW", [])).tolist(), dR_rad_s=np.asarray(par["dR"]).tolist(), rotvec_w_deg=np.degrees(par["rotvec_w"]).tolist(), k_escala=float(par["k"]),
                   az_off_deg=float(par["az_off"]), el_off_deg=float(par["el_off"]), dt_ms=float(par["dt"] * 1e3), n_grupos=len(G["g"]), n_fuentes=int(len(d)), cono_n_a=mejor["g0"]["cono_n_a"],
                   residuo_mediana_m=float(np.median(d)), residuo_p90_m=float(np.percentile(d, 90)),
                   residuo_rms_inliers_m=float(np.sqrt(np.mean(d[d < 0.3] ** 2))) if (d < 0.3).any() else None,
                   frac_menor_0_05=float(np.mean(d < 0.05)), frac_menor_0_15=float(np.mean(d < 0.15)),
                   score_busqueda_psi=mejor["sc"], centro=centro.tolist(),
                   kn_t=par["kn_t"].tolist(), kn_d=par["kn_d"].tolist(),
                   k_ajuste=sorted(usados))                 # sellos k usados en el ajuste: el hold-out no puede usarlos
    return par, informe


def par_desde_informe(inf):
    return dict(a_s=np.array(inf["a_s"]), omega=np.array(inf["omega"]), phi=np.array(inf["phi"]), w_ini=np.array(inf["w_ini"]),
                rotvec_e=np.radians(inf["rotvec_e_deg"]), t_e=np.array(inf["t_e_m"]), dW=np.array(inf["dW_m"]), kdW_t=np.array(inf.get("kdW_t", [])), kdW=np.array(inf.get("kdW", np.zeros((0, 3)))), dV=np.array(inf.get("dV_m_s", np.zeros((len(inf["w_ini"]), 3)))), dR=np.array(inf.get("dR_rad_s", np.zeros((len(inf["w_ini"]), 3)))),
                rotvec_w=np.radians(inf["rotvec_w_deg"]), k=inf["k_escala"], centro=np.array(inf["centro"]), signo=inf["signo_azimut"],
                kn_t=np.array(inf.get("kn_t", [])), kn_d=np.array(inf.get("kn_d", [])), valida=bool(inf.get("valida", True)),
                az_off=inf.get("az_off_deg", 0.0), el_off=inf.get("el_off_deg", 0.0), dt=inf.get("dt_ms", 0.0) * 1e-3)


def calibrar_por_ventana(cfg, tray, ventanas, lims, las_t, las_xyz, log=print, semillas=None):
    """Una extrínseca por ventana (paso del sensor por una calle). Con una sola extrínseca global la ventana
    del paso en sentido contrario no ajustaba (mediana ~1 m), así que cada paso lleva su propia Er, fase y dW.
    Devuelve (lista de par, lista de informes)."""
    pars, infs = [], []
    signos = None if cfg["pcap"]["extrinseca"].get("probar_signo_azimut", False) else (cfg["pcap"]["signo_azimut"],)
    for i, (w, (a, b)) in enumerate(zip(ventanas, lims)):
        m = (las_t >= a) & (las_t <= b)
        log("--- calibración de la ventana %d [%.1f, %.1f]" % (i, a, b))
        par, inf = calibrar(cfg, tray, w["t"], w["az"], w["dist_raw"], las_t[m], las_xyz[m], [a], [b], log=log, signos=signos, semillas=semillas)
        inf["ventana_t"] = [a, b]
        ce = cfg["pcap"]["extrinseca"]
        inf["valida"] = bool(inf["residuo_mediana_m"] < ce["residuo_mediana_max_m"] and inf["residuo_p90_m"] < ce["residuo_p90_max_m"])
        par["valida"] = inf["valida"]
        log("ventana %d: residuo mediana %.3f m (max %.2f), p90 %.3f m (max %.2f) -> %s" % (
            i, inf["residuo_mediana_m"], ce["residuo_mediana_max_m"], inf["residuo_p90_m"], ce["residuo_p90_max_m"],
            "VALIDA" if inf["valida"] else "NO VALIDA (se excluye de P)"))
        pars.append(par)
        infs.append(inf)
    return pars, infs
