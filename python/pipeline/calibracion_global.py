# -*- coding: utf-8 -*-
"""Calibración GLOBAL del escáner por bloque del ensayo (opción pcap.modo_calibracion = "global_bloque").

Por qué: con la lectura original de la orientación de la trayectoria (columnas 5-8 del high_freq como qx qy qz qw)
cada calle necesitaba su propio montaje y muchas no pasaban la validación. Con otra lectura de esas columnas
(pcap.calibracion_global.lectura_cuaternion) y la tabla de elevaciones del fichero LAS-1.EP del RS30, un único
montaje del escáner vale para todas las calles del bloque (experimento del 2026-11-17/18 en IFAPA Cabra:
39 de 45 calles válidas frente a 12 de 45 con la calibración por calle).

Modelo: el de extrinseca.py (a_mundo, residuos_grupos) con una ventana por calle.
  Común a todas las calles : giro escáner->cuerpo (Er), eje del motor a_s, palanca t_e, escala k, az_off, el_off, dt.
  Propio de cada calle      : omega y phi (estado del motor), nudos de fase, dW, dV, dR y nudos de dW
                              (diferencia lenta entre la trayectoria high_freq y la que usó el fabricante para el LAS).
Arranque: la calibración por calle (paso 09, fase "calibrar") de la calle con MENOR p90 de validación del bloque
(si hay válidas, entre las válidas); su montaje se pasa a la nueva lectura y la fase de cada calle se alinea con él
y se vuelve a buscar (±180°). Ajuste conjunto por mínimos cuadrados (misma pérdida soft_l1, topes y pesos a priori que
extrinseca.calibrar) con jacobiano disperso. Validación: grupos de hold-out de cada calle (los mismos que
validacion.holdout: semilla 777, 600 grupos, tramos estables) y el MISMO umbral (mediana < 0,15 m y p90 < 0,30 m).

Entradas: <salida>/09_metricas/_calles/calle_*.json (calibración por calle) y calle_*_grupos.pkl (grupos de ajuste y
de validación, los guarda 09_metricas en la fase "calibrar"). Salida: <salida>/09_metricas/_global/montaje_bloque.json
(con huella) y log en logs/09b_calibracion_global.log.
"""
import utils  # noqa: F401  (primero: entorno)
import glob
import os
import pickle
import struct
import time

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.spatial.transform import Rotation as R

from utils import Paso, dir_salida, guardar_json, cargar_json, huella, huella_codigo
from trayectoria import cargar_trayectoria, cargar_trayectoria_lectura, nombre_lectura
import extrinseca as ex
import pcap_hesai as ph

NOMBRE = "09b_calibracion_global"
VERSION = 1


# --------------------------------------------------------------------------------------
# Tabla de elevaciones del fichero LAS-1.EP (Protocol Buffers)
# --------------------------------------------------------------------------------------
def leer_elevaciones_ep(ruta_ep):
    """Elevación de los 32 canales (grados) = campo 7 de la tabla 7 del registro 0 de LAS-1.EP (en radianes).
    El significado del campo no lo documenta el fabricante: se comprobó que mejora el residuo frente al LAS."""
    def _varint(b, i):
        r = s = 0
        while True:
            c = b[i]; i += 1; r |= (c & 0x7f) << s; s += 7
            if not c & 0x80:
                return r, i

    def _parse(b):
        i, out = 0, []
        while i < len(b):
            k, i = _varint(b, i); f, w = k >> 3, k & 7
            if w == 0: v, i = _varint(b, i)
            elif w == 1: v = struct.unpack('<d', b[i:i + 8])[0]; i += 8
            elif w == 5: v = struct.unpack('<f', b[i:i + 4])[0]; i += 4
            else: n, i = _varint(b, i); v = b[i:i + n]; i += n
            out.append((f, w, v))
        return out

    with open(ruta_ep, "rb") as f:
        b = f.read()
    t7 = _parse([v for f, w, v in _parse(_parse(b)[0][2]) if f == 7][0])
    el = np.degrees([v for f, w, v in t7 if f == 7])
    if len(el) != 32:
        raise ValueError("LAS-1.EP: se esperaban 32 elevaciones y hay %d" % len(el))
    return el


def tabla_elevaciones(cfg):
    cg = cfg["pcap"]["calibracion_global"]
    if cg["elevaciones"] == "EP":
        ruta = cfg["rutas"].get("calib_ep")
        if not ruta or not os.path.exists(ruta):
            raise FileNotFoundError("pcap.calibracion_global.elevaciones = EP pero no se encuentra LAS-1.EP en la carpeta de brutos "
                                    "(SCAN/Scanner1/*.EP); indícalo en entradas.brutos_archivos.calib_ep")
        return leer_elevaciones_ep(ruta)
    if cg["elevaciones"] == "xt32_estandar":
        return ph.ELEV_XT32_DEG.copy()
    raise ValueError("pcap.calibracion_global.elevaciones debe ser 'EP' o 'xt32_estandar'")


def trayectoria_global(cfg):
    return cargar_trayectoria_lectura(cfg, cfg["pcap"]["calibracion_global"]["lectura_cuaternion"],
                                      dir_salida(cfg, "01_qc"))[0]


# --------------------------------------------------------------------------------------
# Grupos por calle (los escribe 09_metricas en la fase "calibrar")
# --------------------------------------------------------------------------------------
def ruta_grupos(ruta_calle):
    return ruta_calle.replace(".json", "_grupos.pkl")


def guardar_grupos(cfg, ruta_calle, inf, w, las_t, las_xyz, hilera):
    """Grupos temporales de ajuste y de validación de una calle (tramos estables), idénticos a los de calibrar
    (semilla del proyecto, n_puntos_ajuste // 300 grupos) y a los de validacion.holdout (semilla 777, 600 grupos)."""
    cg = cfg["pcap"]["calibracion_global"]
    ce = cfg["pcap"]["extrinseca"]
    om_max = ce.get("omega_max_deg_s", 0)
    tray = cargar_trayectoria(cfg, cache_dir=dir_salida(cfg, "01_qc"))[0]
    if om_max:
        est = ex.tasa_giro(tray, las_t) <= om_max
        las_t, las_xyz = las_t[est], las_xyz[est]
    dist_m = w["dist_raw"].astype(np.float64) * cfg["pcap"]["unidad_distancia_m"]
    Gh = ex.grupos_temporales(las_t, las_xyz, w["t"], w["az"], dist_m, cg["n_grupos_holdout"],
                              np.random.default_rng(cg["semilla_holdout"]))
    Gf = ex.grupos_temporales(las_t, las_xyz, w["t"], w["az"], dist_m, max(20, ce["n_puntos_ajuste"] // 300),
                              np.random.default_rng(cfg["proyecto"]["semilla_aleatoria"]))
    if Gh is None or Gf is None:
        return False
    tmp = ruta_grupos(ruta_calle) + ".tmp"
    with open(tmp, "wb") as fo:
        pickle.dump(dict(calle=os.path.basename(ruta_calle), huella_calle=inf["huella"], hilera=str(hilera),
                         ventana=[float(x) for x in inf["ventana_global"]], Gf=Gf, Gh=Gh), fo, protocol=4)
    os.replace(tmp, ruta_grupos(ruta_calle))
    return True


def grupos_validos(ruta_calle, inf):
    r = ruta_grupos(ruta_calle)
    if not os.path.exists(r):
        return False
    try:
        with open(r, "rb") as f:
            return pickle.load(f)["huella_calle"] == inf.get("huella")
    except Exception:
        return False


# --------------------------------------------------------------------------------------
# Ajuste conjunto
# --------------------------------------------------------------------------------------
def _concat(gs):
    cat = lambda k: np.concatenate([g[k] for g in gs])
    return dict(g=gs, t=cat("t"), az=cat("az"), ch=cat("ch"), r=cat("r"), lim=np.cumsum([0] + [len(g["t"]) for g in gs]))


def candidatos_arranque(D, n):
    """Calles de arranque: las n de menor p90 de validación calle a calle (entre las válidas si las hay; si faltan, se completan
    con las no válidas de menor p90)."""
    orden = sorted(D, key=lambda d: (not d["oficial"].get("valida"), d["oficial"]["holdout_p90_m"]))
    return [d["calle"] for d in orden[:max(1, n)]]


def ajustar(cfg, D, tray_of, tray, elev, log, calle_ref=None):
    """D: lista de dict(calle, ventana, Gf, Gh, oficial) ordenada por tiempo. Devuelve (par, lista de resultados por calle)."""
    cg = cfg["pcap"]["calibracion_global"]
    ce = cfg["pcap"]["extrinseca"]
    nw = len(D)
    w_ini = np.array([d["ventana"][0] for d in D]); w_fin = np.array([d["ventana"][1] for d in D])
    NL = cg["n_grupos_lsq_por_calle"]
    G_full = _concat([g for d in D for g in d["Gf"]["g"]])                  # nudos de fase y de dW
    G = _concat([g for d in D for g in d["Gf"]["g"][:NL]])                  # mínimos cuadrados conjuntos
    n_res = [sum(len(g["P"]) for g in d["Gf"]["g"][:NL]) for d in D]
    log("calles: %d | grupos: %d (nudos) / %d (mínimos cuadrados); filas: %d" % (nw, len(G_full["g"]), len(G["g"]), sum(n_res)))

    def er_en_lectura(inf):
        """Er de una calle (calibración por calle, lectura original) expresado con la nueva lectura (media en la ventana)."""
        p = ex.par_desde_informe(inf)
        a, b = inf["ventana_global"]
        ts = np.linspace(a + 3, b - 3, 25)
        M = ex._Rw(p, ts) * tray_of.pose(ts)[1] * R.from_rotvec(p["rotvec_e"])
        return (tray.pose(ts)[1].inv() * M).mean()

    # ---- arranque: calle de referencia = menor p90 de validación (entre las válidas si las hay)
    ref = [d for d in D if d["calle"] == (calle_ref or candidatos_arranque(D, 1)[0])][0]
    inf_r = ref["oficial"]
    log("calle de referencia (arranque): %s (por calle: mediana %.3f, p90 %.3f, %s)" % (
        ref["calle"], inf_r["holdout_mediana_m"], inf_r["holdout_p90_m"], "válida" if inf_r.get("valida") else "no válida"))
    Er_g = er_en_lectura(inf_r)
    a_s = np.array(inf_r["a_s"], float); a_s /= np.linalg.norm(a_s)
    phi, omega, dW, dts = [], [], [], []
    for d in D:
        inf = d["oficial"]
        q = (Er_g.inv() * er_en_lectura(inf)).as_quat()
        delta = 2 * np.arctan2(q[:3] @ a_s, q[3])                         # giro alrededor de a_s (swing-twist)
        phi.append(inf["phi"][0] - delta); omega.append(inf["omega"][0])
        dW.append(np.asarray(inf["dW_m"], float).reshape(-1, 3)[0]); dts.append(inf["dt_ms"])
    par = dict(a_s=a_s, omega=np.array(omega), phi=np.array(phi), w_ini=w_ini, rotvec_e=Er_g.as_rotvec(), t_e=np.zeros(3),
               dW=np.array(dW), dV=np.zeros((nw, 3)), dR=np.zeros((nw, 3)), rotvec_w=np.zeros(3), k=1.0,
               centro=np.concatenate([g["P"] for g in G["g"]]).mean(0), signo=inf_r["signo_azimut"], az_off=0.0,
               el_off=0.0 if cg["elevaciones"] == "EP" else float(np.median([d["oficial"]["el_off_deg"] for d in D])),
               dt=float(np.median(dts)) * 1e-3, elev_deg=elev)
    d0 = ex.residuos_grupos(G, tray, par, 0.5)
    log("arranque: mediana %.3f p90 %.3f" % (np.median(d0), np.percentile(d0, 90)))
    par = ex.escanear_fase_ventanas(G, tray, par, log, etapas=((4.0, 180.0), (1.0, 4.0)), verbose=False)
    d0 = ex.residuos_grupos(G, tray, par, 0.5)
    log("tras buscar la fase de cada calle: mediana %.3f p90 %.3f" % (np.median(d0), np.percentile(d0, 90)))

    # ---- ajuste conjunto (jacobiano disperso: cada calle solo depende de lo común y de lo suyo)
    NC, NP = 10, 11          # común: rotvec_e(3) t_e(3) k az_off el_off dt_ms | por calle: omega phi dW(3) dV(3) dR(3)

    def empaquetar(p):
        x = [p["rotvec_e"], p["t_e"], [p["k"], p["az_off"], p["el_off"], p["dt"] * 1e3]]
        for i in range(nw):
            x += [[p["omega"][i], p["phi"][i]], p["dW"][i], p["dV"][i], p["dR"][i]]
        return np.concatenate([np.ravel(v) for v in x]).astype(float)

    def desempaquetar(x):
        q = dict(par)
        B = x[NC:].reshape(nw, NP)
        q.update(rotvec_e=x[0:3], t_e=x[3:6], k=x[6], az_off=x[7], el_off=x[8], dt=x[9] * 1e-3,
                 omega=B[:, 0].copy(), phi=B[:, 1].copy(), dW=B[:, 2:5].copy(), dV=B[:, 5:8].copy(), dR=B[:, 8:11].copy())
        return q

    x = empaquetar(par)
    libres = np.ones(len(x), bool)
    if cg["palanca"] == "cero":
        libres[3:6] = False
    libres[6] = bool(ce["ajustar_escala_distancia"])
    libres[8] = bool(ce["ajustar_offset_elevacion"])
    B = libres[NC:].reshape(nw, NP)
    B[:, 2:5] = bool(ce["ajustar_traslacion_mundo"]); B[:, 5:8] = bool(ce.get("ajustar_deriva_dW", True))
    B[:, 8:11] = bool(ce.get("ajustar_deriva_rotacion", True))
    libres[NC:] = B.ravel()
    n_fil = sum(n_res) + 3 + 9 * nw                    # residuo por calle + a priori (t_e 3, dW, dV, dR 3*nw cada uno)
    S = lil_matrix((n_fil, len(x)), dtype=np.int8)
    r0 = 0
    for i, n in enumerate(n_res):
        S[r0:r0 + n, :NC] = 1
        S[r0:r0 + n, NC + i * NP:NC + (i + 1) * NP] = 1
        r0 += n
    for j in range(3):
        S[r0 + j, 3 + j] = 1
    r0 += 3
    for bl in (2, 5, 8):                               # dW, dV, dR
        for i in range(nw):
            for j in range(3):
                S[r0 + i * 3 + j, NC + i * NP + bl + j] = 1
        r0 += 3 * nw
    S = S.tocsc()[:, libres]
    for tope in ce["dist_max_correspondencia_m"]:
        def fun(qv, tope=tope):
            xx = x.copy(); xx[libres] = qv
            pp = desempaquetar(xx)
            return np.r_[ex.residuos_grupos(G, tray, pp, tope), 0.05 * pp["t_e"], 0.02 * pp["dW"].ravel(),
                         0.5 * pp["dV"].ravel(), 5.0 * pp["dR"].ravel()]
        r = least_squares(fun, x[libres], loss="soft_l1", f_scale=tope / 3, max_nfev=ce["iteraciones_icp"],
                          diff_step=1e-4, jac_sparsity=S)
        x[libres] = r.x
        dd = ex.residuos_grupos(G, tray, desempaquetar(x), 0.5)
        log("tope %.2f m -> coste %.1f | mediana %.3f p90 %.3f" % (tope, r.cost, np.median(dd), np.percentile(dd, 90)))
    par = desempaquetar(x)
    par = ex.refinar_fase(G_full, tray, par, w_ini, w_fin, log=log)
    if ce.get("nudos_dW_s", 0):
        par = ex.refinar_dW_nudos(G_full, tray, par, w_ini, w_fin, paso=ce["nudos_dW_s"], log=log)

    # ---- validación por calle (hold-out, mismo umbral)
    res = []
    for i, d in enumerate(D):
        hd = ex.residuos_grupos(d["Gh"], tray, par, 2.0)
        fd = ex.residuos_grupos(d["Gf"], tray, par, 2.0)
        med, p90 = float(np.median(hd)), float(np.percentile(hd, 90))
        ok = bool(med < ce["residuo_mediana_max_m"] and p90 < ce["residuo_p90_max_m"])
        o = d["oficial"]
        res.append(dict(calle=d["calle"], ventana=[float(w_ini[i]), float(w_fin[i])], holdout_mediana_m=med, holdout_p90_m=p90,
                        holdout_n=int(len(hd)), frac_menor_0_15=float(np.mean(hd < 0.15)),
                        residuo_rms_inliers_m=float(np.sqrt(np.mean(hd[hd < 0.3] ** 2))) if (hd < 0.3).any() else None,
                        valida=ok, ajuste_mediana_m=float(np.median(fd)), ajuste_p90_m=float(np.percentile(fd, 90)),
                        por_calle_mediana_m=o["holdout_mediana_m"], por_calle_p90_m=o["holdout_p90_m"],
                        por_calle_valida=bool(o["valida"]), phi_deg=float(np.degrees(par["phi"][i])),
                        periodo_s=float(2 * np.pi / abs(par["omega"][i]))))
        log("%s: hold-out %.3f / %.3f -> %s   (por calle %.3f / %.3f %s)" % (
            d["calle"], med, p90, "VALIDA" if ok else "no", o["holdout_mediana_m"], o["holdout_p90_m"],
            "VALIDA" if o["valida"] else "no"))
    return par, res, ref["calle"]


# --------------------------------------------------------------------------------------
# Montaje -> extrínseca de una ventana (para P)
# --------------------------------------------------------------------------------------
def par_a_json(par):
    return dict({k: np.asarray(par[k]).tolist() for k in ("a_s", "omega", "phi", "w_ini", "rotvec_e", "t_e", "dW", "dV", "dR",
                                                          "rotvec_w", "centro", "kn_t", "kn_d", "elev_deg")},
                kdW_t=np.asarray(par.get("kdW_t", [])).tolist(), kdW=np.asarray(par.get("kdW", np.zeros((0, 3)))).tolist(),
                k=float(par["k"]), signo=int(par["signo"]), az_off=float(par["az_off"]), el_off=float(par["el_off"]),
                dt=float(par["dt"]))


def par_ventana(mont, i_calle, valida):
    """Extrínseca de UNA calle (formato de extrinseca.par_desde_informe) a partir del montaje del bloque: lo común +
    lo propio de la calle i (sus nudos de fase y de dW, recortados a su ventana)."""
    j = mont["par"]
    a, b = mont["calles"][i_calle]["ventana"]
    kn_t, kn_d = np.array(j["kn_t"]), np.array(j["kn_d"])
    mk = (kn_t >= a - 1e-6) & (kn_t <= b + 1e-6)
    kdt, kd = np.array(j["kdW_t"]), np.array(j["kdW"]).reshape(-1, 3)
    md = (kdt >= a - 1e-6) & (kdt <= b + 1e-6) if len(kdt) else np.zeros(0, bool)
    i = i_calle
    return dict(a_s=np.array(j["a_s"]), omega=np.array([j["omega"][i]]), phi=np.array([j["phi"][i]]), w_ini=np.array([a]),
                rotvec_e=np.array(j["rotvec_e"]), t_e=np.array(j["t_e"]), dW=np.array([j["dW"][i]]), dV=np.array([j["dV"][i]]),
                dR=np.array([j["dR"][i]]), rotvec_w=np.array(j["rotvec_w"]), k=j["k"], centro=np.array(j["centro"]),
                signo=j["signo"], kn_t=kn_t[mk], kn_d=kn_d[mk], kdW_t=kdt[md] if len(kdt) else kdt,
                kdW=kd[md] if len(kdt) else np.zeros((0, 3)), az_off=j["az_off"], el_off=j["el_off"], dt=j["dt"],
                elev_deg=np.array(j["elev_deg"]), valida=bool(valida))


def ruta_montaje(cfg):
    d = os.path.join(cfg["rutas"]["salida"], "09_metricas", "_global")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, "montaje_bloque.json")


def cargar_montaje(cfg):
    r = ruta_montaje(cfg)
    if not os.path.exists(r):
        raise FileNotFoundError("Falta el montaje global del bloque (%s): ejecuta el pipeline con pcap.modo_calibracion = "
                                "global_bloque desde el paso 09 (run_pipeline.py lo calcula entre las dos fases del 09)" % r)
    return cargar_json(r)


def main(cfg):
    """Ajuste global del bloque con todas las calles calibradas en _calles. Caché con huella."""
    cg = cfg["pcap"]["calibracion_global"]
    calles_dir = os.path.join(cfg["rutas"]["salida"], "09_metricas", "_calles")
    rutas = sorted(f for f in glob.glob(os.path.join(calles_dir, "calle_*.json")))
    with Paso(cfg, NOMBRE, entrada=calles_dir, parametros=dict(calibracion_global=cg, extrinseca=cfg["pcap"]["extrinseca"])) as p:
        D = []
        for f in rutas:
            inf = cargar_json(f)
            if not grupos_validos(f, inf):
                raise RuntimeError("La calle %s no tiene grupos de calibración válidos (%s): repite el paso 09 en fase 'calibrar'"
                                   % (os.path.basename(f), os.path.basename(ruta_grupos(f))))
            with open(ruta_grupos(f), "rb") as fo:
                g = pickle.load(fo)
            D.append(dict(calle=os.path.basename(f), ventana=g["ventana"], Gf=g["Gf"], Gh=g["Gh"], oficial=inf))
        if not D:
            raise RuntimeError("No hay calles calibradas en %s" % calles_dir)
        D.sort(key=lambda d: d["ventana"][0])
        h = huella(dict(v=VERSION, calles=[(d["calle"], d["oficial"]["huella"]) for d in D], cg=cg, ce=cfg["pcap"]["extrinseca"],
                        codigo=huella_codigo("calibracion_global.py", "extrinseca.py", "trayectoria.py")))
        ruta = ruta_montaje(cfg)
        if os.path.exists(ruta) and cfg["pcap"]["extrinseca"].get("reutilizar", True):
            try:
                m = cargar_json(ruta)
                if m.get("huella") == h:
                    p.info("montaje global reutilizado (huella %s): %d/%d calles válidas", h, m["n_validas"], len(m["calles"]))
                    p.salida = ruta
                    return m
            except Exception:
                pass
        t0 = time.time()
        elev = tabla_elevaciones(cfg)
        p.info("lectura del cuaternión %s | elevaciones %s (canales 1, 16, 32: %s) | palanca %s", nombre_lectura(cg["lectura_cuaternion"]),
               cg["elevaciones"], np.round(elev[[0, 15, 31]], 3).tolist(), cg["palanca"])
        tray_of = cargar_trayectoria(cfg, cache_dir=dir_salida(cfg, "01_qc"))[0]
        tray = trayectoria_global(cfg)
        # Varios arranques (pcap.calibracion_global.n_arranques): se ajusta desde cada calle candidata y se queda el montaje con
        # más calles válidas en hold-out (a igualdad, menor mediana global del hold-out).
        arranques, mejor = [], None
        for cref in candidatos_arranque(D, int(cg.get("n_arranques", 1))):
            p.info("=== arranque desde %s", cref)
            par_i, res_i, ref_i = ajustar(cfg, D, tray_of, tray, elev, p.info, calle_ref=cref)
            nv = int(sum(c["valida"] for c in res_i)); med = float(np.median([c["holdout_mediana_m"] for c in res_i]))
            arranques.append(dict(calle_referencia=ref_i, n_validas=nv, holdout_mediana_global_m=med))
            p.info("arranque %s: %d/%d calles válidas, mediana global %.3f m", ref_i, nv, len(res_i), med)
            if mejor is None or (nv, -med) > (mejor[3], -mejor[4]):
                mejor = (par_i, res_i, ref_i, nv, med)
        par, res, ref = mejor[:3]
        p.info("montaje elegido: arranque desde %s", ref)
        Er = R.from_rotvec(par["rotvec_e"])
        m = dict(huella=h, lectura=dict(cg["lectura_cuaternion"], nombre=nombre_lectura(cg["lectura_cuaternion"])),
                 elevaciones=cg["elevaciones"], palanca=cg["palanca"], calle_referencia=ref, arranques=arranques,
                 n_validas=int(sum(c["valida"] for c in res)), n_validas_por_calle=int(sum(c["por_calle_valida"] for c in res)),
                 holdout_mediana_global_m=float(np.median([c["holdout_mediana_m"] for c in res])),
                 comun=dict(rotvec_e_deg=np.degrees(par["rotvec_e"]).tolist(), eje_motor_cuerpo=Er.apply(par["a_s"]).tolist(),
                            t_e_m=np.asarray(par["t_e"]).tolist(), k=float(par["k"]), az_off_deg=float(par["az_off"]),
                            el_off_deg=float(par["el_off"]), dt_ms=float(par["dt"] * 1e3)),
                 calles=res, par=par_a_json(par), tiempo_s=round(time.time() - t0))
        guardar_json(m, ruta)
        p.info("FIN: %d/%d calles válidas con el montaje global (por calle: %d/%d); palanca %s m; dt %.1f ms; el_off %.2f°",
               m["n_validas"], len(res), m["n_validas_por_calle"], len(res), np.round(par["t_e"], 3).tolist(),
               par["dt"] * 1e3, par["el_off"])
        p.extra["calles_validas"] = "%d/%d" % (m["n_validas"], len(res))
        p.salida = ruta
    return m
