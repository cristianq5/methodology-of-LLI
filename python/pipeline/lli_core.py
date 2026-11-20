# -*- coding: utf-8 -*-
"""Núcleo numérico del LLI (sin E/S): MH, MW, envolvente de copa, clasificación de rayos
(Bt / Bi), P, L y LLI. Lo usan 09_metricas, 10_lli y tests\\test_lli_sintetico.py.

Marco de hilera: s (a lo largo del eje), t (perpendicular, + izquierda), h = altura sobre el suelo.
  LLI = MH * MW * (1 - P),  P = Bi / Bt  con  Bi = rayos que ATRAVIESAN la copa  (P = porosidad, L = 1 - P).
"""
import numpy as np


# --------------------------------------------------------------------------------------
# MH y MW por rebanada, y envolvente por celda (rebanada x franja)
# --------------------------------------------------------------------------------------
class Perfil:
    """Resultado por rebanada + envolvente por celda."""
    pass


def perfiles_rebanada(s, t, h, s_ini, n_reb, rebanada=0.1, franja=0.1, umbral=0.3,
                      min_pts_franja=5, min_pts_reb=50, mh_percentil=99.0, mw_percentil_t=(2.5, 97.5),
                      n_franjas=None, env_percentil_t=(1.0, 99.0)):
    """MH, MW (estrictos y con percentil) y envolvente [tmin, tmax] por celda.

    Sólo cuentan los puntos con h >= umbral (0,3 m sobre el suelo).
      MH_max  = max(h) de la rebanada;  MH_p = percentil 'mh_percentil' de h.
      Franja j = [umbral + j*franja, umbral + (j+1)*franja). Anchura de la franja:
        estricta = tmax - tmin;  robusta = p(97,5) - p(2,5) de t.
      MW = media de las anchuras de las franjas con >= min_pts_franja puntos.
    Rebanadas con < min_pts_reb puntos -> NaN.
    """
    ok = h >= umbral
    s, t, h = s[ok], t[ok], h[ok]
    i = np.floor((s - s_ini) / rebanada).astype(np.int64)
    j = np.floor((h - umbral) / franja + 1e-9).astype(np.int64)
    v = (i >= 0) & (i < n_reb) & (j >= 0)
    i, j, s, t, h = i[v], j[v], s[v], t[v], h[v]
    nf = int(j.max()) + 1 if n_franjas is None and len(j) else (n_franjas or 1)
    j = np.minimum(j, nf - 1)
    P = Perfil()
    P.n_reb, P.nf = n_reb, nf
    P.rebanada, P.franja, P.umbral = rebanada, franja, umbral
    P.s_ini = s_ini
    P.n_pts = np.bincount(i, minlength=n_reb)
    # por celda: n, tmin, tmax y percentiles de t
    cel = i * nf + j
    n_cel = np.bincount(cel, minlength=n_reb * nf)
    tmin = np.full(n_reb * nf, np.nan)
    tmax = np.full(n_reb * nf, np.nan)
    tlo = np.full(n_reb * nf, np.nan)
    thi = np.full(n_reb * nf, np.nan)
    elo = np.full(n_reb * nf, np.nan)
    ehi = np.full(n_reb * nf, np.nan)
    if len(cel):
        o = np.argsort(cel, kind="stable")
        cs, ts = cel[o], t[o]
        ini = np.r_[0, np.nonzero(np.diff(cs))[0] + 1]
        fin = np.r_[ini[1:], len(cs)]
        ids = cs[ini]
        tmin[ids] = np.minimum.reduceat(ts, ini)
        tmax[ids] = np.maximum.reduceat(ts, ini)
        for a, b, c in zip(ini, fin, ids):
            if b - a >= min_pts_franja:
                seg = ts[a:b]
                tlo[c], thi[c] = np.percentile(seg, mw_percentil_t)
                elo[c], ehi[c] = np.percentile(seg, env_percentil_t)   # envolvente robusta de P (QA-14)
    P.n_cel = n_cel.reshape(n_reb, nf)
    P.tmin, P.tmax = tmin.reshape(n_reb, nf), tmax.reshape(n_reb, nf)
    P.tlo, P.thi = tlo.reshape(n_reb, nf), thi.reshape(n_reb, nf)
    P.elo, P.ehi = elo.reshape(n_reb, nf), ehi.reshape(n_reb, nf)
    P.ocupada = P.n_cel >= min_pts_franja
    ancho_e = np.where(P.ocupada, P.tmax - P.tmin, np.nan)
    ancho_p = np.where(P.ocupada, P.thi - P.tlo, np.nan)
    P.MW_max = np.nanmean(np.where(np.isfinite(ancho_e), ancho_e, np.nan), axis=1) if nf else np.full(n_reb, np.nan)
    P.MW_p = np.nanmean(ancho_p, axis=1)
    mh_max = np.full(n_reb, np.nan)
    mh_p = np.full(n_reb, np.nan)
    if len(i):
        o = np.argsort(i, kind="stable")
        is_, hs = i[o], h[o]
        ini = np.r_[0, np.nonzero(np.diff(is_))[0] + 1]
        fin = np.r_[ini[1:], len(is_)]
        for a, b in zip(ini, fin):
            mh_max[is_[a]] = hs[a:b].max()
            mh_p[is_[a]] = np.percentile(hs[a:b], mh_percentil)
    P.MH_max, P.MH_p = mh_max, mh_p
    malas = P.n_pts < min_pts_reb
    for a in (P.MH_max, P.MH_p, P.MW_max, P.MW_p):
        a[malas] = np.nan
    return P


def anchuras_franja(P):
    """Anchura de cada celda (rebanada x franja): estricta (tmax - tmin) y robusta (percentiles de t).
    NaN si la franja no está ocupada (< min_pts_franja puntos) o si la rebanada no es válida (< min_pts_reb puntos)."""
    valida = np.isfinite(P.MH_max)[:, None]
    ae = np.where(P.ocupada & valida, P.tmax - P.tmin, np.nan)
    ap = np.where(P.ocupada & valida, P.thi - P.tlo, np.nan)
    return ae, ap


def metricas_copa(P, dist_hileras=np.nan):
    """Métricas de copa por rebanada, a partir de la misma malla rebanada x franja que MH y MW (sin leer otra vez la nube).
      area_seccion_m2   : área de la sección transversal = suma de (anchura de franja x franja) de las franjas ocupadas
                          (estricta; _p con la anchura robusta). Solo copa por encima de 'umbral'.
      vol_copa_m3       : area_seccion_m2 x rebanada (volumen de copa de la rebanada).
      h_base_copa_m     : altura de inicio de copa = borde inferior de la franja ocupada más baja. Si la copa baja de
                          'umbral' (0,3 m), vale 'umbral' (no se mide por debajo).
      prof_copa_m       : MH_max - h_base_copa_m (profundidad de copa).
      anchura_max_m     : anchura de la franja más ancha; h_anchura_max_m = altura (centro) de esa franja.
      TRV_m3_ha         : Tree Row Volume = MH x MW x 10 000 / distancia entre hileras (estricta; TRV_p con MH_p, MW_p).
    Devuelve un dict de arrays (n_reb,)."""
    ae, ap = anchuras_franja(P)
    valida = np.isfinite(P.MH_max)
    occ = np.isfinite(ae)
    n_reb = P.n_reb
    area_e = np.where(valida, np.nansum(np.where(occ, ae, 0.0), axis=1) * P.franja, np.nan)
    area_p = np.where(valida, np.nansum(np.where(np.isfinite(ap), ap, 0.0), axis=1) * P.franja, np.nan)
    hay = occ.any(axis=1) & valida
    j_min = np.where(hay, np.argmax(occ, axis=1), 0)
    h_base = np.where(hay, P.umbral + j_min * P.franja, np.nan)
    j_anch = np.where(hay, np.argmax(np.where(occ, ae, -np.inf), axis=1), 0)
    anch_max = np.where(hay, ae[np.arange(n_reb), j_anch], np.nan)
    h_anch = np.where(hay, P.umbral + (j_anch + 0.5) * P.franja, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        trv = P.MH_max * P.MW_max * 1e4 / dist_hileras
        trv_p = P.MH_p * P.MW_p * 1e4 / dist_hileras
    return dict(area_seccion_m2=area_e, area_seccion_p_m2=area_p,
                vol_copa_m3=area_e * P.rebanada, vol_copa_p_m3=area_p * P.rebanada,
                h_base_copa_m=h_base, prof_copa_m=P.MH_max - h_base,
                anchura_max_m=anch_max, h_anchura_max_m=h_anch,
                TRV_m3_ha=trv, TRV_p_m3_ha=trv_p)


def perfil_anchura(P, grupos=None):
    """Perfil de anchura por altura: para cada franja (h_ini, h_fin) y cada grupo de rebanadas (None = toda la hilera;
    si no, un array (n_reb,) con el número de grupo de cada rebanada, p. ej. el tramo de 1,5 m), sobre las rebanadas válidas:
      n_rebanadas, n_ocupadas, frac_ocupada,
      anchura_media_m      : media de la anchura de la franja en las rebanadas donde está ocupada (estricta),
      anchura_media_p_m    : idem con la anchura robusta,
      anchura_media_con_huecos_m : media contando 0 donde la franja está vacía (= anchura_media_m x frac_ocupada).
    Solo hasta la franja ocupada más alta del grupo. Devuelve una lista de dict (una fila por grupo y franja)."""
    ae, ap = anchuras_franja(P)
    valida = np.isfinite(P.MH_max)
    g = np.zeros(P.n_reb, np.int64) if grupos is None else np.asarray(grupos)
    filas = []
    for k in np.unique(g[valida]) if valida.any() else []:
        m = valida & (g == k)
        e, p_ = ae[m], ap[m]
        occ = np.isfinite(e)
        if not occ.any():
            continue
        j_top = int(np.nonzero(occ.any(axis=0))[0].max())
        for j in range(j_top + 1):
            n_o = int(occ[:, j].sum())
            filas.append(dict(grupo=k.item() if hasattr(k, "item") else k, franja=j,
                              h_ini_m=round(P.umbral + j * P.franja, 3), h_fin_m=round(P.umbral + (j + 1) * P.franja, 3),
                              n_rebanadas=int(m.sum()), n_ocupadas=n_o, frac_ocupada=n_o / int(m.sum()),
                              anchura_media_m=float(np.nanmean(e[:, j])) if n_o else np.nan,
                              anchura_media_p_m=float(np.nanmean(p_[:, j])) if np.isfinite(p_[:, j]).any() else np.nan,
                              anchura_media_con_huecos_m=float(np.nansum(e[:, j]) / int(m.sum()))))
    return filas


def envolvente(P, estadistico="max", relleno=False):
    """(tmin, tmax, ocupada) por celda según config (estricta o percentil).

    relleno=True (2026-11-16, como Escolà et al. 2017): en cada rebanada, las franjas SIN puntos suficientes que quedan
    ENTRE la franja ocupada más baja y la más alta pasan a ser copa, con tmin/tmax interpolados linealmente en altura
    entre las franjas ocupadas vecinas. Así los rayos que cruzan un hueco interior de la copa cuentan en Bt (y en Bi si
    lo atraviesan). Por encima de la franja más alta y por debajo de la más baja no se rellena nada (cielo y suelo fuera).
    """
    if estadistico == "max":
        tmin, tmax, occ = P.tmin, P.tmax, P.ocupada
    else:
        tmin, tmax, occ = P.elo, P.ehi, P.ocupada & np.isfinite(P.elo)
    if not relleno:
        return tmin, tmax, occ
    tmin, tmax, occ = tmin.copy(), tmax.copy(), occ.copy()
    for i in range(occ.shape[0]):
        jj = np.nonzero(occ[i])[0]
        if len(jj) < 2:
            continue
        hueco = np.arange(jj[0], jj[-1] + 1)
        hueco = hueco[~occ[i, hueco]]
        if len(hueco) == 0:
            continue
        tmin[i, hueco] = np.interp(hueco, jj, tmin[i, jj])
        tmax[i, hueco] = np.interp(hueco, jj, tmax[i, jj])
        occ[i, hueco] = True
    return tmin, tmax, occ


# --------------------------------------------------------------------------------------
# Clasificación de rayos
# --------------------------------------------------------------------------------------
def clasificar_rayos(O, D, R, env, P, dtm_fn, paso=0.025, az_max_deg=25.0, el_rango_deg=(-60, 15),
                     origen_t_abs=(0.8, 3.0), alcance_max=30.0, tol_extremo=None):
    """Clasifica rayos respecto al volumen de copa.

    O (N,3): origen en (s, t, z);  D (N,3): dirección unitaria en (s, t, z);
    R (N,): distancia al eco (np.inf = disparo sin eco);  dtm_fn(s, t) -> z del suelo.
    env = (tmin, tmax, ocupada) con forma (n_reb, nf) de 'perfiles_rebanada'.
    Devuelve dict con arrays por rayo: valido (pasa filtros de geometría), entra, atraviesa
    (Bi), lado (+1/-1), i, j (celda de entrada). Un rayo cuenta en Bt si entra en la envolvente
    de copa por su lado, y en Bi si sale por el otro lado sin que su eco caiga dentro.
    """
    tmin, tmax, occ = env
    n_reb, nf = occ.shape
    N = len(R)
    out = dict(valido=np.zeros(N, bool), entra=np.zeros(N, bool), atraviesa=np.zeros(N, bool),
               lado=np.zeros(N, np.int8), i=np.full(N, -1, np.int32), j=np.full(N, -1, np.int32))
    if N == 0:
        return out
    side = np.sign(O[:, 1]).astype(np.int8)
    ta = np.abs(O[:, 1])
    el = np.degrees(np.arcsin(np.clip(D[:, 2], -1, 1)))
    hor = np.hypot(D[:, 0], D[:, 1])
    hacia = (D[:, 1] * side) < -1e-6
    with np.errstate(divide="ignore", invalid="ignore"):
        az = np.degrees(np.arctan2(np.abs(D[:, 0]), -D[:, 1] * side))
    v = (ta >= origen_t_abs[0]) & (ta <= origen_t_abs[1]) & hacia & (az <= az_max_deg) & \
        (el >= el_rango_deg[0]) & (el <= el_rango_deg[1])
    out["valido"] = v
    out["lado"] = side
    if not v.any():
        return out
    iv = np.nonzero(v)[0]
    Rv = np.minimum(R[iv], alcance_max)
    Ov, Dv, sv = O[iv], D[iv], side[iv]
    # rango de t que contiene la envolvente
    tm = np.nanmax(np.abs(np.where(occ, np.stack([tmin, tmax]), np.nan))) if occ.any() else 0.0
    if not np.isfinite(tm) or tm == 0:
        return out
    K = int(np.ceil(2 * tm / paso)) + 1
    tk = (sv[:, None] * (tm - paso * np.arange(K))[None, :]).astype(np.float64)    # de lado propio a lado opuesto
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = (tk - Ov[:, 1:2]) / Dv[:, 1:2]
    sk = Ov[:, 0:1] + lam * Dv[:, 0:1]
    zk = Ov[:, 2:3] + lam * Dv[:, 2:3]
    hk = zk - dtm_fn(sk, tk)
    P0 = P

    def dentro(sk_, tk_, hk_, infl):
        ii = np.floor((sk_ - P0.s_ini) / P0.rebanada).astype(np.int64)
        jj = np.floor((hk_ - P0.umbral) / P0.franja + 1e-9).astype(np.int64)
        okc = (ii >= 0) & (ii < n_reb) & (jj >= 0) & (jj < nf)
        ii2, jj2 = np.where(okc, ii, 0), np.where(okc, jj, 0)
        ins = okc & occ[ii2, jj2] & (tk_ >= tmin[ii2, jj2] - infl) & (tk_ <= tmax[ii2, jj2] + infl)
        return ins, ii2, jj2

    ins, ii, jj = dentro(sk, tk, hk, 0.0)
    ins &= (lam > 0) & (lam <= Rv[:, None])
    # extremo del rayo (eco), con envolvente ligeramente inflada para no perder los ecos de la superficie
    tol = paso if tol_extremo is None else tol_extremo
    fin_t = Ov[:, 1] + Rv * Dv[:, 1]
    fin_s = Ov[:, 0] + Rv * Dv[:, 0]
    fin_h = Ov[:, 2] + Rv * Dv[:, 2] - dtm_fn(fin_s, fin_t)
    hay_eco = np.isfinite(R[iv]) & (R[iv] <= alcance_max)
    fins, fi, fj = dentro(fin_s, fin_t, fin_h, tol)
    fins &= hay_eco
    primero = np.argmax(ins, axis=1)
    alguno = ins.any(axis=1)
    entra = alguno | fins
    i_ent = np.where(alguno, ii[np.arange(len(iv)), primero], fi)
    j_ent = np.where(alguno, jj[np.arange(len(iv)), primero], fj)
    atr = entra & ~fins
    out["entra"][iv] = entra
    out["atraviesa"][iv] = atr
    out["i"][iv] = np.where(entra, i_ent, -1)
    out["j"][iv] = np.where(entra, j_ent, -1)
    return out


def acumular_conteos(cl, n_reb, nf, conteos=None):
    """Bt y Bi por (lado, rebanada, franja). lado 0 = izquierdo (+), 1 = derecho (-)."""
    if conteos is None:
        conteos = dict(Bt=np.zeros((2, n_reb, nf), np.int64), Bi=np.zeros((2, n_reb, nf), np.int64))
    m = cl["entra"]
    if m.any():
        l = (cl["lado"][m] < 0).astype(int)
        np.add.at(conteos["Bt"], (l, cl["i"][m], cl["j"][m]), 1)
        a = cl["atraviesa"][m]
        np.add.at(conteos["Bi"], (l[a], cl["i"][m][a], cl["j"][m][a]), 1)
    return conteos


def p_por_rebanada(Bt, Bi, min_rayos_celda=30, min_franjas=3, modo="suma_rayos", j_min=0):
    """P por lado y rebanada. Sólo franjas con Bt >= min_rayos_celda (fiables).
    modo 'suma_rayos' (ficha: Bi/Bt sobre todos los rayos) | 'media_franjas' (cada franja pesa igual).
    Devuelve P_lado (2, n_reb), n_fiables (2, n_reb), P_media (n_reb) [media de los dos lados; NaN si falta uno]."""
    fi = Bt >= min_rayos_celda
    if j_min > 0:                       # variante de sensibilidad: solo franjas j >= j_min (p. ej. desde 0,6 m)
        fi = fi.copy()
        fi[:, :, :j_min] = False
    nfi = fi.sum(axis=2)
    if modo == "suma_rayos":
        bt = np.where(fi, Bt, 0).sum(axis=2)
        bi = np.where(fi, Bi, 0).sum(axis=2)
        with np.errstate(divide="ignore", invalid="ignore"):
            Pl = np.where(bt > 0, bi / bt, np.nan)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(fi, Bi / np.maximum(Bt, 1), np.nan)
        Pl = np.nanmean(np.where(fi, r, np.nan), axis=2) if fi.any() else np.full(Bt.shape[:2], np.nan)
    Pl = np.where(nfi >= min_franjas, Pl, np.nan)
    Pm = np.where(np.isfinite(Pl).all(axis=0), Pl.mean(axis=0), np.nan)
    return Pl, nfi, Pm


def cobertura_P(Bt, MH, umbral, franja, min_rayos_celda=30, j_min=0):
    """Cobertura vertical de P por lado y rebanada (QA-05): solo cuentan las franjas medidas (>= min_rayos_celda),
    P no se extrapola. Devuelve
      h_max (2, n_reb): techo de la franja fiable más alta (m sobre el suelo; NaN si no hay),
      frac  (2, n_reb): (nº de franjas fiables x franja) / (MH - umbral), recortada a 1,
      frac_min (n_reb): mínimo de los dos lados (lo que cubre la P de la rebanada)."""
    fi = Bt >= min_rayos_celda
    if j_min > 0:
        fi = fi.copy()
        fi[:, :, :j_min] = False
    nf = Bt.shape[2]
    idx = np.where(fi, np.arange(nf)[None, None, :], -1)
    jmax = idx.max(axis=2)
    h_max = np.where(jmax >= 0, umbral + (jmax + 1) * franja, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.clip(fi.sum(axis=2) * franja / (MH - umbral)[None, :], 0, 1)
    frac = np.where(np.isfinite(frac), frac, np.nan)
    return h_max, frac, np.minimum(frac[0], frac[1])


def lli(MH, MW, P):
    """LLI = MH * MW * (1 - P)  [m2]."""
    return MH * MW * (1.0 - P)


def filtro_atipicos(x, n_de=2.5):
    """Máscara de valores dentro de media ± n_de desviaciones típicas (ignora NaN)."""
    x = np.asarray(x, float)
    m, d = np.nanmean(x), np.nanstd(x)
    ok = np.isfinite(x)
    if not np.isfinite(d) or d == 0:
        return ok
    return ok & (np.abs(x - m) <= n_de * d)
