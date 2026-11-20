# -*- coding: utf-8 -*-
"""Recuento de rayos (Bt, Bi) por hilera: vía pcap (disparos reales, incluidos los de distancia 0)
y vía de reserva (trazado desde la trayectoria hasta los puntos del LAS).

Ambas vías pasan por lli_core.clasificar_rayos, de modo que sólo cambia el origen de los rayos.
"""
import numpy as np
import rasterio

import lli_core as lc
import pcap_hesai as ph
import extrinseca as ex


class DTMLookup:
    """Altura de suelo z(x, y) desde el ráster DTM (celda más próxima)."""

    def __init__(self, ruta):
        with rasterio.open(ruta) as r:
            self.a = r.read(1).astype(np.float64)
            self.x0, self.y0 = r.transform.c, r.transform.f
            self.res = r.transform.a
            self.nd = r.nodata
        self.a[self.a == self.nd] = np.nan
        self.zmed = float(np.nanmedian(self.a))

    def z(self, x, y):
        c = np.floor((x - self.x0) / self.res).astype(np.int64)
        r = np.floor((self.y0 - y) / self.res).astype(np.int64)
        ok = (r >= 0) & (r < self.a.shape[0]) & (c >= 0) & (c < self.a.shape[1])
        v = np.full(np.shape(x), self.zmed)
        vv = self.a[np.where(ok, r, 0), np.where(ok, c, 0)]
        v = np.where(ok & np.isfinite(vv), vv, self.zmed)
        return v


def dtm_fn_hilera(dtm, marco):
    """dtm_fn(s, t) para lli_core (convierte (s,t) -> (x,y) con el marco de la hilera)."""
    def f(s, t):
        x = marco.o[0] + s * marco.e[0] + t * marco.n[0]
        y = marco.o[1] + s * marco.e[1] + t * marco.n[1]
        return dtm.z(x, y)
    return f


def rayos_en_marco(marco, origen_w, dir_w):
    """Origen (N,3 mundo) y dirección (N,3 mundo) -> (s, t, z) del marco de hilera."""
    s, t = marco.st(origen_w[:, :2])
    O = np.c_[s, t, origen_w[:, 2]]
    D = np.c_[dir_w[:, :2] @ marco.e, dir_w[:, :2] @ marco.n, dir_w[:, 2]]
    return O, D


def procesar_lote(O, D, R_, hileras, cfg, conteos, estadisticas, el_max_extra=None):
    """Clasifica un lote de rayos (en mundo) contra todas las hileras y acumula conteos."""
    pc = cfg["porosidad"]
    rc = pc["rayos"]
    for h in hileras:
        Os, Ds = rayos_en_marco(h["marco"], O, D)
        # cribado previo: sólo rayos con origen cerca de la hilera
        m = (np.abs(Os[:, 1]) >= pc["origen_t_abs_m"][0]) & (np.abs(Os[:, 1]) <= pc["origen_t_abs_m"][1]) & \
            (Os[:, 0] > h["s_ini"] - pc["margen_s_m"]) & (Os[:, 0] < h["s_fin"] + pc["margen_s_m"])
        if not m.any():
            continue
        idx = np.nonzero(m)[0]
        el = list(rc["elevacion_rango_deg"])
        if el_max_extra is not None:
            el[1] = min(el[1], el_max_extra)
        step = cfg["porosidad"]["rayos_chunk"]
        for a in range(0, len(idx), step):
            ii = idx[a:a + step]
            cl = lc.clasificar_rayos(Os[ii], Ds[ii], R_[ii], h["env"], h["perfil"], h["dtm_fn"],
                                     paso=pc["paso_muestreo_rayo_m"], az_max_deg=rc["azimut_max_respecto_normal_deg"],
                                     el_rango_deg=tuple(el), origen_t_abs=tuple(pc["origen_t_abs_m"]),
                                     alcance_max=rc["alcance_max_m"], tol_extremo=pc["tol_extremo_m"])
            lc.acumular_conteos(cl, h["perfil"].n_reb, h["perfil"].nf, conteos[h["id"]])
            estadisticas[h["id"]]["rayos_candidatos"] += int(len(ii))
            estadisticas[h["id"]]["rayos_validos"] += int(cl["valido"].sum())
            estadisticas[h["id"]]["rayos_entran"] += int(cl["entra"].sum())


def rayos_pcap(cfg, tray, pars, ventanas, hileras, conteos, estadisticas, log=print):
    """Vía principal: disparos decodificados del pcap. ventanas: lista de dict(t, az, dist_raw);
    pars: una extrínseca (par) por ventana."""
    ua = cfg["pcap"]["unidad_distancia_m"]
    dmin, dmax = cfg["pcap"]["dist_min_m"], cfg["pcap"]["dist_max_m"]
    nchunk = cfg["porosidad"]["bloques_por_lote"]      # bloques (disparos de 32 canales) por lote: memoria acotada
    assert len(ventanas) == len(pars), "nº de ventanas (%d) != nº de extrínsecas (%d)" % (len(ventanas), len(pars))
    for iw, (w, par) in enumerate(zip(ventanas, pars)):
        if not par.get("valida", True):
            log("ventana %d EXCLUIDA de P (extrinseca fuera de criterio)" % iw)
            continue
        tb, az, dr = w["t"], w["az"], w["dist_raw"]
        for a in range(0, len(tb), nchunk):
            sl = slice(a, a + nchunk)
            t_, az_, dr_ = tb[sl], az[sl], dr[sl]
            # origen por bloque (posición del escáner en el mundo) y cribado rápido por hilera
            t_ = t_ + par.get("dt", 0.0)
            org = ex.origen_mundo(tray, t_, par)
            utiles = np.zeros(len(t_), bool)
            for h in hileras:
                s, t = h["marco"].st(org[:, :2])
                utiles |= (np.abs(t) >= cfg["porosidad"]["origen_t_abs_m"][0]) & (np.abs(t) <= cfg["porosidad"]["origen_t_abs_m"][1]) &                           (s > h["s_ini"] - cfg["porosidad"]["margen_s_m"]) & (s < h["s_fin"] + cfg["porosidad"]["margen_s_m"])
            if not utiles.any():
                continue
            t_, az_, dr_ = t_[utiles], az_[utiles], dr_[utiles]
            om_max = cfg["pcap"]["extrinseca"].get("omega_max_deg_s", 0)
            if om_max:                                       # sin giros bruscos del cuerpo (la calibración no los cubre)
                est = ex.tasa_giro(tray, t_) <= om_max
                if not est.all():
                    t_, az_, dr_ = t_[est], az_[est], dr_[est]
                if len(t_) == 0:
                    continue
            tabla = par.get("elev_deg")
            tabla = ph.ELEV_XT32_DEG if tabla is None else np.asarray(tabla, float)     # tabla de elevaciones de la calibración
            u = ph.a_direcciones(az_ * par["signo"] + par.get("az_off", 0.0), elev_deg=tabla, offset_elev_deg=par.get("el_off", 0.0))                       # (n,32,3) escáner
            O, dir_w = ex.rayos_mundo(tray, t_, u, par)                   # origen y dirección en el mundo (con el giro del cabezal)
            dist = dr_.astype(np.float64).reshape(-1) * ua * par["k"]
            sin_eco = (dr_.reshape(-1) == 0) | (dist > dmax)
            muy_cerca = (~sin_eco) & (dist < dmin)           # 0 < d < dist_min (operador/equipo): ni Bt ni Bi (QA-15)
            Rr = np.where(sin_eco, np.inf, dist)
            k_ = ~muy_cerca
            procesar_lote(O[k_], dir_w[k_], Rr[k_], hileras, cfg, conteos, estadisticas)
    return conteos


def rayos_las(cfg, tray, pars, lims, t_las, xyz_las, hileras, conteos, estadisticas):
    """Vía de contraste: rayo desde el origen del sensor hasta cada punto del LAS adquirido en las ventanas.
    INDEPENDIENTE de la extrínseca (causa raíz del NaN en 10 hileras: antes se saltaba las ventanas "no válidas" y usaba
    dW/Er de un ajuste fallido): el origen es la posición de la trayectoria post-SLAM (la palanca escáner-cuerpo,
    ~0,1-0,2 m, no se aplica). Sólo se cuentan rayos descendentes u horizontales (exigir_fondo): los ascendentes que
    atraviesan la copa no dejan eco y sesgarían P a la baja. Sin submuestreo."""
    rc = cfg["porosidad"]["rayos"]
    el_max = 0.0 if rc["exigir_fondo"] else None
    paso = cfg["porosidad"]["rayos_chunk"] * 2
    sub = max(1, int(cfg["porosidad"]["rayos_reserva"]["submuestreo_puntos_origen"]))
    for (a0, b0) in lims:
        sel = np.nonzero((t_las >= a0) & (t_las <= b0))[0][::sub]
        for a in range(0, len(sel), paso):
            ii = sel[a:a + paso]
            t_, X = t_las[ii], xyz_las[ii]
            om_max = cfg["pcap"]["extrinseca"].get("omega_max_deg_s", 0)
            if om_max:
                est = ex.tasa_giro(tray, t_) <= om_max
                t_, X = t_[est], X[est]
                if len(t_) == 0:
                    continue
            org, _ = tray.pose(t_)
            v = X - org
            r = np.linalg.norm(v, axis=1)
            ok = (r > cfg["pcap"]["dist_min_m"]) & (r < rc["alcance_max_m"])
            if not ok.any():
                continue
            D = v[ok] / r[ok, None]
            procesar_lote(org[ok], D, r[ok], hileras, cfg, conteos, estadisticas, el_max_extra=el_max)
    return conteos
