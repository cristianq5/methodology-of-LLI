# -*- coding: utf-8 -*-
"""08_segmentacion: puntos de cada hilera = buffer del GPKG (hilera.buffer_m) sobre el tramo de la hilera.

1. Eje inicial = línea del GPKG; eje reajustado por PCA con los puntos de copa (HAG >= reajuste_hag_min_m).
   Longitud: con hilera.extension = "nube" (modo completo) el inicio y el final de la hilera salen de la copa de la nube
   (utils.extension_copa); la línea del GPKG solo identifica la hilera y sirve de referencia. Con "gpkg", la línea.
2. Coordenadas por punto: s (a lo largo, m desde el origen del eje reajustado), t (perpendicular,
   + izquierda), lado = sign(t).
3. Comprobación de que hay puntos por los DOS lados (nube de copa) y de que la trayectoria recorrió las dos calles.
Salida: hileras.npz (un grupo de arrays por hilera), marcos.json y segmentacion.las (dims extra).
"""
import utils  # noqa: F401
import os
import numpy as np
import laspy

from utils import (cargar_config, dir_salida, Paso, cargar_filas, MarcoHilera, reajustar_eje_pca, tramo_s,
                   guardar_json, escribir_las_desde, valores_atributos, extension_desde_nube, extension_copa)
from trayectoria import cargar_trayectoria
from geometria import pasos_sensor

NOMBRE = "08_segmentacion"


def main(cfg=None):
    cfg = cfg or cargar_config()
    out = dir_salida(cfg, NOMBRE)
    ent = os.path.join(dir_salida(cfg, "06_normalizacion"), "recorte_hag.las")
    h = cfg["hilera"]
    u = cfg["normalizacion"]["umbral_altura_m"]
    filas = cargar_filas(cfg)
    # código numérico de cada hilera en la dimensión LAS 'id_hilera' (int16): el id si es un entero; si no,
    # su posición (1, 2, ...) en el archivo vectorial filtrado. La correspondencia se guarda en marcos.json.
    todas_ids = [str(i) for i in cargar_filas(cfg, todas=True)["id_fila"]]
    numericos = all(i.isdigit() and int(i) < 32767 for i in todas_ids)
    codigo = {i: (int(i) if numericos else n + 1) for n, i in enumerate(todas_ids)}
    with Paso(cfg, NOMBRE, entrada=ent, parametros=h) as p:
        las = laspy.read(ent)
        x, y, z = np.asarray(las.x), np.asarray(las.y), np.asarray(las.z)
        hag = np.asarray(las["HeightAboveGround"], dtype=np.float64)
        cls = np.asarray(las.classification)
        p.n_antes = len(x)
        tray, _ = cargar_trayectoria(cfg, cache_dir=dir_salida(cfg, "01_qc"))
        marcos, datos, resumen = {}, {}, []
        idx_total, id_h, s_tot, t_tot, lado_tot = [], [], [], [], []
        for _, f in filas.iterrows():
            m0 = MarcoHilera.de_linea(f.geometry)
            sa, ta = m0.st(np.c_[x, y])
            lo_s, hi_s = tramo_s(cfg, m0)      # QA-01: hilera entera (o zona de búsqueda si la extensión sale de la nube)
            en_buf = (np.abs(ta) <= h["buffer_m"]) & (sa >= lo_s) & (sa <= hi_s) & (cls != 7)
            m = m0
            giro = desp = 0.0
            if h["reajustar_eje"]:
                # el eje se ajusta solo con la copa a lo largo de la línea de referencia (no con objetos de la cabecera)
                cop = en_buf & (hag >= h["reajuste_hag_min_m"]) & (sa >= max(lo_s, 0.0)) & (sa <= min(hi_s, m0.largo))
                if h["reajuste_metodo"] == "pca" and cop.sum() > h["pca_min_puntos"]:
                    m, giro, desp = reajustar_eje_pca(m0, np.c_[x[cop], y[cop]])
            # los límites del tramo se reproyectan al eje reajustado
            sl = m.st(m0.xy(lo_s, 0.0))[0]
            sh = m.st(m0.xy(hi_s, 0.0))[0]
            s, t = m.st(np.c_[x, y])
            ext_info = dict(metodo="gpkg" if not extension_desde_nube(cfg) else "")
            if extension_desde_nube(cfg):
                # longitud real de la hilera = extensión de la copa en la nube; la línea del GPKG solo es referencia
                g0, g1 = m.st(m0.xy(0.0, 0.0))[0], m.st(m0.xy(m0.largo, 0.0))[0]
                cop_e = (np.abs(t) <= h["buffer_m"]) & (s >= sl) & (s <= sh) & (cls != 7) & (hag >= h["extension_hag_min_m"])
                sl, sh, ext_info = extension_copa(cfg, s[cop_e], (sl, sh), (g0, g1))
                if ext_info.get("toca_limite_busqueda"):
                    p.info("AVISO hilera %s: la copa llega al límite de la zona de búsqueda (margen_busqueda_m): %s", f["id_fila"], ext_info)
            sel = (np.abs(t) <= h["buffer_m"]) & (s >= sl) & (s <= sh) & (cls != 7)
            ii = np.nonzero(sel)[0]
            lado = np.sign(t[ii]).astype(np.int8)
            lado[lado == 0] = 1
            clave = f["id_fila"]
            marcos[clave] = dict(marco=m.a_dict(), marco_gpkg=m0.a_dict(), s_rango=[float(sl), float(sh)], extension=ext_info,
                                 atributos=valores_atributos(cfg, f), codigo_las=codigo[clave],
                                 giro_deg=giro, desplazamiento_m=desp)
            datos[clave] = dict(idx=ii, s=s[ii], t=t[ii], hag=hag[ii], x=x[ii], y=y[ii], z=z[ii],
                                gps_time=np.asarray(las.gps_time)[ii], lado=lado)
            # comprobaciones bilaterales
            cop_i = hag[ii] >= u
            n_izq, n_der = int(((lado > 0) & cop_i).sum()), int(((lado < 0) & cop_i).sum())
            pasos = [pp for pp in pasos_sensor(cfg, tray, m, sl, sh)]
            lados_tray = sorted({pp["lado"] for pp in pasos})
            reb = np.floor((s[ii][cop_i] - sl) / cfg["metricas"]["rebanada_m"]).astype(int)
            izq = np.bincount(reb[lado[cop_i] > 0], minlength=int((sh - sl) / .1) + 1)
            der = np.bincount(reb[lado[cop_i] < 0], minlength=int((sh - sl) / .1) + 1)
            ambos = float(np.mean((izq >= cfg["metricas"]["min_puntos_rebanada"] / 2) & (der >= cfg["metricas"]["min_puntos_rebanada"] / 2)))
            resumen.append(dict(id_fila=clave, **valores_atributos(cfg, f), n_puntos=int(len(ii)), n_copa=int(cop_i.sum()),
                                n_copa_izq=n_izq, n_copa_der=n_der, rebanadas_con_ambos_lados_pct=round(100 * ambos, 1),
                                lados_trayectoria=lados_tray, giro_eje_deg=round(giro, 2), desp_eje_m=round(desp, 3),
                                s_rango=[round(sl, 2), round(sh, 2)], largo_m=round(sh - sl, 2),
                                largo_gpkg_m=round(m0.largo, 2), extension=ext_info.get("metodo"),
                                dif_inicio_m=ext_info.get("dif_inicio_m"), dif_final_m=ext_info.get("dif_final_m")))
            p.info("hilera %s: %s", clave, resumen[-1])
            idx_total.append(ii)
            id_h.append(np.full(len(ii), codigo[clave], dtype=np.int16))
            s_tot.append(s[ii].astype(np.float32))
            t_tot.append(t[ii].astype(np.float32))
            lado_tot.append(lado)
        np.savez_compressed(os.path.join(out, "hileras.npz"),
                            **{f"{k}__{c}": v for k, d in datos.items() for c, v in d.items()})
        guardar_json(dict(marcos=marcos, resumen=resumen), os.path.join(out, "marcos.json"))
        idx = np.concatenate(idx_total)
        dims = dict(id_hilera=np.concatenate(id_h), s_hilera=np.concatenate(s_tot), t_hilera=np.concatenate(t_tot),
                    lado=np.concatenate(lado_tot), HAG=hag[idx].astype(np.float32))
        sal = os.path.join(out, "segmentacion.las")
        # una misma nube puede caer en dos buffers: no ocurre con buffer 1,5 m y calles de 3,9 m
        escribir_las_desde(las, idx, sal, dims)
        # PNG de control en planta (QA-17): puntos de copa, eje del GPKG, eje PCA y buffer
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(8, 8), dpi=110)
            for k, d_ in datos.items():
                cop_ = d_["hag"] >= u
                sub_ = np.nonzero(cop_)[0][::max(1, int(cop_.sum() // 40000))]
                ax.scatter(d_["x"][sub_], d_["y"][sub_], s=0.3, c=d_["t"][sub_], cmap="coolwarm", vmin=-1.5, vmax=1.5)
                for nom_, col_ in (("marco_gpkg", "k"), ("marco", "g")):
                    mm = MarcoHilera.de_dict(marcos[k][nom_])
                    s_a, s_b = marcos[k]["s_rango"]
                    xy_ = np.array([mm.xy(s_a, 0.0), mm.xy(s_b, 0.0)])
                    ax.plot(xy_[:, 0], xy_[:, 1], col_ + "-", lw=1.2)
                ax.text(d_["x"].mean(), d_["y"].max(), "hilera %s" % k, fontsize=8)
            ax.set_aspect("equal")
            ax.set_title("08: copa (HAG>=0,3) coloreada por t; archivo vectorial (negro) y eje PCA (verde)")
            fig.savefig(os.path.join(out, "segmentacion_planta.png"))
            plt.close(fig)
        except Exception as e_:   # pragma: no cover
            p.info("PNG omitido: %r", e_)
        p.n_despues = int(len(idx))
        p.salida = out
        p.extra["resumen"] = resumen
    return resumen


if __name__ == "__main__":
    main()
