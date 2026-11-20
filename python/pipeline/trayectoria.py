# -*- coding: utf-8 -*-
"""Trayectoria del sensor en el CRS proyectado del ensayo (posición + orientación interpoladas en el tiempo).

Trayectoria principal: POST_SLAM_DATA\\result_trajectory3d_high_freq.txt
  columnas: t, x, y, z (marco local SLAM), qx, qy, qz, qw, ... (cuerpo -> local).
El marco local se lleva a UTM con un ajuste rígido (Kabsch) frente al PosI (lat/lon/h -> EPSG del ensayo),
que es una georreferenciación independiente. Ese ajuste se refina después con la extrínseca
(extrinseca.py), que es la que decide la calidad real frente al LAS.
"""
import os
import numpy as np
from scipy.spatial.transform import Rotation as R, Slerp

import utils  # noqa: F401  (entorno)


def kabsch(P, Q):
    """R, t tales que Q ~ P @ R.T + t."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = Vt.T @ np.diag([1, 1, d]) @ U.T
    return Rm, qc - Rm @ pc


def leer_posi_utm(ruta, epsg_geo, epsg_utm):
    import pandas as pd
    from pyproj import Transformer
    a = pd.read_csv(ruta, sep=r"\s+", skiprows=2, header=None, usecols=range(7)).values
    tr = Transformer.from_crs(epsg_geo, epsg_utm, always_xy=True)
    x, y, z = tr.transform(a[:, 2], a[:, 1], a[:, 3])
    return a[:, 0], np.c_[x, y, z], a[:, 4:7]


class Trayectoria:
    def __init__(self, t, p, rot):
        self.t, self.p, self.rot = t, p, rot       # rot: scipy Rotation (N,) cuerpo -> UTM
        self._slerp = Slerp(t, rot)

    def pose(self, tq):
        tq = np.clip(tq, self.t[0], self.t[-1])
        p = np.c_[np.interp(tq, self.t, self.p[:, 0]), np.interp(tq, self.t, self.p[:, 1]),
                  np.interp(tq, self.t, self.p[:, 2])]
        return p, self._slerp(tq)

    def hueco(self, tq):
        """Distancia (s) al muestreo de trayectoria más cercano (control de huecos)."""
        i = np.clip(np.searchsorted(self.t, tq), 1, len(self.t) - 1)
        return np.minimum(np.abs(tq - self.t[i - 1]), np.abs(self.t[i] - tq))


def cargar_trayectoria(cfg, cache_dir=None):
    """Devuelve (Trayectoria UTM, info del ajuste). Usa caché npz si existe."""
    cache = os.path.join(cache_dir, "trayectoria_utm.npz") if cache_dir else None
    if cache and os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return Trayectoria(z["t"], z["p"], R.from_quat(z["q"])), dict(z["info"].item()) if "info" in z else {}
    hf = np.loadtxt(cfg["rutas"]["trayectoria_hf"])
    t, ploc, q = hf[:, 0], hf[:, 1:4], hf[:, 4:8]
    tp, pp, _ = leer_posi_utm(cfg["rutas"]["trayectoria_posi"], cfg["trayectoria_modelo"]["epsg_posi"], cfg["crs"]["epsg"])
    m = (t >= tp[0]) & (t <= tp[-1])
    Q = np.c_[[np.interp(t[m], tp, pp[:, k]) for k in range(3)]].T
    Rm, tt = kabsch(ploc[m], Q)
    res = np.linalg.norm(ploc[m] @ Rm.T + tt - Q, axis=1)
    info = dict(rms_vs_posi_m=float(np.sqrt((res ** 2).mean())), p50=float(np.percentile(res, 50)),
                p95=float(np.percentile(res, 95)), p99=float(np.percentile(res, 99)),
                giro_z_deg=float(np.degrees(np.arctan2(Rm[1, 0], Rm[0, 0]))), n=int(m.sum()))
    p_utm = ploc @ Rm.T + tt
    rot = R.from_matrix(Rm) * R.from_quat(q)
    if cache:
        np.savez(cache, t=t, p=p_utm, q=rot.as_quat(), info=np.array(info, dtype=object))
    return Trayectoria(t, p_utm, rot), info


# --------------------------------------------------------------------------------------
# Lectura alternativa de la orientación (calibración global por bloque)
# --------------------------------------------------------------------------------------
def rot_lectura(Q, lectura):
    """Q: (N,4) columnas 5-8 del fichero high_freq -> Rotation con la 'lectura' dada.
    lectura = dict(w="ultimo"|"primero", perm=[i,j,k], signos=[sx,sy,sz]):
      w      : posición de la parte escalar (última o primera de las 4 columnas)
      perm   : qué columna de la parte vectorial va a x, y, z
      signos : signo aplicado a x, y, z.
    La lectura original del pipeline (qx qy qz qw) es dict(w="ultimo", perm=[0,1,2], signos=[1,1,1])."""
    if lectura["w"] == "ultimo":
        vec, esc = Q[:, :3], Q[:, 3]
    else:
        vec, esc = Q[:, 1:], Q[:, 0]
    vec = vec[:, list(lectura["perm"])] * np.asarray(lectura["signos"], float)
    return R.from_quat(np.c_[vec, esc])


def nombre_lectura(lectura):
    return "%s_%s_%s" % (lectura["w"], "".join("xyz"[p] for p in lectura["perm"]),
                         "".join("+" if s > 0 else "-" for s in lectura["signos"]))


def cargar_trayectoria_lectura(cfg, lectura, cache_dir):
    """Misma posición UTM y mismo giro local->UTM (Kabsch frente al PosI) que cargar_trayectoria; solo cambia cómo se
    leen las columnas de orientación del fichero high_freq. Caché: <cache_dir>/trayectoria_utm_<lectura>.npz."""
    tray, info = cargar_trayectoria(cfg, cache_dir=cache_dir)
    cache = os.path.join(cache_dir, "trayectoria_utm_%s.npz" % nombre_lectura(lectura))
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return Trayectoria(z["t"], z["p"], R.from_quat(z["q"])), info
    hf = np.loadtxt(cfg["rutas"]["trayectoria_hf"], usecols=range(8))
    assert len(hf) == len(tray.t) and np.allclose(hf[:, 0], tray.t), "la trayectoria en caché no corresponde al fichero high_freq"
    Q = hf[:, 4:8]
    Rm = tray.rot[0] * R.from_quat(Q[0]).inv()           # rot_original = Rm * R_xyzw(q)  ->  Rm (giro local->UTM)
    chk = (Rm * R.from_quat(Q[::5000])).inv() * tray.rot[::5000]
    assert np.degrees(chk.magnitude()).max() < 1e-6, "no se recupera el giro local->UTM de la trayectoria"
    rot = Rm * rot_lectura(Q, lectura)
    tmp = cache + ".tmp.npz"
    np.savez(tmp, t=tray.t, p=tray.p, q=rot.as_quat())
    os.replace(tmp, cache)
    return Trayectoria(tray.t, tray.p, rot), info
