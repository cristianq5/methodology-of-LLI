# -*- coding: utf-8 -*-
"""Decodificador mínimo de pcap Hesai XT32 (protocolo v6.1, retorno único).

Estructura del paquete UDP (820 B) deducida de b1/b3 [DATO]:
  cabecera 12 B | 6 bloques x (azimut u16 [0,01°] + 32 x (dist u16, refl u8, conf u8)) | cola 28 B
  cola: [10]=modo de retorno, [13:19]=fecha UTC (año-1900, mes, día, h, min, s),
        [19:23]=marca de tiempo en µs (LE), [23:27]=nº de secuencia.
Cada bloque es un disparo distinto (paso de azimut 0,18°); distancia 0 = sin eco.
Los ángulos verticales son la tabla estándar XT32 (canal 1=+15° ... canal 32=-16°) [SUP]:
no se ha podido verificar con la documentación de Hesai (manual no accesible).
Sólo lectura del pcap.
"""
import calendar
import struct
import numpy as np

TAM_REG = 16 + 862          # cabecera de registro pcap + trama Ethernet de 862 B
TAM_PAYLOAD = 820
OFF_ETH = 42                # Ethernet 14 + IP 20 + UDP 8
BLOQUES = 6
CANALES = 32
PASO_BLOQUE_S = 1.0 / 20000.0   # 6 bloques = 300 µs [DATO]
ELEV_XT32_DEG = np.arange(15.0, -16.1, -1.0)   # canal 1..32  [SUP]


def _tiempo_host(f, off):
    f.seek(off)
    h = f.read(16)
    ts, tf, _, _ = struct.unpack('<IIII', h)
    return ts + tf * 1e-6


def _resincronizar(f, off):
    """Desde un offset arbitrario, devuelve el primer inicio de registro válido."""
    f.seek(off)
    buf = f.read(8192)
    pat = struct.pack('<II', 862, 862)          # incl_len y orig_len: bytes 8-15 de la cabecera de registro
    i = 0
    while True:
        i = buf.find(pat, i)
        # i apunta a incl_len (byte 8 del registro); la carga útil UDP empieza en 16 + OFF_ETH desde el inicio del
        # registro, es decir, en i - 8 + 16 + OFF_ETH = i + 8 + OFF_ETH
        if i < 0 or i + 8 + OFF_ETH + 2 > len(buf):
            return None
        if i >= 8 and buf[i + 8 + OFF_ETH:i + 8 + OFF_ETH + 2] == b'\xee\xff':
            return off + i - 8
        i += 1


def buscar_offset(f, t_obj, tam_total):
    """Offset del primer registro con tiempo host >= t_obj (bisección en bytes)."""
    lo, hi = 24, tam_total - TAM_REG
    while hi - lo > 4 * TAM_REG:
        mid = (lo + hi) // 2
        o = _resincronizar(f, mid)
        if o is None:
            hi = mid
            continue
        if _tiempo_host(f, o) < t_obj:
            lo = o
        else:
            hi = mid
    return _resincronizar(f, lo) or lo


def leer_disparos(ruta, t0, t1, margen_host_s=2.0):
    """Decodifica los disparos con tiempo de sensor en [t0, t1].

    Devuelve dict con arrays (un elemento por disparo x canal):
      t (s, base de tiempo del sensor), az_deg, canal (0..31), dist_m (0 = sin eco), refl.
    Tiempo = fecha UTC de la cola + µs + 50 µs x bloque.
    """
    import os
    tam = os.path.getsize(ruta)
    with open(ruta, 'rb') as f:
        off = buscar_offset(f, t0 - margen_host_s - 1.0, tam)
        f.seek(off)
        pls = []
        n_host_fuera = 0
        while True:
            h = f.read(16)
            if len(h) < 16:
                break
            ts, tf, incl, _ = struct.unpack('<IIII', h)
            d = f.read(incl)
            if incl != 862:
                continue
            th = ts + tf * 1e-6
            if th > t1 + margen_host_s + 1.0:
                break
            if th < t0 - margen_host_s - 1.0:      # antes de la ventana: no se acumula en memoria
                continue
            pls.append(d[OFF_ETH:])
    a = np.frombuffer(b''.join(pls), dtype=np.uint8).reshape(-1, TAM_PAYLOAD)
    a = a[(a[:, 0] == 0xEE) & (a[:, 1] == 0xFF)]
    npaq = a.shape[0]
    cola = a[:, 792:820]
    # fecha UTC por paquete (año-1900, mes, día, h, m, s)
    f6 = cola[:, 13:19].astype(np.int64)
    sec = np.empty(npaq, dtype=np.float64)
    uniq, inv = np.unique(f6, axis=0, return_inverse=True)
    sec_u = np.array([calendar.timegm((1900 + u[0], u[1], u[2], u[3], u[4], u[5])) for u in uniq], float)
    sec = sec_u[np.ravel(inv)]
    us = (cola[:, 19].astype(np.int64) | cola[:, 20].astype(np.int64) << 8 |
          cola[:, 21].astype(np.int64) << 16 | cola[:, 22].astype(np.int64) << 24)
    tp = sec + us * 1e-6                                   # (npaq,)
    tp = tp[:, None] + PASO_BLOQUE_S * np.arange(BLOQUES)[None, :]   # (npaq,6)
    cuerpo = a[:, 12:12 + BLOQUES * 130].reshape(npaq, BLOQUES, 130)
    az = (cuerpo[:, :, 0].astype(np.uint16) | cuerpo[:, :, 1].astype(np.uint16) << 8) * 0.01
    u = cuerpo[:, :, 2:].reshape(npaq, BLOQUES, CANALES, 4)
    dist = (u[..., 0].astype(np.uint16) | u[..., 1].astype(np.uint16) << 8)
    refl = u[..., 2]
    ok = (tp >= t0) & (tp <= t1)
    ip, ib = np.nonzero(ok)
    return dict(t=tp[ip, ib], az_deg=az[ip, ib], dist_raw=dist[ip, ib], refl=refl[ip, ib])


def a_direcciones(az_deg, elev_deg=ELEV_XT32_DEG, offset_elev_deg=0.0):
    """Vectores unitarios (N,32,3) en el marco del escáner: x adelante (az=0), y izquierda, z arriba.
    Convención Hesai: azimut creciente en sentido horario visto desde arriba (se ajusta con 'signo_az')."""
    el = np.radians(elev_deg + offset_elev_deg)[None, :]
    az = np.radians(az_deg)[:, None]
    ce = np.cos(el)
    # Hesai: x = r cos(el) sin(az), y = r cos(el) cos(az), z = r sin(el)
    return np.stack([ce * np.sin(az), ce * np.cos(az), np.sin(el) * np.ones_like(az)], axis=-1)
