# pcap_hesai.R: decodificador mínimo de pcap Hesai XT32 (réplica de pcap_hesai.py), readBin vectorizado por bloques.
# Registro pcap: cabecera 16 B + trama Ethernet de 862 B (Ethernet 14 + IP 20 + UDP 8 + payload 820 B).
# Payload: cabecera 12 B | 6 bloques x (azimut u16 [0,01°] + 32 x (dist u16, refl u8, conf u8)) | cola 28 B
# cola: [13:19] fecha UTC (año-1900, mes, día, h, min, s), [19:23] µs (LE). Cada bloque = un disparo (50 µs).
# Distancia 0 = sin eco. Ángulos verticales: tabla estándar XT32 (canal 1 = +15° ... canal 32 = -16°) [SUP].
# Los registros de otro tamaño (paquetes de 554 B) se saltan, igual que en Python.

PCAP_TAM_REG <- 878L          # 16 + 862
PCAP_OFF_ETH <- 42L
PCAP_BLOQUES <- 6L
PCAP_CANALES <- 32L
PCAP_PASO_BLOQUE_S <- 1 / 20000
ELEV_XT32_DEG <- seq(15, -16, by = -1)

.u32 <- function(r, st) as.numeric(r[st]) + 256 * as.numeric(r[st + 1L]) + 65536 * as.numeric(r[st + 2L]) + 16777216 * as.numeric(r[st + 3L])

# días desde 1970-01-01 (algoritmo de Hinnant), vectorizado
.dias_civil <- function(y, m, d) {
  y <- y - (m <= 2)
  era <- floor(y / 400)
  yoe <- y - era * 400
  doy <- floor((153 * (m + ifelse(m > 2, -3, 9)) + 2) / 5) + d - 1
  doe <- yoe * 365 + floor(yoe / 4) - floor(yoe / 100) + doy
  era * 146097 + doe - 719468
}

# tiempo host (s) de la cabecera de registro en offset 'off' (0-based)
.tiempo_host <- function(con, off) {
  seek(con, where = off, origin = "start")
  h <- readBin(con, "raw", 16L)
  if (length(h) < 16L) return(NA_real_)
  .u32(h, 1L) + .u32(h, 5L) * 1e-6
}

# primer inicio de registro válido a partir de 'off' (0-based); NULL si no hay
.resincronizar <- function(con, off) {
  seek(con, where = off, origin = "start")
  r <- readBin(con, "raw", 8192L)
  L <- length(r)
  if (L < 100L) return(NULL)
  k <- which(r[1:(L - 80L)] == as.raw(0x5e))
  for (i in k) {
    if (r[i + 1L] == as.raw(3) && r[i + 2L] == as.raw(0) && r[i + 3L] == as.raw(0) && r[i + 4L] == as.raw(0x5e) &&
        r[i + 5L] == as.raw(3) && r[i + 6L] == as.raw(0) && r[i + 7L] == as.raw(0) &&
        i + 8L + PCAP_OFF_ETH + 1L <= L && r[i + 8L + PCAP_OFF_ETH] == as.raw(0xee) && r[i + 8L + PCAP_OFF_ETH + 1L] == as.raw(0xff) && i >= 9L)
      return(off + i - 9L)          # inicio del registro (0-based): la cabecera de registro empieza 8 B antes de incl_len
  }
  NULL
}

# primer registro con tiempo host >= t_obj (bisección en bytes)
buscar_offset <- function(con, t_obj, tam_total) {
  lo <- 24; hi <- tam_total - PCAP_TAM_REG
  while (hi - lo > 4 * PCAP_TAM_REG) {
    mid <- floor((lo + hi) / 2)
    o <- .resincronizar(con, mid)
    if (is.null(o)) { hi <- mid; next }
    th <- .tiempo_host(con, o)
    if (is.na(th)) { hi <- mid; next }
    if (th < t_obj) lo <- o else hi <- mid
  }
  .resincronizar(con, lo) %||% lo
}

# Registros válidos (862 B) de un trozo crudo: devuelve inicios (1-based) y la posición del primero sin procesar
.parse_trozo <- function(r) {
  L <- length(r); starts <- list(); pos <- 1L
  while (pos + PCAP_TAM_REG - 1L <= L) {
    n <- (L - pos + 1L) %/% PCAP_TAM_REG
    st <- pos + PCAP_TAM_REG * (0:(n - 1L))
    ok <- r[st + 8L] == as.raw(0x5e) & r[st + 9L] == as.raw(3) & r[st + 16L + PCAP_OFF_ETH] == as.raw(0xee) & r[st + 17L + PCAP_OFF_ETH] == as.raw(0xff)
    f <- which(!ok)
    if (!length(f)) { starts[[length(starts) + 1L]] <- st; pos <- st[n] + PCAP_TAM_REG; break }
    if (f[1] > 1L) starts[[length(starts) + 1L]] <- st[seq_len(f[1] - 1L)]
    pos <- st[f[1]]
    # registro de otro tamaño: saltar por su longitud incluida si la cabecera es plausible
    if (pos + 15L > L) break
    incl <- .u32(r, pos + 8L); orig <- .u32(r, pos + 12L)
    if (incl == orig && incl > 0 && incl < 70000) {
      if (pos + 16L + incl > L + 0L) break          # registro incompleto al final del trozo
      pos <- pos + 16L + as.integer(incl)
    } else {
      pos <- pos + 1L                                 # sin sincronía: avanzar un byte (raro)
    }
  }
  list(st = unlist(starts), pos = pos)
}

# Decodifica los disparos con tiempo de sensor en [t0, t1]. Devuelve list(t, az_deg, dist_raw [nb x 32 enteros])
leer_disparos <- function(ruta, t0, t1, margen_host_s = 2.0, trozo_regs = 20000L) {
  tam <- file.size(ruta)
  con <- file(ruta, "rb"); on.exit(close(con))
  off <- buscar_offset(con, t0 - margen_host_s - 1.0, tam)
  T <- list(); A <- list(); D <- list()
  fin <- FALSE
  px <- function(k) 59L + k          # fila (1-based) del byte k del payload dentro de las 878 filas del registro
  filas_dist <- lapply(0:5, function(b) px(14L + 130L * b + 4L * (0:31)))
  while (!fin) {
    seek(con, where = off, origin = "start")
    r <- readBin(con, "raw", PCAP_TAM_REG * trozo_regs)
    if (length(r) < PCAP_TAM_REG) break
    pr <- .parse_trozo(r)
    if (length(pr$st)) {
      st <- pr$st
      th_ultimo <- .u32(r, st[length(st)]) + .u32(r, st[length(st)] + 4L) * 1e-6
      th <- .u32(r, st) + .u32(r, st + 4L) * 1e-6
      keep <- th <= t1 + margen_host_s + 1.0
      if (!all(keep)) { st <- st[keep]; fin <- TRUE }
      if (length(st)) {
        m <- matrix(as.integer(r[outer(0:(PCAP_TAM_REG - 1L), st, "+")]), PCAP_TAM_REG, length(st))
        anio <- 1900 + m[px(792L + 13L), ]; mes <- m[px(792L + 14L), ]; dia <- m[px(792L + 15L), ]
        hh <- m[px(792L + 16L), ]; mi <- m[px(792L + 17L), ]; ss <- m[px(792L + 18L), ]
        sec <- .dias_civil(anio, mes, dia) * 86400 + hh * 3600 + mi * 60 + ss
        us <- m[px(792L + 19L), ] + 256 * m[px(792L + 20L), ] + 65536 * m[px(792L + 21L), ] + 16777216 * m[px(792L + 22L), ]
        tp <- sec + us * 1e-6
        magic <- m[px(0L), ] == 238L & m[px(1L), ] == 255L
        for (b in 0:5) {
          tb <- tp + b * PCAP_PASO_BLOQUE_S
          az <- (m[px(12L + 130L * b), ] + 256L * m[px(13L + 130L * b), ]) * 0.01
          fr <- filas_dist[[b + 1L]]
          dist <- m[fr, , drop = FALSE] + 256L * m[fr + 1L, , drop = FALSE]
          ok <- magic & tb >= t0 & tb <= t1
          if (any(ok)) { T[[length(T) + 1L]] <- tb[ok]; A[[length(A) + 1L]] <- az[ok]; D[[length(D) + 1L]] <- t(dist[, ok, drop = FALSE]) }
        }
      }
    }
    off <- off + pr$pos - 1
    if (length(r) < PCAP_TAM_REG * trozo_regs) fin <- TRUE
  }
  t <- unlist(T); az <- unlist(A); d <- do.call(rbind, D)
  o <- order(t)                      # orden temporal (la decodificación va por bloque dentro del paquete)
  list(t = t[o], az_deg = az[o], dist_raw = d[o, , drop = FALSE])
}

# Direcciones unitarias en el marco del escáner: tres matrices n x 32 (x adelante az=0, y izquierda, z arriba)
a_direcciones <- function(az_deg, elev_deg = ELEV_XT32_DEG, offset_elev_deg = 0.0) {
  el <- (elev_deg + offset_elev_deg) * pi / 180
  az <- az_deg * pi / 180
  ce <- cos(el)
  list(x = outer(sin(az), ce), y = outer(cos(az), ce), z = matrix(sin(el), length(az), length(el), byrow = TRUE))
}
