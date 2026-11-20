# trayectoria.R: trayectoria del sensor en el CRS proyectado del ensayo (réplica de trayectoria.py).
# Posición + orientación (cuaternión x,y,z,w) interpoladas en el tiempo (lineal + slerp).
# El marco local SLAM se lleva a UTM con ajuste rígido (Kabsch) frente al PosI.

# ---- cuaterniones vectorizados (matrices n x 4, orden x y z w) ----------------------
qmul <- function(a, b) {   # producto de Hamilton a*b
  ax <- a[, 1]; ay <- a[, 2]; az <- a[, 3]; aw <- a[, 4]
  bx <- b[, 1]; by <- b[, 2]; bz <- b[, 3]; bw <- b[, 4]
  cbind(aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz)
}
q_de_rotvec <- function(rv) {   # rv: matriz n x 3 (rad)
  rv <- matrix(rv, ncol = 3)
  ang <- sqrt(rowSums(rv^2))
  k <- ifelse(ang < 1e-12, 0.5, sin(ang / 2) / ifelse(ang < 1e-12, 1, ang))
  cbind(rv * k, cos(ang / 2))
}
q_a_matriz <- function(q) {   # lista con los 9 elementos de la matriz de rotación
  n <- sqrt(rowSums(q^2)); x <- q[, 1] / n; y <- q[, 2] / n; z <- q[, 3] / n; w <- q[, 4] / n
  list(m11 = 1 - 2 * (y * y + z * z), m12 = 2 * (x * y - z * w), m13 = 2 * (x * z + y * w),
       m21 = 2 * (x * y + z * w), m22 = 1 - 2 * (x * x + z * z), m23 = 2 * (y * z - x * w),
       m31 = 2 * (x * z - y * w), m32 = 2 * (y * z + x * w), m33 = 1 - 2 * (x * x + y * y))
}
q_aplicar <- function(q, v) {   # v n x 3
  M <- q_a_matriz(q)
  cbind(M$m11 * v[, 1] + M$m12 * v[, 2] + M$m13 * v[, 3],
        M$m21 * v[, 1] + M$m22 * v[, 2] + M$m23 * v[, 3],
        M$m31 * v[, 1] + M$m32 * v[, 2] + M$m33 * v[, 3])
}
q_de_matriz <- function(Rm) {   # matriz 3x3 -> cuaternión (1 x 4)
  tr <- Rm[1, 1] + Rm[2, 2] + Rm[3, 3]
  if (tr > 0) {
    s <- sqrt(tr + 1) * 2
    q <- c((Rm[3, 2] - Rm[2, 3]) / s, (Rm[1, 3] - Rm[3, 1]) / s, (Rm[2, 1] - Rm[1, 2]) / s, 0.25 * s)
  } else if (Rm[1, 1] > Rm[2, 2] && Rm[1, 1] > Rm[3, 3]) {
    s <- sqrt(1 + Rm[1, 1] - Rm[2, 2] - Rm[3, 3]) * 2
    q <- c(0.25 * s, (Rm[1, 2] + Rm[2, 1]) / s, (Rm[1, 3] + Rm[3, 1]) / s, (Rm[3, 2] - Rm[2, 3]) / s)
  } else if (Rm[2, 2] > Rm[3, 3]) {
    s <- sqrt(1 + Rm[2, 2] - Rm[1, 1] - Rm[3, 3]) * 2
    q <- c((Rm[1, 2] + Rm[2, 1]) / s, 0.25 * s, (Rm[2, 3] + Rm[3, 2]) / s, (Rm[1, 3] - Rm[3, 1]) / s)
  } else {
    s <- sqrt(1 + Rm[3, 3] - Rm[1, 1] - Rm[2, 2]) * 2
    q <- c((Rm[1, 3] + Rm[3, 1]) / s, (Rm[2, 3] + Rm[3, 2]) / s, 0.25 * s, (Rm[2, 1] - Rm[1, 2]) / s)
  }
  matrix(q / sqrt(sum(q^2)), 1)
}
q_slerp <- function(q0, q1, f) {
  d <- rowSums(q0 * q1)
  q1 <- q1 * ifelse(d < 0, -1, 1); d <- pmin(abs(d), 1)
  th <- acos(d); sn <- sin(th)
  lin <- sn < 1e-6
  a <- ifelse(lin, 1 - f, sin((1 - f) * th) / ifelse(lin, 1, sn))
  b <- ifelse(lin, f, sin(f * th) / ifelse(lin, 1, sn))
  a * q0 + b * q1
}

kabsch <- function(P, Q) {   # R, t tales que Q ~ P %*% t(R) + t
  pc <- colMeans(P); qc <- colMeans(Q)
  H <- crossprod(sweep(P, 2, pc), sweep(Q, 2, qc))
  s <- svd(H); U <- s$u; V <- s$v
  d <- sign(det(V %*% t(U)))
  Rm <- V %*% diag(c(1, 1, d)) %*% t(U)
  list(R = Rm, t = qc - as.numeric(Rm %*% pc))
}

leer_posi_utm <- function(ruta, epsg_geo, epsg_utm) {
  suppressPackageStartupMessages(library(sf))
  a <- fread(ruta, skip = 2, header = FALSE, select = 1:7)
  pts <- st_as_sf(data.frame(lon = a$V3, lat = a$V2, h = a$V4), coords = c("lon", "lat"), crs = epsg_geo)
  xy <- st_coordinates(st_transform(pts, epsg_utm))
  list(t = a$V1, p = cbind(xy[, 1], xy[, 2], a$V4))
}

# Devuelve lista(t, p, q, info); caché en RDS (cache_dir)
cargar_trayectoria <- function(cfg, cache_dir = NULL) {
  cache <- if (!is.null(cache_dir)) file.path(cache_dir, "trayectoria_utm.rds") else NULL
  if (!is.null(cache) && file.exists(cache)) return(readRDS(cache))
  hf <- as.matrix(fread(cfg$rutas$trayectoria_hf, header = FALSE))
  t <- hf[, 1]; ploc <- hf[, 2:4]; q <- hf[, 5:8]
  pos <- leer_posi_utm(cfg$rutas$trayectoria_posi, cfg$trayectoria_modelo$epsg_posi, cfg$crs$epsg)
  m <- t >= pos$t[1] & t <= pos$t[length(pos$t)]
  Q <- cbind(approx(pos$t, pos$p[, 1], t[m])$y, approx(pos$t, pos$p[, 2], t[m])$y, approx(pos$t, pos$p[, 3], t[m])$y)
  kb <- kabsch(ploc[m, ], Q)
  ajust <- ploc[m, ] %*% t(kb$R) + matrix(kb$t, sum(m), 3, byrow = TRUE)
  res <- sqrt(rowSums((ajust - Q)^2))
  info <- list(rms_vs_posi_m = sqrt(mean(res^2)), p50 = percentil(res, 50), p95 = percentil(res, 95), p99 = percentil(res, 99),
               giro_z_deg = atan2(kb$R[2, 1], kb$R[1, 1]) * 180 / pi, n = sum(m))
  p_utm <- ploc %*% t(kb$R) + matrix(kb$t, nrow(ploc), 3, byrow = TRUE)
  qR <- q_de_matriz(kb$R)
  qtot <- qmul(matrix(qR, nrow(q), 4, byrow = TRUE), q)
  tray <- list(t = t, p = p_utm, q = qtot, info = info)
  if (!is.null(cache)) saveRDS(tray, cache)
  tray
}

# pose interpolada: posición lineal y orientación por slerp (tq recortado al intervalo)
tray_pose <- function(tray, tq) {
  n <- length(tray$t)
  tq <- pmin(pmax(tq, tray$t[1]), tray$t[n])
  i <- pmin(pmax(findInterval(tq, tray$t), 1L), n - 1L)
  f <- (tq - tray$t[i]) / (tray$t[i + 1L] - tray$t[i])
  p <- tray$p[i, , drop = FALSE] * (1 - f) + tray$p[i + 1L, , drop = FALSE] * f
  q <- q_slerp(tray$q[i, , drop = FALSE], tray$q[i + 1L, , drop = FALSE], f)
  list(p = p, q = q)
}

# ---- Lectura alternativa de la orientación (calibración global por bloque; réplica de trayectoria.py) ----------------
# lectura = list(w = "ultimo"|"primero", perm = c(i, j, k) [base 0, como en Python], signos = c(sx, sy, sz)):
#   w: posición de la parte escalar entre las columnas 5-8 del high_freq; perm: columna de la parte vectorial que va a x, y, z;
#   signos: signo aplicado a x, y, z. La lectura original (qx qy qz qw) es list(w = "ultimo", perm = c(0, 1, 2), signos = c(1, 1, 1)).
q_lectura <- function(Q, lectura) {
  if (identical(lectura$w, "ultimo")) { vec <- Q[, 1:3, drop = FALSE]; esc <- Q[, 4] } else { vec <- Q[, 2:4, drop = FALSE]; esc <- Q[, 1] }
  vec <- vec[, as.integer(unlist(lectura$perm)) + 1L, drop = FALSE] * matrix(as.numeric(unlist(lectura$signos)), nrow(Q), 3, byrow = TRUE)
  q <- cbind(vec, esc)
  q / sqrt(rowSums(q^2))
}
nombre_lectura <- function(lectura) {
  sprintf("%s_%s_%s", lectura$w, paste(c("x", "y", "z")[as.integer(unlist(lectura$perm)) + 1L], collapse = ""),
          paste(ifelse(as.numeric(unlist(lectura$signos)) > 0, "+", "-"), collapse = ""))
}
# Misma posición UTM y mismo giro local->UTM (Kabsch frente al PosI) que cargar_trayectoria; solo cambia la lectura de la orientación.
cargar_trayectoria_lectura <- function(cfg, lectura, cache_dir) {
  tray <- cargar_trayectoria(cfg, cache_dir = cache_dir)
  cache <- file.path(cache_dir, sprintf("trayectoria_utm_%s.rds", nombre_lectura(lectura)))
  if (file.exists(cache)) return(readRDS(cache))
  hf <- as.matrix(fread(cfg$rutas$trayectoria_hf, header = FALSE, select = 1:8))
  stopifnot(nrow(hf) == length(tray$t), max(abs(hf[, 1] - tray$t)) < 1e-6)
  Q <- hf[, 5:8]
  qR <- qmul(tray$q[1, , drop = FALSE], q_conj(Q[1, , drop = FALSE] / sqrt(sum(Q[1, ]^2))))     # giro local->UTM (cuaternión)
  chk <- seq(1, nrow(Q), by = 5000)
  qq <- qmul(matrix(qR, length(chk), 4, byrow = TRUE), Q[chk, , drop = FALSE] / sqrt(rowSums(Q[chk, , drop = FALSE]^2)))
  qt <- tray$q[chk, , drop = FALSE]
  cosm <- pmin(abs(rowSums(qq * qt)) / sqrt(rowSums(qq^2) * rowSums(qt^2)), 1)
  stopifnot(max(2 * acos(cosm)) * 180 / pi < 1e-3)        # se recupera el giro local->UTM (ángulo < 0,001°; las normas del fichero no son exactamente 1)
  tray2 <- list(t = tray$t, p = tray$p, q = qmul(matrix(qR, nrow(Q), 4, byrow = TRUE), q_lectura(Q, lectura)), info = tray$info)
  tmp <- paste0(cache, ".tmp"); saveRDS(tray2, tmp); file.rename(tmp, cache)
  tray2
}
