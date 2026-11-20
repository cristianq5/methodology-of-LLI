# extrinseca.R: modelo geométrico de la extrínseca escáner <-> trayectoria y VALIDACIÓN (réplica parcial de extrinseca.py).
#
# QUÉ SE REPLICA EN R: el modelo (giro del cabezal theta(t), Er, Rw, dW lineal+nudos, dt, offsets az/el), la conversión de
# disparos del pcap a rayos en el mundo y la validación del residuo LAS -> eco gemelo (mediana < 0,15 y p90 < 0,30).
# QUÉ NO SE REPLICA: el AJUSTE (ICP + RANSAC de la normal del suelo + búsqueda de fase 360° + mínimos cuadrados con nudos).
#   Es un optimizador no lineal de ~60 nudos de fase y ~36 de dW por ventana (scipy least_squares); reimplementarlo en R
#   sin librerías de optimización/kd-tree y calibrar 15 calles superaba con mucho las 3 h admisibles. Los PARÁMETROS ajustados
#   se leen del intermedio de Python (salidas\09_metricas\_calles\calle_*.json) SOLO para esta parte (decisión declarada).
# La validación sí se recalcula en R con grupos temporales NO usados en el ajuste (hold-out, QA-22).
CUANTO_T <- 1 / 1024
REL_T_GRUPO <- 0.00015

.mat3 <- function(x, nw) { if (is.null(x) || length(x) == 0) return(matrix(0, nw, 3)); matrix(unlist(x), ncol = 3, byrow = !is.matrix(x)) }
par_desde_informe <- function(inf) {
  nw <- length(inf$w_ini)
  kdW <- if (!is.null(inf$kdW) && length(inf$kdW)) matrix(inf$kdW, ncol = 3) else matrix(0, 0, 3)
  list(a_s = as.numeric(inf$a_s), omega = as.numeric(inf$omega), phi = as.numeric(inf$phi), w_ini = as.numeric(inf$w_ini),
       rotvec_e = as.numeric(inf$rotvec_e_deg) * pi / 180, t_e = as.numeric(inf$t_e_m), dW = matrix(inf$dW_m, ncol = 3),
       kdW_t = as.numeric(inf$kdW_t), kdW = kdW,
       dV = if (!is.null(inf$dV_m_s)) matrix(inf$dV_m_s, ncol = 3) else matrix(0, nw, 3),
       dR = if (!is.null(inf$dR_rad_s)) matrix(inf$dR_rad_s, ncol = 3) else matrix(0, nw, 3),
       rotvec_w = as.numeric(inf$rotvec_w_deg) * pi / 180, k = inf$k_escala, centro = as.numeric(inf$centro), signo = inf$signo_azimut,
       kn_t = as.numeric(inf$kn_t), kn_d = as.numeric(inf$kn_d), valida = isTRUE(inf$valida),
       az_off = inf$az_off_deg %||% 0, el_off = inf$el_off_deg %||% 0, dt = (inf$dt_ms %||% 0) * 1e-3)
}

# velocidad angular del cuerpo (grados/s): giro entre pose(t) y pose(t+dt)
tasa_giro <- function(tray, t, dt = 0.1) {
  q1 <- tray_pose(tray, t)$q; q2 <- tray_pose(tray, t + dt)$q
  d <- pmin(abs(rowSums(q1 * q2)), 1)
  2 * acos(d) * 180 / pi / dt
}

.widx <- function(par, t) pmin(pmax(findInterval(t, par$w_ini), 1L), length(par$w_ini))
.theta <- function(par, t) {
  w <- .widx(par, t)
  th <- par$omega[w] * (t - par$w_ini[w]) + par$phi[w]
  if (length(par$kn_t)) th <- th + approx(par$kn_t, par$kn_d, xout = t, rule = 2)$y   # irregularidades de la fase del motor
  th
}
.Rw_q <- function(par, t) {
  w <- .widx(par, t)
  rv <- matrix(par$rotvec_w, length(t), 3, byrow = TRUE) + par$dR[w, , drop = FALSE] * (t - par$w_ini[w])
  q_de_rotvec(rv)
}
.dW <- function(par, t) {
  w <- .widx(par, t)
  out <- par$dW[w, , drop = FALSE] + par$dV[w, , drop = FALSE] * (t - par$w_ini[w])
  if (length(par$kdW_t)) for (j in 1:3) out[, j] <- out[, j] + approx(par$kdW_t, par$kdW[, j], xout = t, rule = 2)$y
  out
}

# Rotación escáner -> mundo por instante: Rw * R(t) * Er * Rot(a_s, -theta)   (cuaternión n x 4)
rot_cuerpo_q <- function(par, t, tray) {
  th <- .theta(par, t)
  spin <- q_de_rotvec(-th %o% par$a_s)
  pose <- tray_pose(tray, t)
  qe <- matrix(q_de_rotvec(matrix(par$rotvec_e, 1)), length(t), 4, byrow = TRUE)
  qmul(.Rw_q(par, t), qmul(pose$q, qmul(qe, spin)))
}

origen_mundo <- function(tray, t, par) {
  pose <- tray_pose(tray, t)
  xr <- q_aplicar(pose$q, matrix(par$t_e, length(t), 3, byrow = TRUE)) + pose$p
  q_aplicar(.Rw_q(par, t), sweep(xr, 2, par$centro)) + matrix(par$centro, length(t), 3, byrow = TRUE) + .dW(par, t)
}

# Direcciones en el mundo (3 matrices n x 32) de los disparos de n bloques; u = a_direcciones(...)
dirs_mundo <- function(tray, t_blq, u, par) {
  M <- q_a_matriz(rot_cuerpo_q(par, t_blq, tray))
  list(x = M$m11 * u$x + M$m12 * u$y + M$m13 * u$z,
       y = M$m21 * u$x + M$m22 * u$y + M$m23 * u$z,
       z = M$m31 * u$x + M$m32 * u$y + M$m33 * u$z)
}

# --------------------------------------------------------------------------------------
# Validación hold-out: residuo de cada punto del LAS al eco gemelo del pcap en grupos NO usados en el ajuste
# --------------------------------------------------------------------------------------
# w: lista(t, az_deg, dist_raw) de la ventana; las_t, las_xyz: puntos LAS (no ruido) de la ventana.
# Grupo = intervalo 1/1024 s del gps_time del LAS (ver NOTAS 4); se evalúan n_grupos al azar con semilla distinta a la del ajuste,
# solo en tramos estables (giro < omega_max). Devuelve residuos (m, tope 0,5) de todos los puntos LAS de los grupos.
residuos_holdout <- function(cfg, tray, par, w, las_t, las_xyz, n_grupos = 500L, semilla = 20260520L, min_pts = 60L, tope = 0.5, excluir = NULL) {
  ua <- cfg$pcap$unidad_distancia_m
  set.seed(semilla)
  k_las <- round(las_t / CUANTO_T)
  tab <- table(k_las)
  ks <- as.numeric(names(tab))[tab >= min_pts]
  if (length(excluir)) ks <- setdiff(ks, as.numeric(excluir))         # sin los grupos usados en el ajuste (k_ajuste)
  om_max <- cfg$pcap$extrinseca$omega_max_deg_s %||% 0
  if (om_max > 0 && length(ks)) ks <- ks[tasa_giro(tray, ks * CUANTO_T) <= om_max]    # solo tramos estables
  if (!length(ks)) return(numeric(0))
  ks <- ks[sample.int(length(ks), min(n_grupos, length(ks)))]
  k_blq <- floor((w$t + REL_T_GRUPO) / CUANTO_T)
  ord_las <- order(k_las)
  res <- list()
  for (k in ks) {
    ib <- which(k_blq == k)
    if (!length(ib)) next
    dm <- w$dist_raw[ib, , drop = FALSE] * ua * par$k
    tt <- w$t[ib] + par$dt
    u <- a_direcciones(w$az_deg[ib] * par$signo + par$az_off, elev_deg = par$elev_deg %||% ELEV_XT32_DEG, offset_elev_deg = par$el_off)
    dr <- dirs_mundo(tray, tt, u, par)
    org <- origen_mundo(tray, tt, par)
    sel <- dm > 0.3
    if (sum(sel) < 20) next
    ex <- as.vector(org[, 1] + dr$x * dm)[as.vector(sel)]; ey <- as.vector(org[, 2] + dr$y * dm)[as.vector(sel)]; ez <- as.vector(org[, 3] + dr$z * dm)[as.vector(sel)]
    il <- which(k_las == k)
    px <- las_xyz[il, 1]; py <- las_xyz[il, 2]; pz <- las_xyz[il, 3]
    # distancia al eco más próximo (fuerza bruta por grupo: ~300 x ~300)
    d2 <- outer(px, ex, "-")^2 + outer(py, ey, "-")^2 + outer(pz, ez, "-")^2
    res[[length(res) + 1L]] <- pmin(sqrt(apply(d2, 1, min)), tope)
  }
  unlist(res)
}


# ======================================================================================
# AJUSTE INDEPENDIENTE DE LA EXTRÍNSECA EN R (réplica de calibrar() de extrinseca.py; NO lee nada de Python)
# Modelo:  X_mundo = Rw (X_ref - c) + c + dW,  X_ref = p(t) + R(t) [ Er Rot(a_s, -theta(t)) (k r u(az, canal)) + Et ]
# Pasos (igual que Python): (1) giro del cabezal a partir de la normal del suelo (RANSAC) -> a_s, omega, phi;
# (2) inclinación de Er por gravedad y guiñada psi por búsqueda de 360°; (3) fase por ventana (barrido +-90° / +-180°);
# (4) mínimos cuadrados robustos (soft_l1) por etapas de tope de correspondencia (minpack.lm::nls.lm con la transformación
#     de residuos de soft_l1); (5) fase del motor a nudos de 1 s; (6) dW a nudos de 5 s; (7) residuo de ajuste y validación hold-out.
# Diferencias deliberadas con Python (documentadas en la comparativa): kd-tree = RANN::nn2 (un solo kd-tree por evaluación, con los
# grupos separados 1e4 m en x); RANSAC con submuestra de 3 000 puntos por trozo de 0,1 s (Python usa todos); optimizador Levenberg-Marquardt
# (minpack.lm) en lugar de scipy least_squares 'trf'; generador aleatorio de R (los grupos elegidos al azar no son los de Python).
# ======================================================================================
suppressPackageStartupMessages({ library(minpack.lm); library(RANN) })

.cross3 <- function(a, b) c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
.n3 <- function(a) a / sqrt(sum(a^2))
q_conj <- function(q) cbind(-q[, 1], -q[, 2], -q[, 3], q[, 4])
q_a_rotvec <- function(q) {
  q <- matrix(q, ncol = 4); q <- q * ifelse(q[, 4] < 0, -1, 1)
  s <- sqrt(rowSums(q[, 1:3, drop = FALSE]^2)); ang <- 2 * atan2(s, q[, 4])
  q[, 1:3, drop = FALSE] * ifelse(s < 1e-12, 2, ang / pmax(s, 1e-12))
}
.rot_min <- function(a, b) {
  a <- .n3(a); b <- .n3(b); v <- .cross3(a, b); s <- sqrt(sum(v^2))
  if (s < 1e-12) return(matrix(c(0, 0, 0, 1), 1))
  q_de_rotvec(matrix(v / s * atan2(s, sum(a * b)), 1))
}
.soft_l1 <- function(r, f) f * sqrt(2 * (sqrt(1 + (r / f)^2) - 1))     # residuo efectivo de la pérdida soft_l1 de scipy

# ---- puntos del eco en el mundo (pose calculada una sola vez) ---------------------------------
mundo_puntos <- function(tray, tt, par, ur) {   # ur: n x 3 = direcciones del escáner * distancia * k
  pose <- tray_pose(tray, tt)
  th <- .theta(par, tt)
  spin <- q_de_rotvec(-th %o% par$a_s)
  qe <- matrix(q_de_rotvec(matrix(par$rotvec_e, 1)), length(tt), 4, byrow = TRUE)
  Rwq <- .Rw_q(par, tt)
  qtot <- qmul(Rwq, qmul(pose$q, qmul(qe, spin)))
  v <- q_aplicar(qtot, ur)
  xr <- q_aplicar(pose$q, matrix(par$t_e, length(tt), 3, byrow = TRUE)) + pose$p
  org <- q_aplicar(Rwq, sweep(xr, 2, par$centro)) + matrix(par$centro, length(tt), 3, byrow = TRUE) + .dW(par, tt)
  org + v
}
puntos_fuente <- function(G, par) {
  tabla <- par$elev_deg %||% ELEV_XT32_DEG               # tabla de elevaciones de la calibración (por defecto la estándar)
  el <- (tabla[G$ch] + par$el_off) * pi / 180
  az <- (G$az * par$signo + par$az_off) * pi / 180
  ce <- cos(el)
  cbind(ce * sin(az), ce * cos(az), sin(el)) * (G$r * par$k)
}
# residuo de cada punto LAS al eco gemelo más próximo de su grupo (tope); un único kd-tree con los grupos desplazados
residuos_grupos <- function(G, tray, par, tope = 0.5, por_grupo = FALSE) {
  W <- mundo_puntos(tray, G$t + par$dt, par, puntos_fuente(G, par))
  W[, 1] <- W[, 1] + G$gid * 1e4
  Q <- G$P; Q[, 1] <- Q[, 1] + G$pgid * 1e4
  d <- pmin(RANN::nn2(W, Q, k = 1)$nn.dists[, 1], tope)
  if (por_grupo) split(d, G$pgid) else d
}
sub_grupos <- function(G, sel) {
  ie <- G$gid %in% sel; ip <- G$pgid %in% sel
  list(t = G$t[ie], az = G$az[ie], ch = G$ch[ie], r = G$r[ie], gid = G$gid[ie], P = G$P[ip, , drop = FALSE], pgid = G$pgid[ip], tm = G$tm[names(G$tm) %in% as.character(sel)])
}

# ---- contexto de una ventana: índices por intervalo de 1/1024 s ---------------------------------------
ctx_ventana <- function(las_t, las_xyz, tb, az, dist_m) {
  list(las_t = las_t, las_xyz = las_xyz, tb = tb, az = az, dist_m = dist_m,
       k_las = round(las_t / CUANTO_T), k_blq = floor((tb + REL_T_GRUPO) / CUANTO_T))
}
grupos_temporales <- function(cx, n_grupos, min_pts = 60L) {
  if (is.null(cx$ind_l)) { cx$ind_l <- split(seq_along(cx$k_las), cx$k_las); cx$ind_b <- split(seq_along(cx$k_blq), cx$k_blq) }
  cnt <- lengths(cx$ind_l)
  ks <- names(cnt)[cnt >= min_pts]
  if (!length(ks)) return(NULL)
  ks <- sample(ks, length(ks))                      # orden aleatorio; se toman los n_grupos primeros con >= 20 ecos
  T <- list(); A <- list(); C <- list(); Rr <- list(); Gi <- list(); Pp <- list(); Pg <- list(); tm <- numeric(0)
  gid <- 0L
  kk <- numeric(0)
  for (k in ks) {
    ib <- cx$ind_b[[k]]
    if (is.null(ib)) next
    dm <- cx$dist_m[ib, , drop = FALSE]
    w <- which(dm > 0.3, arr.ind = TRUE)
    if (nrow(w) < 20) next
    gid <- gid + 1L
    kk[gid] <- as.numeric(k)
    T[[gid]] <- cx$tb[ib][w[, 1]]; A[[gid]] <- cx$az[ib][w[, 1]]; C[[gid]] <- w[, 2]; Rr[[gid]] <- dm[w]; Gi[[gid]] <- rep(gid, nrow(w))
    il <- cx$ind_l[[k]]; Pp[[gid]] <- cx$las_xyz[il, , drop = FALSE]; Pg[[gid]] <- rep(gid, length(il))
    tm[gid] <- mean(T[[gid]])
    if (gid >= n_grupos) break
  }
  if (gid == 0L) return(NULL)
  names(tm) <- seq_len(gid)
  list(t = unlist(T), az = unlist(A), ch = unlist(C), r = unlist(Rr), gid = unlist(Gi), P = do.call(rbind, Pp), pgid = unlist(Pg), tm = tm, k = kk)
}

# ---- paso 1: giro del cabezal a partir de la normal del suelo -----------------------------------------
.plano_ransac <- function(P, it = 300L, tol = 0.03) {
  n_ <- nrow(P); mejor <- 0; bn <- NULL; ba <- NULL
  for (i in seq_len(it)) {
    j <- sample.int(n_, 3L)
    a <- P[j[1], ]; nn <- .cross3(P[j[2], ] - a, P[j[3], ] - a); l <- sqrt(sum(nn^2))
    if (l < 1e-6) next
    nn <- nn / l
    cnt <- sum(abs((P[, 1] - a[1]) * nn[1] + (P[, 2] - a[2]) * nn[2] + (P[, 3] - a[3]) * nn[3]) < tol)
    if (cnt > mejor) { mejor <- cnt; bn <- nn; ba <- a }
  }
  if (is.null(bn)) return(list(n = c(0, 0, 1), inl = 0))
  if (sum(bn * ba) > 0) bn <- -bn                 # normal hacia el sensor
  list(n = bn, inl = mejor / n_)
}
estimar_giro <- function(tb, az, dist_m, w_ini, w_fin, signo, log = function(...) NULL, max_pts = 3000L) {
  Ns <- list(); Ts <- numeric(0); Ws <- integer(0)
  for (iw in seq_along(w_ini)) {
    for (t0 in seq(w_ini[iw] + 0.5, w_fin[iw] - 0.6, by = 0.1)) {
      i0 <- findInterval(t0, tb, left.open = TRUE) + 1L; i1 <- findInterval(t0 + 0.1, tb, left.open = TRUE)
      if (i1 < i0) next
      m <- i0:i1
      d <- dist_m[m, , drop = FALSE]
      ok <- which(d > 0.5 & d < 12)
      if (length(ok) < 2000) next
      if (length(ok) > max_pts) ok <- sample(ok, max_pts)
      u <- a_direcciones(az[m] * signo)
      P <- cbind((u$x * d)[ok], (u$y * d)[ok], (u$z * d)[ok])
      pl <- .plano_ransac(P)
      if (pl$inl > 0.35) { Ns[[length(Ns) + 1]] <- pl$n; Ts <- c(Ts, t0 + 0.05); Ws <- c(Ws, iw) }
    }
  }
  Ns <- do.call(rbind, Ns)
  stopifnot(nrow(Ns) > 20)
  cc <- colMeans(Ns)
  a_s <- svd(sweep(Ns, 2, cc), nu = 0)$v[, 3]
  if (a_s[2] < 0) a_s <- -a_s
  e1 <- .n3(.cross3(a_s, c(1, 0, 0))); e2 <- .cross3(a_s, e1)
  ang <- atan2(Ns %*% e2, Ns %*% e1)[, 1]
  omega <- numeric(length(w_ini)); phi <- numeric(length(w_ini))
  for (iw in seq_along(w_ini)) {
    m <- Ws == iw; t <- Ts[m]; th <- ang[m]
    th <- th[1] + c(0, cumsum((diff(th) + pi) %% (2 * pi) - pi))        # np.unwrap
    k <- rep(TRUE, length(t)); co <- c(0, 0)
    for (i in 1:3) {
      A <- cbind(t - w_ini[iw], 1)
      co <- qr.coef(qr(A), th)
      r <- th - A %*% co
      k <- abs(r) < 20 * pi / 180
      if (sum(k) < 10) break
      t <- t[k]; th <- th[k]
    }
    omega[iw] <- co[1]; phi[iw] <- co[2]
    log(sprintf("ventana %d: omega=%.4f rad/s (periodo %.3f s), n=%d", iw, co[1], 2 * pi / abs(co[1]), sum(m)))
  }
  th_all <- omega[Ws] * (Ts - w_ini[Ws]) + phi[Ws]
  nr <- q_aplicar(q_de_rotvec(-th_all %o% a_s), Ns)
  m_ref <- .n3(colMeans(nr))
  list(a_s = a_s, omega = omega, phi = phi, m_ref = m_ref, cono_n_a = mean(Ns %*% a_s), n_normales = nrow(Ns))
}

# ---- fase por ventana (barrido) -----------------------------------------------------------------------
escanear_fase_ventanas <- function(G, tray, par, log = function(...) NULL, etapas = list(c(3, 90), c(0.5, 3)), verbose = TRUE) {
  par$phi <- as.numeric(par$phi)
  wg <- .widx(par, G$tm)
  for (i in seq_along(par$phi)) {
    sel <- as.integer(names(G$tm))[wg == i]
    if (length(sel) < 3) next
    Gs <- sub_grupos(G, sel); p0 <- par$phi[i]; sc_m <- 9
    for (et in etapas) {
      mejor <- c(9, p0)
      for (dphi in seq(-et[2], et[2], by = et[1]) * pi / 180) {
        par$phi[i] <- p0 + dphi
        sc <- mean(residuos_grupos(Gs, tray, par, 0.5))
        if (sc < mejor[1]) mejor <- c(sc, p0 + dphi)
      }
      p0 <- mejor[2]; par$phi[i] <- p0; sc_m <- mejor[1]
    }
    if (verbose) log(sprintf("ventana %d: fase ajustada a %.1f deg (score %.3f, %d grupos)", i, p0 * 180 / pi, sc_m, length(sel)))
  }
  par
}

# ---- fase del motor a nudos (búsqueda 1D por nudo) ------------------------------------------------------
refinar_fase <- function(G, tray, par, w_ini, w_fin, paso_nudo = 1.0, log = function(...) NULL, pasadas = 2L) {
  kn <- unlist(mapply(function(a, b) seq(a, b + 1e-6, by = paso_nudo), w_ini, w_fin, SIMPLIFY = FALSE))
  par$kn_t <- kn; par$kn_d <- numeric(length(kn))
  gt <- G$tm
  for (ps in seq_len(pasadas) - 1L) {
    for (i in seq_along(kn)) {
      sel <- as.integer(names(gt))[abs(gt - kn[i]) < paso_nudo]
      if (length(sel) < 4) next
      Gs <- sub_grupos(G, sel); base <- par$kn_d[i]; mejor <- c(9, base)
      for (dd in seq(-6, 6, by = if (ps == 0) 0.5 else 0.25) * pi / 180) {
        par$kn_d[i] <- base + dd
        sc <- mean(residuos_grupos(Gs, tray, par, 0.3))
        if (sc < mejor[1]) mejor <- c(sc, par$kn_d[i])
      }
      par$kn_d[i] <- mejor[2]
    }
    d <- residuos_grupos(G, tray, par, 0.5)
    log(sprintf("fase pasada %d: mediana %.3f m, p90 %.3f m", ps + 1L, median(d), percentil(d, 90)))
  }
  par
}

# ---- dW lineal a trozos (nudos) ---------------------------------------------------------------------------
refinar_dW_nudos <- function(G, tray, par, w_ini, w_fin, paso = 5.0, pasadas = 2L, log = function(...) NULL) {
  kn <- unlist(mapply(function(a, b) seq(a, b + 1e-6, by = paso), w_ini, w_fin, SIMPLIFY = FALSE))
  par$kdW_t <- kn; par$kdW <- matrix(0, length(kn), 3)
  gt <- G$tm
  for (ps in seq_len(pasadas)) {
    for (i in seq_along(kn)) {
      sel <- as.integer(names(gt))[abs(gt - kn[i]) < paso]
      if (length(sel) < 6) next
      Gs <- sub_grupos(G, sel)
      f <- function(q) { p2 <- par; p2$kdW[i, ] <- q; .soft_l1(c(residuos_grupos(Gs, tray, p2, 0.5), 0.3 * q), 0.1) }
      r <- tryCatch(nls.lm(par = par$kdW[i, ], fn = f, control = nls.lm.control(maxiter = 8, epsfcn = 1e-4, ftol = 1e-6, ptol = 1e-6)), error = function(e) NULL)
      if (!is.null(r)) par$kdW[i, ] <- r$par
    }
    d <- residuos_grupos(G, tray, par, 0.5)
    log(sprintf("nudos dW pasada %d: mediana %.3f m, p90 %.3f m", ps, median(d), percentil(d, 90)))
  }
  par
}

# ---- calibración completa ----------------------------------------------------------------------------------
# tb, az, dist_raw: disparos del pcap de la ventana; las_t/las_xyz: puntos LAS no ruido de la ventana; semillas: lista de list(rotvec_e, omega)
calibrar <- function(cfg, tray, tb, az, dist_raw, las_t, las_xyz, w_ini, w_fin, log = function(...) NULL, semillas = NULL) {
  ce <- cfg$pcap$extrinseca
  set.seed(cfg$proyecto$semilla_aleatoria)
  dist_m <- dist_raw * cfg$pcap$unidad_distancia_m
  w_ini <- as.numeric(w_ini); w_fin <- as.numeric(w_fin); nw <- length(w_ini)
  om_max <- ce$omega_max_deg_s %||% 0
  if (om_max > 0) {
    est <- tasa_giro(tray, las_t) <= om_max
    log(sprintf("tramos estables (giro <= %.0f deg/s): %.1f %% de los puntos del LAS de la ventana", om_max, 100 * mean(est)))
    las_t <- las_t[est]; las_xyz <- las_xyz[est, , drop = FALSE]
  }
  centro <- colMeans(las_xyz)
  cx <- ctx_ventana(las_t, las_xyz, tb, az, dist_m)
  signo <- cfg$pcap$signo_azimut
  g0 <- estimar_giro(tb, az, dist_m, w_ini, w_fin, signo, log)
  G <- grupos_temporales(cx, 25L)
  usados <- G$k                                       # sellos k de todos los grupos del ajuste (se excluyen del hold-out)
  up_b <- colMeans(do.call(rbind, lapply(seq_len(nw), function(i) q_aplicar(q_conj(tray_pose(tray, 0.5 * (w_ini[i] + w_fin[i]))$q), matrix(c(0, 0, 1), 1)))))
  up_b <- .n3(up_b)
  Er0 <- .rot_min(g0$m_ref, up_b)
  base_par <- list(a_s = g0$a_s, omega = g0$omega, phi = g0$phi, w_ini = w_ini, rotvec_e = c(0, 0, 0), t_e = c(0, 0, 0), dW = matrix(0, nw, 3),
                   dV = matrix(0, nw, 3), dR = matrix(0, nw, 3), kdW_t = numeric(0), kdW = matrix(0, 0, 3), kn_t = numeric(0), kn_d = numeric(0),
                   rotvec_w = c(0, 0, 0), k = 1.0, centro = centro, signo = signo, az_off = 0.0, el_off = 0.0, dt = 0.0, valida = TRUE)
  best <- list(sc = 9, par = NULL)
  for (psi in seq(0, 358, by = 2) * pi / 180) {
    Er <- qmul(q_de_rotvec(matrix(up_b * psi, 1)), Er0)
    par <- base_par; par$rotvec_e <- as.numeric(q_a_rotvec(Er))
    sc <- mean(residuos_grupos(G, tray, par, 0.5))
    if (sc < best$sc) best <- list(sc = sc, par = par)
  }
  log(sprintf("signo_az=%+d: mejor guiñada -> score %.3f (cono n.a=%.3f)", signo, best$sc, g0$cono_n_a))
  par <- best$par
  par$omega <- rep(mean(par$omega), nw)
  if (length(semillas)) {
    Gs <- grupos_temporales(cx, 40L)
    usados <- c(usados, Gs$k)
    cands <- list(list(nom = "busqueda_ciega", par = par))
    for (k in seq_along(semillas)) cands[[length(cands) + 1]] <- list(nom = sprintf("semilla_%d", k), par = modifyList(par, list(rotvec_e = semillas[[k]]$rotvec_e, omega = rep(semillas[[k]]$omega, nw))))
    mej <- list(sc = 9)
    for (cd in cands) {
      c_ <- escanear_fase_ventanas(Gs, tray, cd$par, log, etapas = list(c(4, 180), c(1, 4)), verbose = FALSE)
      sc <- mean(residuos_grupos(Gs, tray, c_, 0.5))
      if (sc < mej$sc) mej <- list(sc = sc, nom = cd$nom, par = c_)
    }
    log(sprintf("arranque elegido: %s (score %.3f; %d candidatos)", mej$nom, mej$sc, length(cands)))
    par <- mej$par
  }
  Gp <- grupos_temporales(cx, 90L)
  usados <- c(usados, Gp$k)
  par <- escanear_fase_ventanas(Gp, tray, par, log)
  G <- grupos_temporales(cx, max(20L, as.integer(ce$n_puntos_ajuste %/% 300)))
  usados <- c(usados, G$k)
  # vector: rotvec_e(3) t_e(3) k(1) omega(1) phi(nw) az_off el_off dt_ms(3) dW(3nw) dV(3nw) dR(3nw)
  empaquetar <- function(p) c(p$rotvec_e, p$t_e, p$k, mean(p$omega), p$phi, p$az_off, p$el_off, p$dt * 1e3, as.vector(t(p$dW)), as.vector(t(p$dV)), as.vector(t(p$dR)))
  desempaq <- function(x, q = par) {
    q$rotvec_e <- x[1:3]; q$t_e <- x[4:6]; q$k <- x[7]; q$omega <- rep(x[8], nw); q$phi <- x[9:(8 + nw)]
    q$az_off <- x[9 + nw]; q$el_off <- x[10 + nw]; q$dt <- x[11 + nw] * 1e-3
    q$dW <- matrix(x[(12 + nw):(11 + 4 * nw)], nw, 3, byrow = TRUE); q$dV <- matrix(x[(12 + 4 * nw):(11 + 7 * nw)], nw, 3, byrow = TRUE)
    q$dR <- matrix(x[(12 + 7 * nw):(11 + 10 * nw)], nw, 3, byrow = TRUE); q
  }
  x <- empaquetar(par)
  libres <- rep(TRUE, length(x))
  libres[4:6] <- FALSE; libres[7] <- isTRUE(ce$ajustar_escala_distancia); libres[10 + nw] <- isTRUE(ce$ajustar_offset_elevacion); libres[11 + nw] <- TRUE
  libres[(12 + nw):length(x)] <- isTRUE(ce$ajustar_traslacion_mundo)
  libres[(12 + 4 * nw):(11 + 7 * nw)] <- isTRUE(ce$ajustar_deriva_dW %||% TRUE)
  libres[(12 + 7 * nw):(11 + 10 * nw)] <- isTRUE(ce$ajustar_deriva_rotacion %||% TRUE)
  for (tope in ce$dist_max_correspondencia_m) {
    fun <- function(q) {
      xx <- x; xx[libres] <- q; pp <- desempaq(xx)
      .soft_l1(c(residuos_grupos(G, tray, pp, tope), 0.05 * pp$t_e, 0.02 * as.vector(pp$dW), 0.5 * as.vector(pp$dV), 5.0 * as.vector(pp$dR)), tope / 3)
    }
    r <- nls.lm(par = x[libres], fn = fun, control = nls.lm.control(maxiter = ce$iteraciones_icp, epsfcn = 1e-8, ftol = 1e-8, ptol = 1e-8))
    x[libres] <- r$par
    log(sprintf("tope %.2f m -> coste %.4f", tope, 0.5 * sum(r$fvec^2)))
  }
  par <- desempaq(x); par$signo <- signo
  log(sprintf("tras ajuste global: mediana %.3f m", median(residuos_grupos(G, tray, par, 2.0))))
  par <- refinar_fase(G, tray, par, w_ini, w_fin, log = log)
  if ((ce$nudos_dW_s %||% 0) > 0) par <- refinar_dW_nudos(G, tray, par, w_ini, w_fin, paso = ce$nudos_dW_s, log = log)
  d <- residuos_grupos(G, tray, par, 2.0)
  informe <- list(signo_azimut = signo, a_s = par$a_s, omega = par$omega, phi = par$phi, w_ini = w_ini,
                  rotvec_e_deg = par$rotvec_e * 180 / pi, t_e_m = par$t_e, dW_m = par$dW, dV_m_s = par$dV, kdW_t = par$kdW_t, kdW = par$kdW,
                  dR_rad_s = par$dR, rotvec_w_deg = par$rotvec_w * 180 / pi, k_escala = par$k, az_off_deg = par$az_off, el_off_deg = par$el_off,
                  dt_ms = par$dt * 1e3, n_grupos = length(unique(G$gid)), n_fuentes = length(d), cono_n_a = g0$cono_n_a,
                  ajuste_mediana_m = median(d), ajuste_p90_m = percentil(d, 90), residuo_mediana_m = median(d), residuo_p90_m = percentil(d, 90),
                  frac_menor_0_05 = mean(d < 0.05), frac_menor_0_15 = mean(d < 0.15), score_busqueda_psi = best$sc,
                  centro = centro, kn_t = par$kn_t, kn_d = par$kn_d, k_ajuste = sort(unique(usados)))
  list(par = par, informe = informe, cx = cx)
}

# ---- calibración de una calle con validación hold-out (QA-22) -----------------------------------------------
# Devuelve el informe (lista, formato del json de Python + holdout_*, ajuste_*, valida = hold-out)
calibrar_calle <- function(cfg, tray, w, lim, las_t, las_xyz, semillas = NULL, log = function(...) NULL) {
  t0 <- Sys.time()
  ce <- cfg$pcap$extrinseca
  r <- calibrar(cfg, tray, w$t, w$az_deg, w$dist_raw, las_t, las_xyz, lim[1], lim[2], log = log, semillas = semillas)
  inf <- r$informe
  par <- par_desde_informe(inf)
  rh <- residuos_holdout(cfg, tray, par, w, las_t, las_xyz, n_grupos = 600L, semilla = 777L + round(lim[1]) %% 1000L,
                         excluir = inf$k_ajuste)
  inf$holdout_mediana_m <- median(rh); inf$holdout_p90_m <- percentil(rh, 90); inf$holdout_n <- length(rh)
  inf$valida_ajuste <- inf$ajuste_mediana_m < ce$residuo_mediana_max_m && inf$ajuste_p90_m < ce$residuo_p90_max_m
  inf$valida <- isTRUE(inf$holdout_mediana_m < ce$residuo_mediana_max_m && inf$holdout_p90_m < ce$residuo_p90_max_m)
  inf$ventana_t <- lim; inf$ventana_global <- lim
  inf$tiempo_s <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  inf$fecha <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  inf$origen <- "ajuste propio en R (extrinseca.R)"
  log(sprintf("calle [%.1f, %.1f]: ajuste mediana %.3f p90 %.3f | hold-out mediana %.3f p90 %.3f -> %s", lim[1], lim[2], inf$ajuste_mediana_m, inf$ajuste_p90_m,
              inf$holdout_mediana_m, inf$holdout_p90_m, if (inf$valida) "VALIDA" else "NO VALIDA"))
  inf
}
