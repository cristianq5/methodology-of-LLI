# calibracion_global.R: calibración GLOBAL del escáner por bloque del ensayo (réplica de calibracion_global.py).
# Opción pcap.modo_calibracion = "global_bloque" de config/parametros.yaml (la misma que usa Python).
#
# Un único montaje del escáner (giro escáner->cuerpo Er, eje del motor a_s, palanca t_e, escala k, az_off, el_off, dt) para
# todas las calles del bloque; por calle solo el estado del motor (omega, phi, nudos de fase) y la corrección de la
# trayectoria (dW, dV, dR, nudos de dW). Requiere la lectura de la orientación pcap.calibracion_global.lectura_cuaternion
# y la tabla de elevaciones del fichero LAS-1.EP (pcap.calibracion_global.elevaciones = "EP").
# Arranque: la calibración PROPIA de R por calle (paso 09, fase "calibrar") de la calle con menor p90 de validación del bloque
# (entre las válidas si las hay), pasada a la nueva lectura; la fase de cada calle se alinea con ella y se vuelve a buscar.
# Ajuste conjunto por mínimos cuadrados robustos (soft_l1, mismos topes y pesos a priori que calibrar) con
# minpack.lm::nls.lm y un jacobiano por diferencias finitas AGRUPADAS (cada calle solo depende de lo común y de lo suyo:
# ~21 evaluaciones del residuo por jacobiano en lugar de ~175), equivalente al jacobiano disperso de Python.
# Validación: grupos de hold-out de cada calle (semilla 777 + segundo de inicio %% 1000, como calibrar_calle) y el MISMO umbral.
# Diferencias con Python (declaradas): generador aleatorio de R (otros grupos), Levenberg-Marquardt de minpack.lm en lugar de
# scipy 'trf', residuos soft_l1 por transformación (como el resto de extrinseca.R).
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
for (f in c("utils.R", "trayectoria.R", "pcap_hesai.R", "extrinseca.R")) source(file.path(DIR_SCRIPTS, f))

# ---- tabla de elevaciones de LAS-1.EP (Protocol Buffers) ------------------------------------------------
.pb_varint <- function(b, i) {
  r <- 0; s <- 0
  repeat { c <- as.integer(b[i]); i <- i + 1L; r <- r + bitwAnd(c, 127L) * 2^s; s <- s + 7; if (c < 128L) break }
  list(v = r, i = i)
}
.pb_parse <- function(b) {
  i <- 1L; n <- length(b); out <- list()
  while (i <= n) {
    v <- .pb_varint(b, i); k <- v$v; i <- v$i
    f <- k %/% 8; w <- k %% 8
    if (w == 0) { v <- .pb_varint(b, i); val <- v$v; i <- v$i }
    else if (w == 1) { val <- readBin(b[i:(i + 7L)], "double", size = 8, endian = "little"); i <- i + 8L }
    else if (w == 5) { val <- readBin(b[i:(i + 3L)], "double", size = 4, endian = "little"); i <- i + 4L }
    else { v <- .pb_varint(b, i); len <- v$v; i <- v$i; val <- if (len > 0) b[i:(i + len - 1)] else raw(0); i <- as.integer(i + len) }
    out[[length(out) + 1L]] <- list(f = f, w = w, v = val)
  }
  out
}
# Elevación de los 32 canales (grados) = campo 7 de la tabla 7 del registro 0 (radianes). Mismo lector que Python.
leer_elevaciones_ep <- function(ruta) {
  b <- readBin(ruta, "raw", file.info(ruta)$size)
  r0 <- .pb_parse(b)[[1]]$v
  t7 <- Filter(function(z) z$f == 7, .pb_parse(r0))[[1]]$v
  el <- vapply(Filter(function(z) z$f == 7, .pb_parse(t7)), function(z) as.numeric(z$v), 0) * 180 / pi
  if (length(el) != 32) stop(sprintf("LAS-1.EP: se esperaban 32 elevaciones y hay %d", length(el)))
  el
}
tabla_elevaciones <- function(cfg) {
  cg <- cfg$pcap$calibracion_global
  if (identical(cg$elevaciones, "EP")) {
    r <- cfg$rutas$calib_ep
    if (is.null(r) || !file.exists(r)) stop("pcap.calibracion_global.elevaciones = EP pero no se encuentra LAS-1.EP (SCAN/Scanner1/*.EP); indícalo en entradas.brutos_archivos.calib_ep")
    return(leer_elevaciones_ep(r))
  }
  if (identical(cg$elevaciones, "xt32_estandar")) return(ELEV_XT32_DEG)
  stop("pcap.calibracion_global.elevaciones debe ser 'EP' o 'xt32_estandar'")
}
trayectoria_global <- function(cfg) cargar_trayectoria_lectura(cfg, cfg$pcap$calibracion_global$lectura_cuaternion, dir_salida(cfg, "01_qc"))
modo_global <- function(cfg) {
  m <- cfg$pcap$modo_calibracion %||% "por_calle"
  if (!m %in% c("por_calle", "global_bloque")) stop("pcap.modo_calibracion debe ser 'por_calle' o 'global_bloque'")
  identical(m, "global_bloque")
}
dir_calles_bloque <- function(cfg) { c <- cfg; c$.bloque <- NULL; file.path(dir_salida(c, "09_metricas"), "_calles") }
.huella <- function(obj) {                         # hash corto y estable de un objeto (validación de la caché del montaje)
  f <- tempfile(); on.exit(unlink(f))
  writeLines(as.character(toJSON(obj, auto_unbox = TRUE, digits = NA)), f)
  substr(unname(tools::md5sum(f)), 1, 16)
}

# ---- grupos por calle (los guarda 09_metricas en la fase "calibrar") -------------------------------------
ruta_grupos <- function(ruta_calle) sub("\\.json$", "_grupos.rds", ruta_calle)
guardar_grupos <- function(cfg, tray, ruta_calle, inf, w, lim, las_t, las_xyz, hilera) {
  ce <- cfg$pcap$extrinseca; cg <- cfg$pcap$calibracion_global
  om_max <- ce$omega_max_deg_s %||% 0
  if (om_max > 0) { est <- tasa_giro(tray, las_t) <= om_max; las_t <- las_t[est]; las_xyz <- las_xyz[est, , drop = FALSE] }
  cx <- ctx_ventana(las_t, las_xyz, w$t, w$az_deg, w$dist_raw * cfg$pcap$unidad_distancia_m)
  set.seed(cg$semilla_holdout + round(lim[1]) %% 1000L)
  Gh <- grupos_temporales(cx, as.integer(cg$n_grupos_holdout))
  set.seed(cfg$proyecto$semilla_aleatoria)
  Gf <- grupos_temporales(cx, max(20L, as.integer(ce$n_puntos_ajuste %/% 300)))
  if (is.null(Gh) || is.null(Gf)) return(FALSE)
  tmp <- paste0(ruta_grupos(ruta_calle), ".tmp")
  saveRDS(list(calle = basename(ruta_calle), fecha_calle = inf$fecha, hilera = as.character(hilera %||% "prueba"),
               ventana = as.numeric(unlist(inf$ventana_global)), Gf = Gf, Gh = Gh), tmp)
  file.rename(tmp, ruta_grupos(ruta_calle))
  TRUE
}
grupos_validos <- function(ruta_calle, inf) {
  r <- ruta_grupos(ruta_calle)
  if (!file.exists(r)) return(FALSE)
  g <- tryCatch(readRDS(r), error = function(e) NULL)
  !is.null(g) && identical(g$fecha_calle, inf$fecha)
}

# ---- utilidades de grupos y cuaterniones ---------------------------------------------------------------
.ids <- function(G) as.integer(names(G$tm))
.primeros <- function(G, n) sub_grupos(G, head(.ids(G), n))
.concat_G <- function(Gs) {
  off <- 0L; L <- list()
  for (G in Gs) {
    ids <- .ids(G); m <- match(G$gid, ids); mp <- match(G$pgid, ids)        # renumeración 1..n por grupo
    tm <- G$tm; names(tm) <- as.character(seq_along(ids) + off)
    L[[length(L) + 1L]] <- list(t = G$t, az = G$az, ch = G$ch, r = G$r, gid = m + off, P = G$P, pgid = mp + off, tm = tm)
    off <- off + length(ids)
  }
  list(t = unlist(lapply(L, `[[`, "t")), az = unlist(lapply(L, `[[`, "az")), ch = unlist(lapply(L, `[[`, "ch")), r = unlist(lapply(L, `[[`, "r")),
       gid = unlist(lapply(L, `[[`, "gid")), P = do.call(rbind, lapply(L, `[[`, "P")), pgid = unlist(lapply(L, `[[`, "pgid")),
       tm = unlist(lapply(L, `[[`, "tm")))
}
q_media <- function(Q) { v <- eigen(crossprod(Q), symmetric = TRUE)$vectors[, 1]; if (v[4] < 0) v <- -v; matrix(v, 1) }

# ---- ajuste conjunto ------------------------------------------------------------------------------------
# calles de arranque: las n de menor p90 de validación calle a calle (primero las válidas)
candidatos_arranque <- function(D, n) {
  o <- order(!vapply(D, function(d) isTRUE(d$oficial$valida), TRUE), vapply(D, function(d) d$oficial$holdout_p90_m, 0))
  head(o, max(1L, n))
}
ajustar_global <- function(cfg, D, tray_of, tray, elev, log, iref = NULL) {
  cg <- cfg$pcap$calibracion_global; ce <- cfg$pcap$extrinseca
  nw <- length(D)
  w_ini <- vapply(D, function(d) d$ventana[1], 0); w_fin <- vapply(D, function(d) d$ventana[2], 0)
  NL <- as.integer(cg$n_grupos_lsq_por_calle)
  G_full <- .concat_G(lapply(D, function(d) d$Gf))
  Gl <- lapply(D, function(d) .primeros(d$Gf, NL))
  G <- .concat_G(Gl)
  n_res <- vapply(Gl, function(g) nrow(g$P), 0L)
  log(sprintf("calles: %d | grupos: %d (nudos) / %d (mínimos cuadrados); filas: %d", nw, length(G_full$tm), length(G$tm), sum(n_res)))
  er_en_lectura <- function(inf) {
    p <- par_desde_informe(inf); ab <- as.numeric(unlist(inf$ventana_global))
    ts <- seq(ab[1] + 3, ab[2] - 3, length.out = 25)
    qe <- matrix(q_de_rotvec(matrix(p$rotvec_e, 1)), 25, 4, byrow = TRUE)
    M <- qmul(.Rw_q(p, ts), qmul(tray_pose(tray_of, ts)$q, qe))
    q_media(qmul(q_conj(tray_pose(tray, ts)$q), M))
  }
  # arranque: calle con menor p90 de validación (entre las válidas si las hay)
  if (is.null(iref)) iref <- candidatos_arranque(D, 1L)
  inf_r <- D[[iref]]$oficial
  log(sprintf("calle de referencia (arranque): %s (por calle: mediana %.3f, p90 %.3f, %s)", D[[iref]]$calle, inf_r$holdout_mediana_m,
              inf_r$holdout_p90_m, if (isTRUE(inf_r$valida)) "válida" else "no válida"))
  Er_g <- er_en_lectura(inf_r)
  a_s <- as.numeric(inf_r$a_s); a_s <- a_s / sqrt(sum(a_s^2))
  phi <- omega <- dts <- numeric(nw); dW <- matrix(0, nw, 3)
  for (i in seq_len(nw)) {
    inf <- D[[i]]$oficial
    q <- qmul(q_conj(Er_g), er_en_lectura(inf))
    delta <- 2 * atan2(sum(q[1, 1:3] * a_s), q[1, 4])                     # giro alrededor de a_s (swing-twist)
    phi[i] <- as.numeric(inf$phi)[1] - delta; omega[i] <- as.numeric(inf$omega)[1]
    dW[i, ] <- matrix(as.numeric(unlist(inf$dW_m)), ncol = 3)[1, ]; dts[i] <- inf$dt_ms
  }
  par <- list(a_s = a_s, omega = omega, phi = phi, w_ini = w_ini, rotvec_e = as.numeric(q_a_rotvec(Er_g)), t_e = c(0, 0, 0),
              dW = dW, dV = matrix(0, nw, 3), dR = matrix(0, nw, 3), kdW_t = numeric(0), kdW = matrix(0, 0, 3), kn_t = numeric(0), kn_d = numeric(0),
              rotvec_w = c(0, 0, 0), k = 1.0, centro = colMeans(G$P), signo = inf_r$signo_azimut, az_off = 0.0,
              el_off = if (identical(cg$elevaciones, "EP")) 0.0 else median(vapply(D, function(d) d$oficial$el_off_deg, 0)),
              dt = median(dts) * 1e-3, elev_deg = elev, valida = TRUE)
  d0 <- residuos_grupos(G, tray, par, 0.5)
  log(sprintf("arranque: mediana %.3f p90 %.3f", median(d0), percentil(d0, 90)))
  par <- escanear_fase_ventanas(G, tray, par, log, etapas = list(c(4, 180), c(1, 4)), verbose = FALSE)
  d0 <- residuos_grupos(G, tray, par, 0.5)
  log(sprintf("tras buscar la fase de cada calle: mediana %.3f p90 %.3f", median(d0), percentil(d0, 90)))

  NC <- 10L; NP <- 11L    # común: rotvec_e(3) t_e(3) k az_off el_off dt_ms | por calle: omega phi dW(3) dV(3) dR(3)
  empaquetar <- function(p) c(p$rotvec_e, p$t_e, p$k, p$az_off, p$el_off, p$dt * 1e3,
                              as.vector(rbind(p$omega, p$phi, t(p$dW), t(p$dV), t(p$dR))))
  desempaq <- function(x) {
    q <- par; B <- matrix(x[-(1:NC)], nw, NP, byrow = TRUE)
    q$rotvec_e <- x[1:3]; q$t_e <- x[4:6]; q$k <- x[7]; q$az_off <- x[8]; q$el_off <- x[9]; q$dt <- x[10] * 1e-3
    q$omega <- B[, 1]; q$phi <- B[, 2]; q$dW <- B[, 3:5, drop = FALSE]; q$dV <- B[, 6:8, drop = FALSE]; q$dR <- B[, 9:11, drop = FALSE]
    q
  }
  x <- empaquetar(par)
  libres <- rep(TRUE, length(x))
  if (identical(cg$palanca, "cero")) libres[4:6] <- FALSE
  libres[7] <- isTRUE(ce$ajustar_escala_distancia); libres[9] <- isTRUE(ce$ajustar_offset_elevacion)
  B <- matrix(libres[-(1:NC)], nw, NP, byrow = TRUE)
  B[, 3:5] <- isTRUE(ce$ajustar_traslacion_mundo); B[, 6:8] <- isTRUE(ce$ajustar_deriva_dW %||% TRUE); B[, 9:11] <- isTRUE(ce$ajustar_deriva_rotacion %||% TRUE)
  libres[-(1:NC)] <- as.vector(t(B))
  # filas de cada columna (estructura dispersa): residuo de su calle + sus a priori
  lim_r <- c(0L, cumsum(n_res)); nr <- sum(n_res)
  fil_col <- vector("list", length(x)); tipo_col <- integer(length(x))
  for (j in 1:NC) fil_col[[j]] <- seq_len(nr + 3L + 9L * nw)            # lo común afecta a todo (las a priori de t_e incluidas)
  for (i in seq_len(nw)) for (j in seq_len(NP)) {
    col <- NC + (i - 1L) * NP + j
    f <- (lim_r[i] + 1L):lim_r[i + 1L]
    if (j >= 3) { bl <- (j - 3L) %/% 3L; comp <- (j - 3L) %% 3L                # dW, dV, dR -> filas a priori
      f <- c(f, nr + 3L + bl * 3L * nw + (i - 1L) * 3L + comp + 1L) }
    fil_col[[col]] <- f; tipo_col[col] <- j
  }
  il <- which(libres)
  grupos_col <- c(as.list(il[il <= NC]), lapply(seq_len(NP), function(j) il[tipo_col[il] == j]))
  grupos_col <- Filter(length, grupos_col)
  for (tope in ce$dist_max_correspondencia_m) {
    fun <- function(q) {
      xx <- x; xx[libres] <- q; pp <- desempaq(xx)
      .soft_l1(c(residuos_grupos(G, tray, pp, tope), 0.05 * pp$t_e, 0.02 * as.vector(t(pp$dW)), 0.5 * as.vector(t(pp$dV)), 5.0 * as.vector(t(pp$dR))), tope / 3)
    }
    jac <- function(q) {                                     # diferencias finitas agrupadas (paso relativo 1e-4, como diff_step de Python)
      f0 <- fun(q); J <- matrix(0, length(f0), length(q))
      xi <- x; xi[libres] <- q
      h <- 1e-4 * pmax(1, abs(xi))
      h[10] <- max(h[10], 0.05)       # dt (ms): con tiempos GPS ~1,8e9 s la resolución es ~4e-7 s; un paso de 1e-4 ms no cambia nada
      pos <- match(seq_along(x), il)
      for (gc_ in grupos_col) {
        x1 <- xi; x1[gc_] <- x1[gc_] + h[gc_]
        f1 <- fun(x1[libres])
        for (j in gc_) { rr <- fil_col[[j]]; J[rr, pos[j]] <- (f1[rr] - f0[rr]) / h[j] }
      }
      J
    }
    r <- nls.lm(par = x[libres], fn = fun, jac = jac, control = nls.lm.control(maxiter = ce$iteraciones_icp, ftol = 1e-8, ptol = 1e-8))
    x[libres] <- r$par
    dd <- residuos_grupos(G, tray, desempaq(x), 0.5)
    log(sprintf("tope %.2f m -> coste %.1f | mediana %.3f p90 %.3f", tope, 0.5 * sum(r$fvec^2), median(dd), percentil(dd, 90)))
  }
  par <- desempaq(x)
  par <- refinar_fase(G_full, tray, par, w_ini, w_fin, log = log)
  if ((ce$nudos_dW_s %||% 0) > 0) par <- refinar_dW_nudos(G_full, tray, par, w_ini, w_fin, paso = ce$nudos_dW_s, log = log)
  res <- list()
  for (i in seq_len(nw)) {
    d <- D[[i]]; o <- d$oficial
    hd <- residuos_grupos(d$Gh, tray, par, 2.0); fd <- residuos_grupos(d$Gf, tray, par, 2.0)
    med <- median(hd); p90 <- percentil(hd, 90)
    ok <- med < ce$residuo_mediana_max_m && p90 < ce$residuo_p90_max_m
    res[[i]] <- list(calle = d$calle, ventana = c(w_ini[i], w_fin[i]), holdout_mediana_m = med, holdout_p90_m = p90, holdout_n = length(hd),
                     frac_menor_0_15 = mean(hd < 0.15), residuo_rms_inliers_m = if (any(hd < 0.3)) sqrt(mean(hd[hd < 0.3]^2)) else NA,
                     valida = ok, ajuste_mediana_m = median(fd), ajuste_p90_m = percentil(fd, 90),
                     por_calle_mediana_m = o$holdout_mediana_m, por_calle_p90_m = o$holdout_p90_m, por_calle_valida = isTRUE(o$valida),
                     phi_deg = par$phi[i] * 180 / pi, periodo_s = 2 * pi / abs(par$omega[i]))
    log(sprintf("%s: hold-out %.3f / %.3f -> %s   (por calle %.3f / %.3f %s)", d$calle, med, p90, if (ok) "VALIDA" else "no",
                o$holdout_mediana_m, o$holdout_p90_m, if (isTRUE(o$valida)) "VALIDA" else "no"))
  }
  list(par = par, res = res, ref = D[[iref]]$calle)
}

# ---- montaje -> extrínseca de UNA calle (para P) --------------------------------------------------------
ruta_montaje <- function(cfg) { d <- file.path(dirname(dir_calles_bloque(cfg)), "_global"); dir.create(d, showWarnings = FALSE, recursive = TRUE); file.path(d, "montaje_bloque.json") }
cargar_montaje <- function(cfg) {
  r <- ruta_montaje(cfg)
  if (!file.exists(r)) stop(sprintf("Falta el montaje global del bloque (%s): ejecuta run_pipeline.R desde el paso 09 con pcap.modo_calibracion = global_bloque", r))
  fromJSON(r, simplifyVector = FALSE)
}
par_ventana <- function(mont, i, valida) {
  j <- mont$par; ab <- as.numeric(unlist(mont$calles[[i]]$ventana)); a <- ab[1]; b <- ab[2]
  v <- function(z) as.numeric(unlist(z))
  m3 <- function(z) matrix(v(z), ncol = 3, byrow = TRUE)
  kn_t <- v(j$kn_t); kn_d <- v(j$kn_d); mk <- kn_t >= a - 1e-6 & kn_t <= b + 1e-6
  kdt <- v(j$kdW_t); kd <- if (length(kdt)) m3(j$kdW) else matrix(0, 0, 3); md <- kdt >= a - 1e-6 & kdt <= b + 1e-6
  list(a_s = v(j$a_s), omega = v(j$omega)[i], phi = v(j$phi)[i], w_ini = a, rotvec_e = v(j$rotvec_e), t_e = v(j$t_e),
       dW = m3(j$dW)[i, , drop = FALSE], dV = m3(j$dV)[i, , drop = FALSE], dR = m3(j$dR)[i, , drop = FALSE], rotvec_w = v(j$rotvec_w),
       k = j$k, centro = v(j$centro), signo = j$signo, kn_t = kn_t[mk], kn_d = kn_d[mk], kdW_t = kdt[md], kdW = kd[md, , drop = FALSE],
       az_off = j$az_off, el_off = j$el_off, dt = j$dt, elev_deg = v(j$elev_deg), valida = isTRUE(valida))
}

main_09b <- function(cfg = cargar_config()) {
  cg <- cfg$pcap$calibracion_global
  calles <- dir_calles_bloque(cfg)
  con_paso(cfg, "09b_calibracion_global", entrada = calles, parametros = list(calibracion_global = cg, extrinseca = cfg$pcap$extrinseca), fn = function(p) {
    fs <- sort(list.files(calles, pattern = "^calle_[0-9]+_[0-9]+[.]json$", full.names = TRUE))
    if (!length(fs)) stop(sprintf("No hay calles calibradas en %s", calles))
    D <- list()
    for (f in fs) {
      inf <- fromJSON(f, simplifyVector = TRUE)
      if (!grupos_validos(f, inf)) stop(sprintf("La calle %s no tiene grupos válidos (%s): repite el paso 09 en fase 'calibrar'", basename(f), basename(ruta_grupos(f))))
      g <- readRDS(ruta_grupos(f))
      D[[length(D) + 1L]] <- list(calle = basename(f), ventana = g$ventana, Gf = g$Gf, Gh = g$Gh, oficial = inf)
    }
    D <- D[order(vapply(D, function(d) d$ventana[1], 0))]
    huella <- .huella(list(calles = lapply(D, function(d) c(d$calle, d$oficial$fecha)), cg = cg, ce = cfg$pcap$extrinseca,
                              codigo = vapply(c("calibracion_global.R", "extrinseca.R", "trayectoria.R"), function(f) .huella(readLines(file.path(DIR_SCRIPTS, f), warn = FALSE)), "")))
    ruta <- ruta_montaje(cfg)
    if (file.exists(ruta) && !isTRUE(cfg$recalibrar)) {
      m <- tryCatch(fromJSON(ruta, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(m) && identical(m$huella, huella)) { p$info("montaje global reutilizado (huella %s): %d/%d calles válidas", huella, m$n_validas, length(m$calles)); p$salida <- ruta; return(invisible(m)) }
    }
    t0 <- Sys.time()
    elev <- tabla_elevaciones(cfg)
    p$info("lectura del cuaternión %s | elevaciones %s (canales 1, 16, 32: %s) | palanca %s", nombre_lectura(cg$lectura_cuaternion), cg$elevaciones,
           paste(round(elev[c(1, 16, 32)], 3), collapse = ", "), cg$palanca)
    tray_of <- cargar_trayectoria(cfg, cache_dir = dir_salida(cfg, "01_qc"))
    tray <- trayectoria_global(cfg)
    # varios arranques (pcap.calibracion_global.n_arranques): se queda el montaje con más calles válidas (a igualdad, menor mediana)
    arranques <- list(); a <- NULL; mej <- c(-1, Inf)
    for (ir in candidatos_arranque(D, as.integer(cg$n_arranques %||% 1L))) {
      p$info("=== arranque desde %s", D[[ir]]$calle)
      ai <- ajustar_global(cfg, D, tray_of, tray, elev, function(s) p$info("%s", s), iref = ir)
      nv <- sum(vapply(ai$res, function(z) z$valida, TRUE)); md <- median(vapply(ai$res, function(z) z$holdout_mediana_m, 0))
      arranques[[length(arranques) + 1L]] <- list(calle_referencia = ai$ref, n_validas = nv, holdout_mediana_global_m = md)
      p$info("arranque %s: %d/%d calles válidas, mediana global %.3f m", ai$ref, nv, length(ai$res), md)
      if (nv > mej[1] || (nv == mej[1] && md < mej[2])) { a <- ai; mej <- c(nv, md) }
    }
    p$info("montaje elegido: arranque desde %s", a$ref)
    par <- a$par; res <- a$res
    Er <- q_de_rotvec(matrix(par$rotvec_e, 1))
    m <- list(huella = huella, lectura = c(cg$lectura_cuaternion, list(nombre = nombre_lectura(cg$lectura_cuaternion))), elevaciones = cg$elevaciones,
              palanca = cg$palanca, calle_referencia = a$ref, arranques = arranques, n_validas = sum(vapply(res, function(z) z$valida, TRUE)),
              n_validas_por_calle = sum(vapply(res, function(z) z$por_calle_valida, TRUE)),
              holdout_mediana_global_m = median(vapply(res, function(z) z$holdout_mediana_m, 0)),
              comun = list(rotvec_e_deg = par$rotvec_e * 180 / pi, eje_motor_cuerpo = as.numeric(q_aplicar(Er, matrix(par$a_s, 1))), t_e_m = par$t_e,
                           k = par$k, az_off_deg = par$az_off, el_off_deg = par$el_off, dt_ms = par$dt * 1e3),
              calles = res,
              par = list(a_s = par$a_s, omega = par$omega, phi = par$phi, w_ini = par$w_ini, rotvec_e = par$rotvec_e, t_e = par$t_e,
                         dW = as.vector(t(par$dW)), dV = as.vector(t(par$dV)), dR = as.vector(t(par$dR)), rotvec_w = par$rotvec_w, centro = par$centro,
                         kn_t = par$kn_t, kn_d = par$kn_d, kdW_t = par$kdW_t, kdW = as.vector(t(par$kdW)), elev_deg = par$elev_deg,
                         k = par$k, signo = par$signo, az_off = par$az_off, el_off = par$el_off, dt = par$dt),
              tiempo_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    guardar_json(m, ruta)
    p$info("FIN: %d/%d calles válidas con el montaje global (por calle: %d/%d); palanca %s m; dt %.1f ms; el_off %.2f°", m$n_validas, length(res),
           m$n_validas_por_calle, length(res), paste(round(par$t_e, 3), collapse = ", "), par$dt * 1e3, par$el_off)
    p$extra$calles_validas <- sprintf("%d/%d", m$n_validas, length(res))
    p$salida <- ruta
    invisible(m)
  })
}
