# porosidad.R: recuento de rayos (Bt, Bi) por hilera: vía pcap (disparos reales, incluidos los de distancia 0) y vía de
# contraste (trazado desde la trayectoria hasta los puntos del LAS). Réplica de porosidad.py.
# Ambas vías pasan por clasificar_rayos (lli_core.R); solo cambia el origen de los rayos.

dtm_fn_hilera <- function(dtm, marco) {
  dtm$zmed <- median(dtm$a, na.rm = TRUE)
  function(s, t) {
    x <- marco$o[1] + s * marco$e[1] + t * marco$n[1]
    y <- marco$o[2] + s * marco$e[2] + t * marco$n[2]
    dtm_z(dtm, x, y)
  }
}

# Clasifica un lote de rayos (origen O n x 3, dirección D n x 3, distancia R, en el mundo) contra la hilera h y acumula en 'acc'
# (lista de entorno con $cont y $est). Devuelve invisible.
procesar_lote <- function(O, D, R_, h, cfg, accs, el_max_extra = NULL) {
  if (is.environment(accs)) accs <- list(accs)
  pc <- cfg$porosidad; rc <- pc$rayos
  s <- marco_st(h$marco, O[, 1], O[, 2])
  Os <- cbind(s$s, s$t, O[, 3])
  m <- abs(Os[, 2]) >= pc$origen_t_abs_m[1] & abs(Os[, 2]) <= pc$origen_t_abs_m[2] &
    Os[, 1] > h$s_ini - pc$margen_s_m & Os[, 1] < h$s_fin + pc$margen_s_m
  if (!any(m)) return(invisible())
  idx <- which(m)
  Dm <- cbind(D[idx, 1] * h$marco$e[1] + D[idx, 2] * h$marco$e[2], D[idx, 1] * h$marco$n[1] + D[idx, 2] * h$marco$n[2], D[idx, 3])
  Os <- Os[idx, , drop = FALSE]; Rr <- R_[idx]
  el <- rc$elevacion_rango_deg
  if (!is.null(el_max_extra)) el[2] <- min(el[2], el_max_extra)
  step <- pc$rayos_chunk
  for (a in seq(1L, length(idx), by = step)) {
    ii <- a:min(a + step - 1L, length(idx))
    cl <- clasificar_rayos(Os[ii, , drop = FALSE], Dm[ii, , drop = FALSE], Rr[ii], h$env, h$perfil, h$dtm_fn,
                           paso = pc$paso_muestreo_rayo_m, az_max_deg = rc$azimut_max_respecto_normal_deg, el_rango_deg = el,
                           origen_t_abs = pc$origen_t_abs_m, alcance_max = rc$alcance_max_m, tol_extremo = pc$tol_extremo_m)
    for (acc in accs) {
      acc$cont <- acumular_conteos(cl, h$perfil$n_reb, h$perfil$nf, acc$cont)
      acc$est$rayos_candidatos <- acc$est$rayos_candidatos + length(ii)
      acc$est$rayos_validos <- acc$est$rayos_validos + sum(cl$valido)
      acc$est$rayos_entran <- acc$est$rayos_entran + sum(cl$entra)
    }
  }
  invisible()
}

nuevo_acc <- function(n_reb, nf) {
  a <- new.env()
  a$cont <- list(Bt = array(0, c(2, n_reb, nf)), Bi = array(0, c(2, n_reb, nf)))
  a$est <- list(rayos_candidatos = 0, rayos_validos = 0, rayos_entran = 0)
  a
}

# Vía principal: disparos del pcap. ventanas: lista de list(t, az_deg, dist_raw); pars: una extrínseca por ventana;
# 'acc_valida' recibe solo las ventanas válidas (P OFICIAL) y 'acc_todas' todas (alternativa no oficial). valida_ventana: vector lógico.
rayos_pcap <- function(cfg, tray, pars, ventanas, valida_ventana, h, acc_valida, acc_todas, log = function(...) NULL) {
  ua <- cfg$pcap$unidad_distancia_m; dmin <- cfg$pcap$dist_min_m; dmax <- cfg$pcap$dist_max_m
  nchunk <- cfg$porosidad$bloques_por_lote
  om_max <- cfg$pcap$extrinseca$omega_max_deg_s %||% 0
  pc <- cfg$porosidad
  stopifnot(length(ventanas) == length(pars))
  for (iw in seq_along(ventanas)) {
    w <- ventanas[[iw]]; par <- pars[[iw]]
    tb <- w$t; az <- w$az_deg; dr <- w$dist_raw
    n <- length(tb)
    for (a in seq(1L, n, by = nchunk)) {
      sl <- a:min(a + nchunk - 1L, n)
      t_ <- tb[sl] + par$dt
      org <- origen_mundo(tray, t_, par)
      st <- marco_st(h$marco, org[, 1], org[, 2])
      util <- abs(st$t) >= pc$origen_t_abs_m[1] & abs(st$t) <= pc$origen_t_abs_m[2] & st$s > h$s_ini - pc$margen_s_m & st$s < h$s_fin + pc$margen_s_m
      if (!any(util)) next
      keep <- which(util)
      t_ <- t_[keep]; az_ <- az[sl][keep]; dr_ <- dr[sl, , drop = FALSE][keep, , drop = FALSE]; org <- org[keep, , drop = FALSE]
      if (om_max > 0) {
        est <- tasa_giro(tray, t_) <= om_max
        if (!all(est)) { t_ <- t_[est]; az_ <- az_[est]; dr_ <- dr_[est, , drop = FALSE]; org <- org[est, , drop = FALSE] }
        if (!length(t_)) next
      }
      u <- a_direcciones(az_ * par$signo + par$az_off, elev_deg = par$elev_deg %||% ELEV_XT32_DEG, offset_elev_deg = par$el_off)
      dd <- dirs_mundo(tray, t_, u, par)
      nb <- length(t_)
      dist <- dr_ * ua * par$k
      sin_eco <- dr_ == 0 | dist > dmax
      muy_cerca <- !sin_eco & dist < dmin
      Rr <- ifelse(sin_eco, Inf, dist)
      k_ <- !as.vector(muy_cerca)
      O <- cbind(rep(org[, 1], 32L), rep(org[, 2], 32L), rep(org[, 3], 32L))[k_, , drop = FALSE]
      D <- cbind(as.vector(dd$x), as.vector(dd$y), as.vector(dd$z))[k_, , drop = FALSE]
      Rk <- as.vector(Rr)[k_]
      procesar_lote(O, D, Rk, h, cfg, if (valida_ventana[iw]) list(acc_todas, acc_valida) else list(acc_todas))
    }
    log(sprintf("ventana %d procesada (%s bloques)", iw, format(n, scientific = FALSE)))
  }
}

# Vía de contraste: rayo desde el origen del sensor (posición de la trayectoria, sin extrínseca) hasta cada punto LAS de las ventanas.
rayos_las <- function(cfg, tray, lims, t_las, xyz_las, h, acc) {
  rc <- cfg$porosidad$rayos
  el_max <- if (isTRUE(rc$exigir_fondo)) 0.0 else NULL
  paso <- cfg$porosidad$rayos_chunk * 2L
  sub <- max(1L, as.integer(cfg$porosidad$rayos_reserva$submuestreo_puntos_origen))
  om_max <- cfg$pcap$extrinseca$omega_max_deg_s %||% 0
  for (lim in lims) {
    sel <- which(t_las >= lim[1] & t_las <= lim[2]); sel <- sel[seq(1L, length(sel), by = sub)]
    for (a in seq(1L, length(sel), by = paso)) {
      ii <- sel[a:min(a + paso - 1L, length(sel))]
      t_ <- t_las[ii]; X <- xyz_las[ii, , drop = FALSE]
      if (om_max > 0) { est <- tasa_giro(tray, t_) <= om_max; t_ <- t_[est]; X <- X[est, , drop = FALSE]; if (!length(t_)) next }
      org <- tray_pose(tray, t_)$p
      v <- X - org
      r <- sqrt(rowSums(v^2))
      ok <- r > cfg$pcap$dist_min_m & r < rc$alcance_max_m
      if (!any(ok)) next
      procesar_lote(org[ok, , drop = FALSE], v[ok, , drop = FALSE] / r[ok], r[ok], h, cfg, acc, el_max_extra = el_max)
    }
  }
}
