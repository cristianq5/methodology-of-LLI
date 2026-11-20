# geometria.R: recorte por tramo de hilera y pasos del sensor por cada calle (réplica de geometria.py)

# Rectángulo (marco, s0, s1, ancho) como lista; contención analítica (equivale al polígono de shapely)
rect_recorte <- function(cfg, filas) {
  ancho <- cfg$recorte$buffer_recorte_m
  lapply(seq_len(nrow(filas)), function(i) {
    m <- marco_de_linea(sf::st_geometry(filas)[[i]]); ss <- tramo_s(cfg, m)
    list(marco = m, s0 = ss[1], s1 = ss[2], ancho = ancho)
  })
}
dentro_rects <- function(rects, x, y) {
  ok <- rep(FALSE, length(x))
  for (r in rects) { st <- marco_st(r$marco, x, y); ok <- ok | (st$s >= r$s0 & st$s <= r$s1 & abs(st$t) <= r$ancho) }
  ok
}
bbox_rects <- function(rects) {
  xy <- do.call(rbind, lapply(rects, function(r) rbind(marco_xy(r$marco, r$s0, -r$ancho), marco_xy(r$marco, r$s1, -r$ancho),
                                                       marco_xy(r$marco, r$s1, r$ancho), marco_xy(r$marco, r$s0, r$ancho))))
  c(min(xy[, 1]), min(xy[, 2]), max(xy[, 1]), max(xy[, 2]))
}

segmentar_pasos <- function(t_tray, s, t, s0, s1, pad_s, t_abs, hueco_s, min_muestras = 20, min_tramo_s = 3.0, suav_s = 1.0) {
  n <- length(s)
  dt_med <- median(diff(t_tray))
  k <- max(1L, as.integer(round(suav_s / dt_med / 2)))
  v <- numeric(n)
  v[(k + 1):(n - k)] <- s[(2 * k + 1):n] - s[1:(n - 2 * k)]
  m <- s > s0 - pad_s & s < s1 + pad_s & abs(t) >= t_abs[1] & abs(t) <= t_abs[2]
  idx <- which(m)
  if (!length(idx)) return(list())
  lado <- sign(t[idx]); sent <- ifelse(v[idx] >= 0, 1, -1)
  corte <- c(diff(t_tray[idx]) > hueco_s | diff(lado) != 0 | diff(sent) != 0)
  grp <- cumsum(c(0, corte))
  out <- list()
  for (g in split(idx, grp)) {
    if (length(g) < min_muestras || abs(s[g[length(g)]] - s[g[1]]) < min_tramo_s) next
    out[[length(out) + 1]] <- list(t0 = t_tray[g[1]], t1 = t_tray[g[length(g)]], lado = as.integer(sign(median(t[g]))),
                                   sentido = as.integer(sign(s[g[length(g)]] - s[g[1]])), t_med = median(t[g]),
                                   s_ini = s[g[1]], s_fin = s[g[length(g)]])
  }
  out
}

pasos_sensor <- function(cfg, tray, marco, s0, s1) {
  r <- cfg$recorte
  st <- marco_st(marco, tray$p[, 1], tray$p[, 2])
  segmentar_pasos(tray$t, st$s, st$t, s0, s1, r$ventana_tiempo_s_pad_m, c(r$ventana_t_abs_min_m, r$ventana_tiempo_t_max_m), r$ventana_hueco_s)
}

.CACHE_FILAS <- new.env()
.todas <- function(cfg) {
  clave <- paste(cfg$rutas$filas_gpkg, cfg$recorte$modo)
  if (is.null(.CACHE_FILAS[[clave]])) {
    c <- cfg; c$.bloque <- NULL; c$.filas_bloque <- NULL; c$.bloques_sel <- NULL
    .CACHE_FILAS[[clave]] <- cargar_filas(c, todas = TRUE)
  }
  .CACHE_FILAS[[clave]]
}

ventanas_globales <- function(cfg, tray) {
  pasos <- list(); fl <- .todas(cfg)
  for (i in seq_len(nrow(fl))) {
    m <- marco_de_linea(sf::st_geometry(fl)[[i]]); ss <- tramo_s(cfg, m)
    for (p in pasos_sensor(cfg, tray, m, ss[1], ss[2])) {
      p$fila <- fl$id_fila[i]
      # rumbo en el MUNDO: el sentido es relativo a la línea de cada hilera y las vecinas pueden estar digitalizadas al revés
      p$rumbo <- p$sentido * m$e
      pasos[[length(pasos) + 1]] <- p
    }
  }
  pasos <- pasos[order(vapply(pasos, function(p) p$t0, 0))]
  glob <- list()
  for (p in pasos) {
    hecho <- FALSE
    for (j in seq_along(glob)) {
      g <- glob[[j]]
      ov <- min(g$t1, p$t1) - max(g$t0, p$t0)
      if (sum(g$rumbo * p$rumbo) > 0 && ov > 0.5 * min(g$t1 - g$t0, p$t1 - p$t0) && sign(g$t_med) * sign(p$t_med) != 0 && abs(g$t0 - p$t0) < 30) {
        g$t0 <- min(g$t0, p$t0); g$t1 <- max(g$t1, p$t1)
        g$sirve[[length(g$sirve) + 1]] <- list(fila = p$fila, lado = p$lado, sentido = p$sentido)
        glob[[j]] <- g; hecho <- TRUE; break
      }
    }
    if (!hecho) glob[[length(glob) + 1]] <- list(t0 = p$t0, t1 = p$t1, sentido = p$sentido, rumbo = p$rumbo, t_med = p$t_med,
                                             sirve = list(list(fila = p$fila, lado = p$lado, sentido = p$sentido)))
  }
  glob
}

# Devuelve list(ventanas = lista de c(a,b), pasos = lista de pasos por fila)
ventanas_tiempo <- function(cfg, tray, filas) {
  pad_t <- cfg$recorte$ventana_tiempo_pad_s
  ids <- as.character(filas$id_fila)
  glob <- ventanas_globales(cfg, tray)
  glob <- glob[order(vapply(glob, function(g) g$t0, 0))]
  ven <- list()
  for (i in seq_along(glob)) {
    g <- glob[[i]]; a <- g$t0 - pad_t; b <- g$t1 + pad_t
    if (i > 1 && a < ven[[length(ven)]][2]) {
      mid <- 0.5 * (glob[[i - 1]]$t1 + g$t0)
      ven[[length(ven)]][2] <- mid; a <- mid
    }
    ven[[length(ven) + 1]] <- c(a, b)
  }
  sel <- list(); pasos <- list()
  for (i in seq_along(glob)) {
    g <- glob[[i]]
    sirve <- Filter(function(z) z$fila %in% ids, g$sirve)
    if (length(sirve)) {
      sel[[length(sel) + 1]] <- ven[[i]]
      for (z in sirve) pasos[[length(pasos) + 1]] <- list(t0 = g$t0, t1 = g$t1, lado = z$lado, sentido = z$sentido, t_med = g$t_med, fila = z$fila)
    }
  }
  list(ventanas = sel, pasos = pasos)
}
