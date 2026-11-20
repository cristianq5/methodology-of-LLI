# 05_dtm.R: ráster de terreno y QA de la altura del sensor (réplica de 05_dtm.py, incluida la corrección de hierba QA-11).
# DTM = IDW 0,25 m con el suelo de las CALLES (se excluye +-excluir_bajo_hilera_m del eje GPKG de cada hilera), luego relleno de huecos.
# Los puntos de suelo > rechazo_sobre_calle_m por encima del DTM se reclasifican a clase 1 (hierba) en recorte_suelo_calles.las.
#
# DIFERENCIAS con Python (que usa PDAL writers.gdal tipo idw + GDAL fillnodata):
#  - IDW propio (potencia 2) con los puntos a distancia <= res*sqrt(2) del centro de celda (igual que PDAL);
#  - celda vacía -> IDW con las celdas con dato en una ventana de idw_window celdas (equivale al window_size de PDAL);
#  - segundo relleno hasta rellenar_huecos_dist_max_m (GDAL fillnodata): aquí IDW con ventana de esa distancia (GDAL busca en
#    4 direcciones; es una aproximación), y por último vecino más próximo.
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R")); source(file.path(DIR_SCRIPTS, "trayectoria.R")); source(file.path(DIR_SCRIPTS, "geometria.R"))

# IDW de celdas vacías a partir de las celdas con dato (matriz), ventana en celdas
rellenar_idw <- function(a, win) {
  nr <- nrow(a); nc <- ncol(a)
  vac <- is.na(a)
  if (!any(vac)) return(a)
  sw <- matrix(0, nr, nc); sz <- matrix(0, nr, nc)
  A0 <- a; A0[vac] <- 0; M <- (!vac) * 1
  for (dy in -win:win) for (dx in -win:win) {
    d2 <- dx * dx + dy * dy
    if (d2 == 0 || d2 > win * win) next
    w <- 1 / d2
    ri <- max(1, 1 - dy):min(nr, nr - dy); ci <- max(1, 1 - dx):min(nc, nc - dx)   # la celda (i,j) toma el valor de (i+dy, j+dx)
    if (!length(ri) || !length(ci) || ri[1] > ri[length(ri)] || ci[1] > ci[length(ci)]) next
    m <- M[ri + dy, ci + dx, drop = FALSE]
    sw[ri, ci] <- sw[ri, ci] + w * m
    sz[ri, ci] <- sz[ri, ci] + w * m * A0[ri + dy, ci + dx]
  }
  rel <- ifelse(vac & sw > 0, sz / pmax(sw, 1e-300), a)
  matrix(rel, nr, nc)
}

dtm_idw <- function(x, y, z, res, win_pdal, win_fill) {
  x0 <- min(x); y1 <- max(y)
  nc <- floor((max(x) - x0) / res) + 1L; nr <- floor((y1 - min(y)) / res) + 1L
  ci <- floor((x - x0) / res); ri <- floor((y1 - y) / res)
  rad <- res * sqrt(2)
  acc <- list()
  for (dy in -1:1) for (dx in -1:1) {
    cj <- ci + dx; rj <- ri + dy
    ok <- cj >= 0 & cj < nc & rj >= 0 & rj < nr
    cx <- x0 + (cj + 0.5) * res; cy <- y1 - (rj + 0.5) * res
    d2 <- (x - cx)^2 + (y - cy)^2
    k <- ok & d2 <= rad^2
    if (any(k)) acc[[length(acc) + 1]] <- data.table(cel = rj[k] * nc + cj[k] + 1L, w = 1 / pmax(d2[k], 1e-12), wz = z[k] / pmax(d2[k], 1e-12))
  }
  ac <- rbindlist(acc)[, .(w = sum(w), wz = sum(wz)), by = cel]
  a <- rep(NA_real_, nr * nc); a[ac$cel] <- ac$wz / ac$w
  a <- matrix(a, nr, nc, byrow = TRUE)
  bruto <- a
  vacio0 <- is.na(a)
  a <- rellenar_idw(a, win_pdal)
  a <- rellenar_idw(a, win_fill)
  sin <- which(is.na(a), arr.ind = TRUE)
  if (nrow(sin)) {
    don <- which(!is.na(a), arr.ind = TRUE)
    for (k in seq_len(nrow(sin))) { j <- which.min((don[, 1] - sin[k, 1])^2 + (don[, 2] - sin[k, 2])^2); a[sin[k, 1], sin[k, 2]] <- a[don[j, 1], don[j, 2]] }
  }
  list(a = a, bruto = bruto, x0 = x0, y0 = y1, res = res, pct_vacio = 100 * mean(vacio0), nr = nr, nc = nc)
}

escribir_tif <- function(a, x0, y0, res, ruta, epsg, nodata = -9999) {
  suppressPackageStartupMessages(library(terra))
  r <- rast(nrows = nrow(a), ncols = ncol(a), xmin = x0, xmax = x0 + ncol(a) * res, ymax = y0, ymin = y0 - nrow(a) * res, crs = paste0("EPSG:", epsg))
  v <- a; v[is.na(v)] <- nodata
  values(r) <- as.vector(t(v))
  writeRaster(r, ruta, overwrite = TRUE, datatype = "FLT4S", NAflag = nodata)
}

main_05 <- function(cfg = cargar_config()) {
  suppressPackageStartupMessages(library(lidR))
  ent <- file.path(dir_salida(cfg, "04_suelo"), "recorte_suelo.las")
  out <- dir_salida(cfg, "05_dtm")
  d <- cfg$dtm; sc <- cfg$suelo
  con_paso(cfg, "05_dtm", entrada = ent, parametros = c(d, list(banda_excluida_m = sc$excluir_bajo_hilera_m, rechazo_m = sc$rechazo_sobre_calle_m)), fn = function(p) {
    las <- readLAS(ent, select = "*")
    x <- las@data$X; y <- las@data$Y; z <- las@data$Z; cls <- las@data$Classification
    suelo <- cls == 2L
    p$n_antes <- sum(suelo)
    banda <- rep(FALSE, length(x))
    fl <- cargar_filas(cfg, todas = TRUE)
    for (i in seq_len(nrow(fl))) {
      m <- marco_de_linea(sf::st_geometry(fl)[[i]]); st <- marco_st(m, x, y)
      banda <- banda | (abs(st$t) < sc$excluir_bajo_hilera_m & st$s > -1 & st$s < m$largo + 1)
    }
    calle <- suelo & !banda
    p$info("suelo CSF: %d puntos; en la banda bajo la hilera: %d; de calle: %d", sum(suelo), sum(suelo & banda), sum(calle))
    res <- d$resolucion_m
    dm <- dtm_idw(x[calle], y[calle], z[calle], res, d$idw_window, as.integer(round(d$rellenar_huecos_dist_max_m / res)))
    ruta_dtm <- file.path(out, "dtm.tif")
    escribir_tif(dm$a, dm$x0, dm$y0, res, ruta_dtm, cfg$crs$epsg)
    escribir_tif(dm$bruto, dm$x0, dm$y0, res, file.path(out, "dtm_bruto.tif"), cfg$crs$epsg)
    p$extra$celdas_sin_dato_inicial_pct <- round(dm$pct_vacio, 1)
    dtm <- list(a = dm$a, x0 = dm$x0, y0 = dm$y0, res = res)
    zd <- dtm_z(dtm, x, y, fuera_nan = TRUE)
    dz <- z - zd
    en_banda <- suelo & banda & is.finite(dz)
    esc_antes <- if (any(en_banda)) median(dz[en_banda]) else NaN
    recl <- en_banda & dz > sc$rechazo_sobre_calle_m
    cls2 <- cls; cls2[recl] <- 1L
    suelo2 <- cls2 == 2L
    en_banda2 <- suelo2 & banda & is.finite(dz)
    esc_despues <- if (any(en_banda2)) median(dz[en_banda2]) else NaN
    p95_antes <- if (any(en_banda)) percentil(dz[en_banda], 95) else NaN
    las@data$Classification <- cls2
    writeLAS(las, file.path(out, "recorte_suelo_calles.las"))
    qa <- list(puntos_suelo_csf = sum(suelo), puntos_suelo_bajo_hilera = sum(suelo & banda), reclasificados_hierba = sum(recl),
               pct_reclasificados_de_banda = round(100 * sum(recl) / max(1, sum(suelo & banda)), 1),
               escalon_mediano_bajo_hilera_antes_m = round(esc_antes, 3), escalon_p95_antes_m = round(p95_antes, 3),
               escalon_mediano_bajo_hilera_despues_m = round(esc_despues, 3), banda_m = sc$excluir_bajo_hilera_m, rechazo_m = sc$rechazo_sobre_calle_m)
    guardar_json(qa, file.path(out, "qa_suelo.json"))
    p$info("QA suelo: %s", toJSON(qa, auto_unbox = TRUE, digits = 4))
    p$n_despues <- sum(suelo2)
    ext <- c(dm$x0, dm$x0 + dm$nc * res, dm$y0 - dm$nr * res, dm$y0)
    png_matriz(dm$a, file.path(out, "dtm.png"), "DTM (m) desde el suelo de las calles - UTM 30N", extent = ext,
               cmap = "terrain", cbar = "z (m)")
    # QA altura del sensor sobre el DTM
    tray <- cargar_trayectoria(cfg, cache_dir = dir_salida(cfg, "01_qc"))
    pasos <- ventanas_tiempo(cfg, tray, cargar_filas(cfg))$pasos
    qa2 <- lapply(pasos, function(pa) {
      mm <- tray$t >= pa$t0 & tray$t <= pa$t1
      h <- tray$p[mm, 3] - dtm_z(dtm, tray$p[mm, 1], tray$p[mm, 2], fuera_nan = TRUE); h <- h[is.finite(h)]
      list(fila = pa$fila, lado = pa$lado, t_med_m = round(pa$t_med, 2), n = length(h),
           h_p5 = round(percentil(h, 5), 2), h_med = round(median(h), 2), h_p95 = round(percentil(h, 95), 2))
    })
    meds <- vapply(qa2, function(q) q$h_med, 0)
    guardar_json(list(pasos = qa2, rango_esperado = cfg$qc$altura_sensor_sobre_dtm_m, mediana_global = median(meds),
                      dentro_rango = median(meds) >= cfg$qc$altura_sensor_sobre_dtm_m[1] && median(meds) <= cfg$qc$altura_sensor_sobre_dtm_m[2]),
                 file.path(out, "qa_altura_sensor_sobre_dtm.json"))
    p$salida <- ruta_dtm
  })
}
if (sys.nframe() == 0) main_05()
