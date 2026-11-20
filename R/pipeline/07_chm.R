# 07_chm.R: modelo de altura de copa (máximo de HAG por celda de 0,1 m, sin ruido, HAG >= umbral) (réplica de 07_chm.py).
# Python: PDAL writers.gdal tipo max (radio res*sqrt2). R: binning propio con el mismo radio (origen en el mínimo x / máximo y de los puntos) y GeoTIFF con terra.
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R"))

main_07 <- function(cfg = cargar_config()) {
  suppressPackageStartupMessages({library(lidR); library(terra)})
  ent <- file.path(dir_salida(cfg, "06_normalizacion"), "recorte_hag.las")
  out <- dir_salida(cfg, "07_chm")
  cc <- cfg$chm; u <- cfg$normalizacion$umbral_altura_m
  con_paso(cfg, "07_chm", entrada = ent, parametros = c(cc, list(umbral = u)), fn = function(p) {
    las <- readLAS(ent, select = "*")
    d <- las@data
    d <- d[Classification != 7L & HeightAboveGround >= u]
    p$n_antes <- p$n_despues <- nrow(d)
    res <- cc$resolucion_m
    x0 <- min(d$X); y1 <- max(d$Y)
    nc <- floor((max(d$X) - x0) / res) + 1L; nr <- floor((y1 - min(d$Y)) / res) + 1L
    ci <- floor((d$X - x0) / res); ri <- floor((y1 - d$Y) / res)
    # PDAL asigna cada punto a TODAS las celdas cuyo centro está a <= res*sqrt(2) (radio por defecto del escritor GDAL): el máximo
    # de cada celda incluye así los puntos de las celdas vecinas (dilatación). Se replica ese criterio.
    rad <- res * sqrt(2); acc <- list()
    for (dy in -1:1) for (dx in -1:1) {
      cj <- ci + dx; rj <- ri + dy
      cx <- x0 + (cj + 0.5) * res; cy <- y1 - (rj + 0.5) * res
      k <- cj >= 0 & cj < nc & rj >= 0 & rj < nr & (d$X - cx)^2 + (d$Y - cy)^2 <= rad^2
      if (any(k)) acc[[length(acc) + 1L]] <- data.table(cel = rj[k] * nc + cj[k] + 1L, h = d$HeightAboveGround[k])
    }
    g <- rbindlist(acc)[, .(h = max(h)), by = cel]
    a <- rep(NA_real_, nr * nc); a[g$cel] <- g$h
    a <- matrix(a, nr, nc, byrow = TRUE)
    r <- rast(nrows = nr, ncols = nc, xmin = x0, xmax = x0 + nc * res, ymax = y1, ymin = y1 - nr * res, crs = paste0("EPSG:", cfg$crs$epsg))
    v <- a; v[is.na(v)] <- -9999
    values(r) <- as.vector(t(v))
    sal <- file.path(out, "chm.tif")
    writeRaster(r, sal, overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
    p$extra$CHM_max_m <- max(a, na.rm = TRUE)
    p$extra$CHM_dim <- c(nr, nc)
    png_matriz(a, file.path(out, "chm.png"), "CHM (m, 0,1 m)", extent = c(x0, x0 + nc * res, y1 - nr * res, y1),
               cmap = "viridis", cbar = "HAG máx (m)")
    p$salida <- sal
  })
}
if (sys.nframe() == 0) main_07()
