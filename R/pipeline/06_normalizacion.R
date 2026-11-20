# 06_normalizacion.R: altura sobre el suelo (HAG) con el DTM de 05 (réplica de 06_normalizacion.py).
# Python: PDAL filters.hag_dem (valor de la celda del DTM, zero_ground = TRUE). R: z - valor de celda del DTM (sin interpolar,
# igual que PDAL); los puntos de suelo (clase 2) llevan HAG = 0. Se conservan TODOS los puntos (también < 0,3 m).
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R"))

main_06 <- function(cfg = cargar_config()) {
  suppressPackageStartupMessages(library(lidR))
  ent <- file.path(dir_salida(cfg, "05_dtm"), "recorte_suelo_calles.las")
  dtm_f <- file.path(dir_salida(cfg, "05_dtm"), "dtm.tif")
  out <- dir_salida(cfg, "06_normalizacion")
  con_paso(cfg, "06_normalizacion", entrada = c(ent, dtm_f), parametros = cfg$normalizacion, fn = function(p) {
    las <- readLAS(ent, select = "*")
    p$n_antes <- npoints(las)
    dtm <- leer_dtm(dtm_f)
    zd <- dtm_z(dtm, las@data$X, las@data$Y, fuera_nan = TRUE)
    hag <- las@data$Z - zd
    hag[las@data$Classification == 2L] <- 0
    p$extra$puntos_fuera_del_DTM <- sum(!is.finite(hag))
    hag[!is.finite(hag)] <- 0
    las <- add_lasattribute(las, hag, "HeightAboveGround", "altura sobre el suelo (m)")
    p$n_despues <- npoints(las)
    u <- cfg$normalizacion$umbral_altura_m
    p$extra$HAG_min_max <- round(range(hag), 3)
    p$extra$`puntos_HAG>=umbral` <- sprintf("%d (%.1f %%)", sum(hag >= u), 100 * mean(hag >= u))
    p$extra$HAG_p50_p99 <- round(c(percentil(hag, 50), percentil(hag, 99)), 3)
    sal <- file.path(out, "recorte_hag.las")
    writeLAS(las, sal)
    p$salida <- sal
  })
}
if (sys.nframe() == 0) main_06()
