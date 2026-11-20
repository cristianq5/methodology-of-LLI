# 03_filtrado_sor.R: Statistical Outlier Removal (réplica de 03_filtrado_sor.py).
# Python: PDAL filters.outlier (k=10, multiplicador 2.5). R: lidR::classify_noise(sor(k, m)).
# Los puntos de ruido se MARCAN con clase 7 (lidR usa 18 internamente; se reasigna a 7), no se borran ni se submuestrea.
# DIFERENCIA: lidR::sor calcula la distancia media a los k vecinos con su propio árbol espacial (3D); el umbral es
# media + m*DE de esas distancias medias (igual que PDAL), pero PDAL cuenta k vecinos distintos del punto y lidR puede
# incluir el propio punto en la consulta, lo que desplaza ligeramente el umbral (documentado en la comparativa).
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R"))

main_03 <- function(cfg = cargar_config()) {
  suppressPackageStartupMessages(library(lidR))
  ent <- file.path(dir_salida(cfg, "02_recorte"), "recorte.las")
  out <- dir_salida(cfg, "03_filtrado_sor")
  s <- cfg$sor
  stopifnot(is.null(s$submuestreo_voxel))
  con_paso(cfg, "03_filtrado_sor", entrada = ent, parametros = s, fn = function(p) {
    set_lidr_threads(max(1L, parallel::detectCores() %/% 3L))
    las <- readLAS(ent, select = "*")
    p$n_antes <- npoints(las)
    las <- classify_noise(las, sor(k = s$mean_k, m = s$multiplicador))
    cl <- las@data$Classification
    ruido <- cl == 18L
    las@data$Classification[ruido] <- 7L
    p$n_despues <- sum(!ruido)
    p$extra$marcados_ruido <- sum(ruido)
    sal <- file.path(out, "recorte_sor.las")
    writeLAS(las, sal)
    p$salida <- sal
  })
}
if (sys.nframe() == 0) main_03()
