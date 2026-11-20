# 04_suelo.R: clasificación de suelo ignorando el ruido (clase 7) (réplica de 04_suelo.py).
# Python: PDAL filters.csf (o smrf). R: lidR::classify_ground(csf(...)) con el paquete RCSF (mismo algoritmo CSF de
# Zhang et al. 2016, misma librería C++ de base que PDAL) sobre los puntos no ruido; suelo = clase 2.
# Si metodo = "smrf" el paso se para: lidR no implementa SMRF (en Python se usa PDAL filters.smrf).
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R"))

main_04 <- function(cfg = cargar_config()) {
  suppressPackageStartupMessages(library(lidR))
  ent <- file.path(dir_salida(cfg, "03_filtrado_sor"), "recorte_sor.las")
  out <- dir_salida(cfg, "04_suelo")
  s <- cfg$suelo
  con_paso(cfg, "04_suelo", entrada = ent, parametros = s, fn = function(p) {
    las <- readLAS(ent, select = "*")
    n <- npoints(las); p$n_antes <- n
    noruido <- las@data$Classification != 7L
    sub <- lidR::LAS(las@data[noruido], las@header, crs = sf::st_crs(las), check = FALSE)
    if (s$metodo == "csf") {
      c <- s$csf
      alg <- csf(sloop_smooth = isTRUE(c$slope_smooth), class_threshold = c$umbral, cloth_resolution = c$resolucion,
                 rigidness = c$rigidez, iterations = c$iteraciones, time_step = c$paso_tiempo)
    } else {
      # lidR no implementa SMRF (solo csf, pmf y mcc) y no se sustituye por otro algoritmo para no dar un suelo distinto
      # al de Python con el mismo nombre: se para con un mensaje claro.
      stop(sprintf("suelo.metodo = '%s' no está disponible en R (lidR no tiene SMRF). Usa 'csf' o la versión de Python (PDAL filters.smrf).",
                   s$metodo))
    }
    sub <- classify_ground(sub, alg, last_returns = FALSE)
    cls <- las@data$Classification
    cls[noruido] <- ifelse(sub@data$Classification == 2L, 2L, cls[noruido])
    las@data$Classification <- cls
    p$n_despues <- sum(cls == 2L)
    p$extra$puntos_suelo <- sprintf("%d (%.1f %%)", p$n_despues, 100 * p$n_despues / n)
    sal <- file.path(out, "recorte_suelo.las")
    writeLAS(las, sal)
    p$salida <- sal
  })
}
if (sys.nframe() == 0) main_04()
