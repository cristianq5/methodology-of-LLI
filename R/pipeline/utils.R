# -*- coding: utf-8 -*-
# =====================================================================================
# utils.R: utilidades comunes del flujo LLI en R (réplica de python/pipeline/utils.py)
#
# Misma configuración que Python (parametros.yaml + ficha del ensayo). Las salidas van a <salida>_R
# (NO se pisan las de Python) y los logs a <salida>_R/logs/R_NN_nombre[_hilera_X].log.
# Marco de hilera: s (a lo largo del eje), t (perpendicular, + izquierda), h = HAG.
# =====================================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
})
options(lidR.progress = FALSE, warn = 1, digits = 10)
Sys.setlocale("LC_COLLATE", "C")
setDTthreads(0)

# ---- localizar el directorio R/pipeline -------------------------------------------
if (!exists("DIR_SCRIPTS")) {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  DIR_SCRIPTS <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}
# Los scripts viven en <repo>/R/pipeline; los parámetros del método, compartidos con Python, en <repo>/config.
RAIZ_REPO <- dirname(dirname(DIR_SCRIPTS))
RUTA_PARAMETROS <- file.path(RAIZ_REPO, "config", "parametros.yaml")
VAR_ENTORNO_ENSAYO <- "LLI_ENSAYO"
`%||%` <- function(a, b) if (is.null(a)) b else a

# Ficheros de la captura del CHCNAV RS30 dentro de la carpeta de brutos (igual que ARCHIVOS_BRUTOS de utils.py)
ARCHIVOS_BRUTOS <- list(
  trayectoria_posi = list(carpeta = "SLAM_DATA", patron = "\\.PosI$", obligatorio = TRUE),
  trayectoria_hf   = list(carpeta = "POST_SLAM_DATA", patron = "^result_trajectory3d_high_freq\\.txt$", obligatorio = TRUE),
  rtk_gga          = list(carpeta = file.path("GPS", "Rtk"), patron = "^rtk\\.gga$", obligatorio = FALSE),
  pcap             = list(carpeta = file.path("SCAN", "Scanner1"), patron = "\\.pcap$", obligatorio = TRUE),
  calib_ep         = list(carpeta = file.path("SCAN", "Scanner1"), patron = "\\.EP$", obligatorio = FALSE)   # tabla de elevaciones (modo global_bloque)
)

# ---- Configuración y rutas ----------------------------------------------------------
fusionar <- function(base, extra) {                   # diccionarios anidados: 'extra' encima de 'base'
  for (k in names(extra)) {
    base[[k]] <- if (is.list(extra[[k]]) && !is.null(names(extra[[k]])) && is.list(base[[k]])) fusionar(base[[k]], extra[[k]]) else extra[[k]]
  }
  base
}
ruta_abs <- function(base_dir, r) {                   # relativas: desde la carpeta de la ficha del ensayo
  if (is.null(r) || !nzchar(r)) return(NULL)
  r <- path.expand(r)
  if (!grepl("^([A-Za-z]:|[\\\\/])", r)) r <- file.path(base_dir, r)
  normalizePath(r, mustWork = FALSE)
}
buscar_unico <- function(carpeta, def, nombre) {
  c <- sort(list.files(file.path(carpeta, def$carpeta), pattern = def$patron, full.names = TRUE, ignore.case = TRUE))
  if (length(c) == 1) return(normalizePath(c))
  if (!length(c)) {
    if (def$obligatorio) stop(sprintf("No se encuentra '%s' en %s. Indícalo en entradas.brutos_archivos.%s", nombre, file.path(carpeta, def$carpeta), nombre))
    return(NULL)
  }
  stop(sprintf("Hay %d candidatos para '%s': %s. Indica uno en entradas.brutos_archivos.%s", length(c), nombre, paste(basename(c), collapse = ", "), nombre))
}
cargar_parametros <- function(ruta = NULL) yaml::read_yaml(ruta %||% RUTA_PARAMETROS)

# Configuración efectiva = parametros.yaml + ficha del ensayo + 'ajustes' (igual que cargar_config de utils.py).
# Las salidas de R van a <salidas.carpeta>_R (no pisan las de Python); 'salida_python' es la de Python
# (la usa calibracion = "python", que lee las calle_*.json de Python).
cargar_config <- function(ensayo = NULL, parametros = NULL, modo = NULL, calibracion = NULL) {
  ensayo <- ensayo %||% (if (nzchar(Sys.getenv(VAR_ENTORNO_ENSAYO))) Sys.getenv(VAR_ENTORNO_ENSAYO) else NULL)
  if (is.null(ensayo)) stop("Falta la ficha del ensayo: run_pipeline.R --ensayo <ficha.yaml> (plantilla: config/ensayo_plantilla.yaml)")
  ensayo <- normalizePath(ensayo, mustWork = TRUE)
  e <- yaml::read_yaml(ensayo)
  for (sec in c("ensayo", "entradas", "croquis")) if (is.null(e[[sec]])) stop(sprintf("La ficha %s no tiene la sección '%s'", ensayo, sec))
  base <- dirname(ensayo)
  cfg <- fusionar(cargar_parametros(ruta_abs(base, e$parametros) %||% parametros), e$ajustes %||% list())
  ent <- e$entradas
  rutas <- list(las_entrada = ruta_abs(base, ent$las), brutos = ruta_abs(base, ent$brutos))
  for (nombre in names(ARCHIVOS_BRUTOS)) {
    r <- ruta_abs(base, ent$brutos_archivos[[nombre]])
    if (is.null(r)) {
      if (is.null(rutas$brutos)) stop(sprintf("Falta entradas.brutos (o entradas.brutos_archivos.%s)", nombre))
      r <- buscar_unico(rutas$brutos, ARCHIVOS_BRUTOS[[nombre]], nombre)
    }
    rutas[nombre] <- list(r)
  }
  cq <- fusionar(list(capa = NULL, campo_id = "id_fila", atributos = list(), filtro = list(), normalizar_texto = TRUE,
                      reemplazos = list(), id_compuesto = FALSE), e$croquis)
  cq$ruta <- ruta_abs(base, cq$ruta)
  rutas$filas_gpkg <- cq$ruta; rutas["filas_capa"] <- list(cq$capa)
  for (k in c("las_entrada", "filas_gpkg", "trayectoria_posi", "trayectoria_hf", "pcap"))
    if (is.null(rutas[[k]]) || !file.exists(rutas[[k]])) stop(sprintf("No existe la entrada '%s': %s", k, rutas[[k]] %||% "(vacía)"))
  sal <- e$salidas %||% list()
  rutas$salida_python <- ruta_abs(base, sal$carpeta %||% "salidas")
  rutas$salida <- paste0(sub("[\\\\/]+$", "", rutas$salida_python), "_R")
  rutas$logs <- ruta_abs(base, sal$logs) %||% file.path(rutas$salida, "logs")
  cfg$pcap$ruta <- rutas$pcap; rutas$pcap <- NULL
  cfg$rutas <- rutas; cfg$croquis <- cq
  epsg <- e$crs$epsg
  if (is.null(epsg)) {
    suppressPackageStartupMessages(library(lidR))
    epsg <- sf::st_crs(readLASheader(rutas$las_entrada))$epsg
    if (is.null(epsg) || is.na(epsg)) stop("El LAS no trae CRS en la cabecera: indica crs.epsg en la ficha del ensayo")
  }
  cfg$crs$epsg <- as.integer(epsg)
  ej <- e$ejecucion_ensayo %||% list()
  cfg$recorte$modo <- modo %||% ej$modo %||% "prueba"
  cfg$recorte$filas_prueba_id <- as.character(unlist(ej$hileras_prueba %||% list()))
  cfg$recorte$tramo_prueba_s_m <- as.numeric(unlist(ej$tramo_prueba_s_m %||% list(0, 1e9)))
  if (!cfg$recorte$modo %in% c("prueba", "completo")) stop("ejecucion_ensayo.modo debe ser 'prueba' o 'completo'")
  cfg$proyecto$nombre <- as.character(e$ensayo$nombre %||% "ensayo")
  cfg$ensayo <- e$ensayo
  cfg$.ruta_ensayo <- ensayo
  cfg$calibracion <- calibracion %||% "propia"        # "propia" (ajuste de la extrínseca en R) | "python" (lee calle_*.json de Python)
  cfg
}

PASOS_GLOBALES <- c("01_qc")

dir_salida <- function(cfg, nombre) {
  base <- cfg$rutas$salida
  if (cfg$recorte$modo == "prueba") base <- file.path(base, cfg$subdir_prueba %||% "prueba")
  d <- file.path(base, nombre)
  if (!is.null(cfg[[".bloque"]]) && !(nombre %in% PASOS_GLOBALES) && cfg$recorte$modo != "prueba")
    d <- file.path(d, paste0("hilera_", cfg[[".bloque"]]))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

cfg_bloque <- function(cfg, bloque) {
  c <- cfg
  if (!is.null(bloque)) { c$.bloque <- as.character(bloque); c$.filas_bloque <- as.character(bloque) }
  c
}

lista_bloques <- function(cfg) {
  if (cfg$recorte$modo == "prueba") return(list(NULL))
  ids <- cargar_filas(cfg, todas = TRUE)$id_fila
  if (!is.null(cfg$.bloques_sel)) ids <- ids[ids %in% as.character(cfg$.bloques_sel)]
  as.list(ids)
}

# TRUE si la longitud de la hilera se toma de la nube (hilera$extension = "nube") y no de la línea del archivo vectorial
extension_desde_nube <- function(cfg) identical(cfg$hilera$extension %||% "gpkg", "nube") && cfg$recorte$modo != "prueba"

# tramo [s0, s1] de la hilera, en m a lo largo de la línea del archivo vectorial (QA-01; réplica de utils.tramo_s).
# Con hilera$extension = "nube" (modo completo) es la ZONA DE BÚSQUEDA: la línea ampliada margen_busqueda_m por cada extremo.
tramo_s <- function(cfg, marco) {
  ext <- cfg$hilera$recorte_extremos_m
  lo <- ext; hi <- marco$largo - ext
  if (cfg$recorte$modo == "prueba") {
    ab <- cfg$recorte$tramo_prueba_s_m
    return(c(max(ab[1], lo), min(ab[2], hi)))
  }
  if (extension_desde_nube(cfg)) { mg <- cfg$hilera$margen_busqueda_m; return(c(-mg, marco$largo + mg)) }
  c(lo, hi)
}

# Inicio y final reales de la hilera a partir de la copa (réplica de utils.extension_copa de Python).
# Copa: tramo de extension_celda_m con >= max(extension_min_puntos, extension_frac_min x mediana de los tramos de la línea) puntos.
# 1. Dentro de la línea ref = (g0, g1): primer y último tramo con copa (los huecos internos no cortan la hilera).
# 2. Si la copa llega al extremo de la línea: barrido fuera en pasos de extension_paso_m hasta margen_busqueda_m; el límite es
#    el último paso con copa pegada a la hilera. Luego recorte_extremos_m.
extension_copa <- function(cfg, s_copa, lim, ref) {
  h <- cfg$hilera; c_ <- h$extension_celda_m; paso <- h$extension_paso_m
  M <- min(h$margen_busqueda_m, ref[1] - lim[1], lim[2] - ref[2])
  s_ <- sort(as.numeric(s_copa))
  contar <- function(a, b) findInterval(b, s_, left.open = TRUE) - findInterval(a, s_, left.open = TRUE)   # puntos en [a, b)
  g0 <- ref[1]; g1 <- ref[2]
  bordes <- seq(g0, g1 - 1e-9, by = c_)
  cnt <- vapply(bordes, function(a) contar(a, min(a + c_, g1)), 0)
  tipico <- if (any(cnt > 0)) median(cnt[cnt > 0]) else 0
  umbral <- max(h$extension_min_puntos, (h$extension_frac_min %||% 0) * tipico)
  es_copa <- function(a, b) contar(a, b) >= umbral * (b - a) / c_
  occ <- which(vapply(bordes, function(a) es_copa(a, min(a + c_, g1)), TRUE))
  if (!length(occ)) return(list(s0 = g0, s1 = g1, info = list(metodo = "gpkg (sin copa en la nube)")))
  s0 <- bordes[occ[1]]; s1 <- min(bordes[occ[length(occ)]] + c_, g1); eps <- 1e-6
  pasos <- round(paso * seq_len(round(M / paso)), 3)
  bf <- list(); bi <- list()
  if (s1 >= g1 - c_ + eps) for (m in pasos) {
    a <- g1 + m - paso; b <- g1 + m; ok <- es_copa(a, b); bf[[length(bf) + 1]] <- list(m, contar(a, b), ok)
    if (!ok) break
    s1 <- b
  }
  if (s0 <= g0 + c_ - eps) for (m in pasos) {
    a <- g0 - m; b <- g0 - m + paso; ok <- es_copa(a, b); bi[[length(bi) + 1]] <- list(m, contar(a, b), ok)
    if (!ok) break
    s0 <- a
  }
  ext <- h$recorte_extremos_m
  info <- list(metodo = "nube", s_copa = round(c(s0, s1), 2), s_gpkg = round(ref, 2), dif_inicio_m = round(g0 - s0, 2),
               dif_final_m = round(s1 - g1, 2), largo_nube_m = round(s1 - s0, 2), largo_gpkg_m = round(g1 - g0, 2),
               barrido_inicio = bi, barrido_final = bf,
               toca_limite_busqueda = (g0 - s0 >= M - eps) || (s1 - g1 >= M - eps), recorte_extremos_m = ext,
               puntos_tramo_tipico = round(tipico, 1), umbral_puntos_tramo = round(umbral, 1))
  list(s0 = s0 + ext, s1 = s1 - ext, info = info)
}

# ---- JSON ---------------------------------------------------------------------------
guardar_json <- function(obj, ruta) {
  write_json(obj, ruta, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "string", null = "null")
}
cargar_json <- function(ruta) fromJSON(ruta, simplifyVector = TRUE)

# ---- Log por paso -------------------------------------------------------------------
# con_paso(): escribe logs\R_NN_nombre[_hilera_X].log con entrada, parámetros, puntos antes/después,
# tiempo y salida. fn(p) recibe un entorno 'p' con p$info(), p$n_antes, p$n_despues, p$salida, p$extra
con_paso <- function(cfg, nombre, entrada = "", parametros = list(), fn) {
  suf <- if (!is.null(cfg[[".bloque"]])) paste0("_hilera_", cfg[[".bloque"]]) else ""
  dir.create(cfg$rutas$logs, showWarnings = FALSE, recursive = TRUE)
  ruta <- file.path(cfg$rutas$logs, paste0("R_", nombre, suf, ".log"))
  if (cfg$recorte$modo == "prueba") ruta <- file.path(cfg$rutas$logs, paste0("R_", nombre, "_", cfg$subdir_prueba %||% "prueba", ".log"))
  p <- new.env()
  p$n_antes <- NA; p$n_despues <- NA; p$salida <- ""; p$extra <- list()
  escribir <- function(txt) {
    l <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), txt)
    cat(l, "\n", file = ruta, append = TRUE, sep = "")
    cat(sprintf("[R_%s] %s\n", nombre, txt))
    flush.console()
  }
  p$info <- function(...) escribir(sprintf(...))
  p$txt <- function(t) escribir(t)
  cat("", file = ruta)
  t0 <- Sys.time()
  escribir("INICIO")
  escribir(paste0("entrada: ", paste(entrada, collapse = " | ")))
  escribir(paste0("parametros: ", toJSON(parametros, auto_unbox = TRUE, null = "null", digits = NA)))
  res <- tryCatch(fn(p), error = function(e) {
    escribir(paste0("FALLO: ", conditionMessage(e)))
    escribir(paste(capture.output(traceback()), collapse = " "))
    e
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  escribir(sprintf("puntos antes: %s | despues: %s", format(p$n_antes, scientific = FALSE), format(p$n_despues, scientific = FALSE)))
  for (k in names(p$extra)) escribir(sprintf("%s: %s", k, paste(format(p$extra[[k]]), collapse = " ")))
  escribir(paste0("salida: ", p$salida))
  escribir(sprintf("tiempo: %.1f s", dt))
  p$tiempo_s <- dt
  if (inherits(res, "error")) { escribir("FIN CON ERROR"); stop(res) }
  escribir("FIN")
  invisible(res)
}

# ---- Hileras del archivo vectorial del ensayo -------------------------------------------------
# Igual que cargar_filas de utils.py: capa de LÍNEAS (una por hilera), id = croquis.campo_id (o compuesto),
# atributos de texto normalizados y 'reemplazos' en memoria, 'filtro', reproyección al EPSG de trabajo.
atributos_croquis <- function(cfg) setdiff(as.character(unlist(cfg$croquis$atributos)), "id_fila")
valores_atributos <- function(cfg, g, i) {
  a <- atributos_croquis(cfg)
  stats::setNames(lapply(a, function(k) { v <- sf::st_drop_geometry(g)[[k]][i]; if (is.factor(v)) as.character(v) else v }), a)
}
texto_campo <- function(v) {
  if (is.numeric(v)) ifelse(is.finite(v) & v == round(v), format(round(v), scientific = FALSE, trim = TRUE), as.character(v))
  else trimws(as.character(v))
}
cargar_filas <- function(cfg, todas = FALSE) {
  suppressPackageStartupMessages(library(sf))
  cq <- cfg$croquis
  g <- if (!is.null(cfg$rutas$filas_capa)) st_read(cfg$rutas$filas_gpkg, layer = cfg$rutas$filas_capa, quiet = TRUE)
       else st_read(cfg$rutas$filas_gpkg, quiet = TRUE)
  filtro <- cq$filtro %||% list()
  campos <- c(cq$campo_id, atributos_croquis(cfg), names(filtro))
  faltan <- setdiff(campos, names(g))
  if (length(faltan)) stop(sprintf("El archivo vectorial no tiene los campos %s. Disponibles: %s", paste(faltan, collapse = ", "),
                                   paste(setdiff(names(g), attr(g, "sf_column")), collapse = ", ")))
  tipos <- unique(as.character(st_geometry_type(g)))
  if (!all(tipos %in% c("LINESTRING", "MULTILINESTRING"))) stop(sprintf("El archivo vectorial debe ser de líneas; tiene %s", paste(tipos, collapse = ", ")))
  if (any(tipos == "MULTILINESTRING")) {
    g <- st_line_merge(g)
    if (!all(st_geometry_type(g) == "LINESTRING")) stop("Hay hileras multiparte que no se pueden unir en una sola línea")
  }
  ep <- st_crs(g)$epsg
  if (!is.null(ep) && !is.na(ep) && ep != cfg$crs$epsg) g <- st_transform(g, cfg$crs$epsg)
  norm <- if (isTRUE(cq$normalizar_texto %||% TRUE)) tolower else identity
  for (a in unique(c(atributos_croquis(cfg), names(filtro)))) {
    if (is.character(g[[a]]) || is.factor(g[[a]])) g[[a]] <- norm(trimws(as.character(g[[a]])))
    rp <- cq$reemplazos[[a]]
    if (!is.null(rp)) for (k in names(rp)) g[[a]][g[[a]] == k] <- rp[[k]]
  }
  for (a in names(filtro)) {
    vals <- norm(texto_campo(unlist(filtro[[a]])))
    g <- g[norm(texto_campo(g[[a]])) %in% vals, ]
  }
  if (!nrow(g)) stop("Ninguna hilera del archivo vectorial cumple el filtro")
  ids <- texto_campo(g[[cq$campo_id]])
  if (isTRUE(cq$id_compuesto)) ids <- do.call(paste, c(list(ids), lapply(atributos_croquis(cfg), function(a) texto_campo(g[[a]])), sep = "_"))
  g$id_fila <- gsub("[^0-9A-Za-z_.-]", "_", ids)
  if (anyDuplicated(g$id_fila)) stop(sprintf("El id de hilera '%s' se repite (%s): añade un filtro o pon croquis.id_compuesto: true",
                                             cq$campo_id, paste(unique(g$id_fila[duplicated(g$id_fila)]), collapse = ", ")))
  if (!todas) {
    if (!is.null(cfg$.filas_bloque)) g <- g[g$id_fila %in% as.character(cfg$.filas_bloque), ]
    else if (cfg$recorte$modo == "prueba") {
      g <- g[g$id_fila %in% as.character(cfg$recorte$filas_prueba_id), ]
      if (!nrow(g)) stop("Modo prueba: ninguna de las hileras_prueba está en el archivo vectorial filtrado")
    }
  }
  g
}

# ---- Marco de hilera (s, t) ---------------------------------------------------------
marco_nuevo <- function(origen_xy, direccion_xy, largo = NULL) {
  e <- direccion_xy / sqrt(sum(direccion_xy^2))
  list(o = as.numeric(origen_xy), e = e, n = c(-e[2], e[1]), largo = largo)
}
marco_de_linea <- function(geom) {
  cc <- sf::st_coordinates(geom)[, 1:2, drop = FALSE]
  d <- cc[nrow(cc), ] - cc[1, ]
  marco_nuevo(cc[1, ], d, sqrt(sum(d^2)))
}
marco_st <- function(m, x, y) {
  qx <- x - m$o[1]; qy <- y - m$o[2]
  list(s = qx * m$e[1] + qy * m$e[2], t = qx * m$n[1] + qy * m$n[2])
}
marco_xy <- function(m, s, t) cbind(m$o[1] + s * m$e[1] + t * m$n[1], m$o[2] + s * m$e[2] + t * m$n[2])
marco_a_lista <- function(m) list(origen = m$o, e = m$e, largo = m$largo)
marco_de_lista <- function(d) marco_nuevo(d$origen, d$e, d$largo)

# Reajuste del eje por PCA (equivale al SVD de numpy sobre los puntos centrados)
reajustar_eje_pca <- function(marco, x, y) {
  cx <- mean(x); cy <- mean(y)
  X <- cbind(x - cx, y - cy)
  d <- svd(crossprod(X), nu = 0)$v[, 1]               # 1.er eje principal (matriz 2x2: rápido)
  if (sum(d * marco$e) < 0) d <- -d
  o2 <- c(cx, cy) + d * sum((marco$o - c(cx, cy)) * d)
  m2 <- marco_nuevo(o2, d, marco$largo)
  giro <- acos(max(-1, min(1, sum(d * marco$e)))) * 180 / pi
  list(marco = m2, giro = giro, desp = sqrt(sum((o2 - marco$o)^2)))
}

# ---- LAS ----------------------------------------------------------------------------
leer_las <- function(ruta, filtro = "", select = "*") {
  suppressPackageStartupMessages(library(lidR))
  readLAS(ruta, select = select, filter = filtro)
}
epsg_las <- function(las, epsg) { suppressPackageStartupMessages(library(sf)); sf::st_crs(las) <- sf::st_crs(epsg); las }

# ---- Utilidades numéricas -----------------------------------------------------------
nanmean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
nanmedian <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
# Percentil con interpolación lineal (igual que numpy.percentile por defecto = quantile type 7)
percentil <- function(x, p) as.numeric(quantile(x, p / 100, type = 7, names = FALSE))
pct_f <- function(x, d = 3) round(x, d)

# ---- Figuras (PNG sencillos, base R) ------------------------------------------------
# Paletas equivalentes a las de matplotlib que usa Python (utils.png_matriz).
paleta_mpl <- function(nombre, n = 256) {
  if (nombre == "terrain") {   # anclas de matplotlib "terrain"
    anc <- c(0, 0.15, 0.25, 0.50, 0.75, 1)
    cols <- rgb(c(0.2, 0.0, 0.0, 1.0, 0.5, 1.0), c(0.2, 0.6, 0.8, 1.0, 0.36, 1.0), c(0.6, 1.0, 0.4, 0.6, 0.33, 1.0))
    f <- lapply(1:3, function(k) approxfun(anc, col2rgb(cols)[k, ] / 255))
    v <- seq(0, 1, length.out = n)
    return(rgb(f[[1]](v), f[[2]](v), f[[3]](v)))
  }
  hcl.colors(n, "viridis")
}
# Mismo formato que Python: 7 x 6 pulgadas a 110 ppp, aspecto 1:1 (metros iguales en x e y), título, ejes en UTM y
# barra de color con su etiqueta. extent = c(xmin, xmax, ymin, ymax); a: fila 1 = norte (como imshow origin="upper").
png_matriz <- function(a, ruta, titulo = "", extent = NULL, cbar = "", cmap = "viridis", vmin = NULL, vmax = NULL) {
  a[!is.finite(a)] <- NA
  if (is.null(extent)) extent <- c(0, ncol(a), 0, nrow(a))
  rg <- range(a, na.rm = TRUE)
  if (!all(is.finite(rg))) rg <- c(0, 1)
  lo <- if (is.null(vmin)) rg[1] else vmin
  hi <- if (is.null(vmax)) rg[2] else vmax
  if (hi <= lo) hi <- lo + 1e-9
  col <- paleta_mpl(cmap)
  png(ruta, width = 770, height = 660, res = 110)
  layout(matrix(1:2, 1), widths = c(5, 1))
  par(mar = c(3, 4.5, 3, 0.5), mgp = c(2, 0.6, 0), las = 1, cex.axis = 0.85)
  rx <- (extent[2] - extent[1]) / ncol(a); ry <- (extent[4] - extent[3]) / nrow(a)
  image(x = extent[1] + (seq_len(ncol(a)) - 0.5) * rx, y = extent[3] + (seq_len(nrow(a)) - 0.5) * ry,
        z = t(a[nrow(a):1, , drop = FALSE]), zlim = c(lo, hi), col = col, asp = 1, useRaster = TRUE,
        xlab = "", ylab = "", main = titulo, font.main = 1, xlim = extent[1:2], ylim = extent[3:4])
  box()
  par(mar = c(3, 0.5, 3, 4))
  zz <- seq(lo, hi, length.out = length(col))
  image(x = 1, y = zz, z = matrix(zz, 1), col = col, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  axis(4, las = 1, cex.axis = 0.85); box()
  mtext(cbar, side = 4, line = 2.8, las = 0, cex = 0.85)
  dev.off()
}
png_lineas <- function(series, ruta, titulo = "", xlabel = "", ylabel = "") {
  png(ruta, width = 880, height = 500, res = 110)
  par(mar = c(4, 4, 3, 1))
  fin <- function(v) { v <- v[is.finite(v)]; if (length(v)) range(v) else c(0, 1) }   # series todas NaN (hilera sin P): eje por defecto
  yl <- fin(unlist(lapply(series, function(s) s[[2]])))
  xl <- fin(unlist(lapply(series, function(s) s[[1]])))
  plot(NA, xlim = xl, ylim = yl, main = titulo, xlab = xlabel, ylab = ylabel)
  grid()
  for (i in seq_along(series)) lines(series[[i]][[1]], series[[i]][[2]], col = i, lwd = 1.2)
  legend("topright", legend = vapply(series, function(s) s[[3]], ""), col = seq_along(series), lty = 1, cex = 0.7)
  dev.off()
}
copiar_config <- function(cfg, destino) {              # ficha del ensayo tal cual + configuración efectiva
  if (!is.null(cfg$.ruta_ensayo)) file.copy(cfg$.ruta_ensayo, file.path(destino, "ficha_ensayo_usada.yaml"), overwrite = TRUE)
  yaml::write_yaml(cfg[!startsWith(names(cfg), ".")], file.path(destino, "config_efectiva.yaml"))
}

# ---- Raster DTM: lectura rápida como matriz ----------------------------------------
# Devuelve lista(a = matriz [fila, col] con NA, x0 = xmin, y0 = ymax, res)
leer_dtm <- function(ruta) {
  suppressPackageStartupMessages(library(terra))
  r <- terra::rast(ruta)
  a <- terra::as.matrix(r, wide = TRUE)
  a[!is.finite(a) | a == -9999] <- NA
  list(a = a, x0 = terra::xmin(r), y0 = terra::ymax(r), res = terra::res(r)[1])
}
# valor del DTM (celda) en (x,y); fuera de rango o sin dato -> mediana global (zmed)
dtm_z <- function(dtm, x, y, fuera_nan = FALSE) {
  cc <- floor((x - dtm$x0) / dtm$res) + 1L
  rr <- floor((dtm$y0 - y) / dtm$res) + 1L
  nr <- nrow(dtm$a); nc <- ncol(dtm$a)
  ok <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
  v <- rep(NA_real_, length(x))
  v[ok] <- dtm$a[(cc[ok] - 1L) * nr + rr[ok]]
  if (!fuera_nan) {
    if (is.null(dtm$zmed)) dtm$zmed <- median(dtm$a, na.rm = TRUE)
    v[is.na(v)] <- dtm$zmed
  }
  v
}
