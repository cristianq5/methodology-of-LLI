# -*- coding: utf-8 -*-
# =====================================================================================
# comprobar_entradas.R: comprueba, ANTES de lanzar el pipeline, que están todos los datos de un ensayo y que se pueden leer
# (réplica de python/pipeline/comprobar_entradas.py: mismas comprobaciones y mismos niveles OK / AVISO / ERROR).
#
#     Rscript comprobar_entradas.R --ensayo ficha.yaml [--parametros parametros.yaml]
#
# Solo lee (nunca modifica las entradas) y tarda segundos: lee cabeceras y unas pocas líneas o paquetes, no las nubes enteras.
# Revisa la ficha, la nube LAS, cada fichero de la carpeta de brutos del RS30 (PosI, high_freq, rtk.gga, pcap, LAS-1.EP), que
# los tiempos de todos se solapen, el archivo vectorial de hileras y la carpeta de salida.
# Cada comprobación sale como OK, AVISO (se puede lanzar, pero conviene mirarlo) o ERROR (el pipeline fallaría).
# Informe: <salidas.carpeta>_R/logs/comprobar_entradas.json (la carpeta de salida de R, como el resto del flujo R, para no pisar el
# de Python; si la ficha fija salidas.logs, carpeta compartida con Python, el fichero se llama R_comprobar_entradas.json).
# Termina con código 1 si hay algún ERROR. Qué debe contener la carpeta de entradas: docs/ENTRADAS.md.
#
# Diferencias con Python (declaradas):
#  - gps_time del LAS: Python lee el primer trozo de 1 000 000 puntos con laspy. lidR no permite leer "los N primeros puntos"
#    (sus filtros son por atributo, no por índice), así que aquí se leen los MISMOS 1 000 000 primeros registros con readBin a partir
#    de la cabecera (readLASheader: desplazamiento a los puntos, longitud de registro y formato; gps_time en el byte 20 de los
#    formatos 1, 3, 4 y 5 y en el 22 de los formatos 6 a 10). En un LAZ (comprimido) eso no es posible sin descomprimir la nube
#    entera: se comprueba solo la cabecera y se da un AVISO (sin rango de tiempo del LAS no se comprueba su solape).
#  - El filtro del archivo vectorial se escribe en JSON ({"bloque":1,...}) en vez del dict de Python.
# =====================================================================================
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R"))
source(file.path(DIR_SCRIPTS, "calibracion_global.R"))      # leer_elevaciones_ep (y pcap_hesai.R: .tiempo_host, .u32, PCAP_*)

MAGIA_PCAP <- c(0xA1B2C3D4, 0xA1B23C4D)        # pcap clásico (µs / ns), little endian
N_MUESTRA_LAS <- 1000000                       # puntos del primer trozo (como chunk_iterator(1_000_000) de Python)

# ---- informe -------------------------------------------------------------------------
INF <- new.env(); INF$filas <- list()
.add <- function(nivel, que, detalle = "") {
  detalle <- paste(as.character(detalle), collapse = ", ")
  INF$filas[[length(INF$filas) + 1L]] <- list(nivel = nivel, comprobacion = que, detalle = detalle)
  cat(formatC(nivel, width = -5), "  ", formatC(que, width = -42), " ", detalle, "\n", sep = "")
  flush.console()
}
inf_ok <- function(q, d = "") .add("OK", q, d)
inf_aviso <- function(q, d = "") .add("AVISO", q, d)
inf_error <- function(q, d = "") .add("ERROR", q, d)
n_nivel <- function(n) sum(vapply(INF$filas, function(f) f$nivel == n, TRUE))

# ---- ficheros de texto ---------------------------------------------------------------
# Primera línea de datos (tras 'saltar' líneas) y última línea no vacía, sin leer el fichero entero.
primera_ultima_linea <- function(ruta, saltar = 0) {
  con <- file(ruta, "rb"); on.exit(close(con))
  l <- readLines(con, n = saltar + 1, warn = FALSE)
  prim <- if (length(l) > saltar) trimws(iconv(l[saltar + 1], "latin1", "UTF-8")) else ""
  tam <- file.info(ruta)$size
  seek(con, where = max(0, tam - 4096), origin = "start")
  r <- readBin(con, "raw", 4096L)
  r[r == as.raw(0)] <- as.raw(0x20)
  ls_ <- strsplit(iconv(rawToChar(r), "latin1", "UTF-8"), "\r\n|\n|\r")[[1]]
  ls_ <- ls_[nzchar(trimws(ls_))]
  if (!length(ls_)) stop("fichero sin líneas")
  list(prim = prim, ult = trimws(ls_[length(ls_)]))
}
# Valores numéricos del principio de la línea (se para en el primer campo de texto, p. ej. 'Fixed' del PosI).
numeros <- function(linea) {
  out <- numeric(0)
  for (x in strsplit(trimws(gsub(",", " ", linea)), "[[:space:]]+")[[1]]) {
    v <- suppressWarnings(as.numeric(x))
    if (!nzchar(x) || is.na(v)) break
    out <- c(out, v)
  }
  out
}
# Fichero de trayectoria en texto: nº de columnas numéricas y rango de tiempo (primera columna).
comprobar_texto_tiempo <- function(nombre, ruta, saltar, min_cols) {
  r <- tryCatch(primera_ultima_linea(ruta, saltar), error = function(e) e)
  if (inherits(r, "error")) { inf_error(paste0(nombre, ": lectura"), sprintf("%s (%s)", ruta, conditionMessage(r))); return(NULL) }
  na <- numeros(r$prim); nb <- numeros(r$ult)
  if (length(na) < min_cols) {
    inf_error(paste0(nombre, ": columnas"), sprintf("%d columnas numéricas; se esperan al menos %d", length(na), min_cols)); return(NULL)
  }
  t0 <- na[1]; t1 <- if (length(nb)) nb[1] else NA_real_
  if (!isTRUE(t1 > t0)) { inf_error(paste0(nombre, ": tiempo"), sprintf("el tiempo no avanza (%.3f -> %.3f)", t0, t1)); return(NULL) }
  inf_ok(nombre, sprintf("%d columnas; t %.1f -> %.1f (%.0f s)", length(na), t0, t1, t1 - t0))
  c(t0, t1)
}

# ---- pcap ------------------------------------------------------------------------------
# Inicios (offset 0-based) de los registros pcap de una trama del XT32 (862 B, payload que empieza por EE FF) en [off, off + n).
registros_xt32 <- function(con, off, n_bytes) {
  seek(con, where = off, origin = "start")
  r <- readBin(con, "raw", n_bytes)
  L <- length(r)
  k <- which(r == as.raw(0x5e))                       # incl_len y orig_len = 862 (5E 03 00 00 5E 03 00 00), bytes 8-15 de la cabecera
  k <- k[k + 7L <= L]
  pat <- as.raw(c(0x5e, 0x03, 0x00, 0x00, 0x5e, 0x03, 0x00, 0x00))
  for (d in 1:7) k <- k[r[k + d] == pat[d + 1L]]
  i0 <- k - 1L                                        # posición 0-based del patrón en el trozo
  p <- i0 + 8L + PCAP_OFF_ETH                         # payload = registro + 16 B de cabecera + 42 B Ethernet/IP/UDP
  ok <- i0 >= 8L & p + 2L <= L & i0 - 8L + PCAP_TAM_REG <= L
  ok[ok] <- r[p[ok] + 1L] == as.raw(0xee) & r[p[ok] + 2L] == as.raw(0xff)
  off + i0[ok] - 8
}
# Cabecera global del pcap, registros de trama del XT32 al principio y al final y rango de tiempo (cabecera de cada registro).
comprobar_pcap <- function(ruta) {
  res <- tryCatch({
    tam <- file.info(ruta)$size
    con <- file(ruta, "rb"); on.exit(close(con))
    b4 <- readBin(con, "raw", 4L); magia <- .u32(b4, 1L); hex <- toupper(paste(rev(as.character(b4)), collapse = ""))
    if (!magia %in% MAGIA_PCAP) list(err = c("pcap: formato", sprintf("no es un pcap clásico little endian (magia 0x%s)", hex)))
    else {
      ini <- registros_xt32(con, 24, 2^20)
      fin <- registros_xt32(con, max(24, tam - 2^20), 2^20)
      if (!length(ini) || !length(fin)) list(err = c("pcap: paquetes del XT32", "no hay registros de 862 B con payload EE FF en el primer o el último MB"))
      else list(tam = tam, t0 = .tiempo_host(con, ini[1]), t1 = .tiempo_host(con, fin[length(fin)]))
    }
  }, error = function(e) list(err = c("pcap: lectura", conditionMessage(e))))
  if (!is.null(res$err)) { inf_error(res$err[1], res$err[2]); return(NULL) }
  inf_ok("pcap", sprintf("%.2f GB; paquetes del Hesai XT32; t (host) %.1f -> %.1f (%.0f s)", res$tam / 2^30, res$t0, res$t1, res$t1 - res$t0))
  c(res$t0, res$t1)
}

# ---- LAS -------------------------------------------------------------------------------
# gps_time de los primeros n registros, leídos con readBin a partir de la cabecera (no se lee la nube entera).
gps_primeros <- function(ruta, phb, n) {
  fmt <- phb[["Point Data Format ID"]]; lr <- phb[["Point Data Record Length"]]
  npts <- as.numeric(phb[["Number of point records"]])
  nr <- as.integer(min(n, npts))
  if (nr <= 0) return(numeric(0))
  off_g <- if (fmt >= 6) 22L else 20L
  con <- file(ruta, "rb"); on.exit(close(con))
  seek(con, where = phb[["Offset to point data"]], origin = "start")
  r <- readBin(con, "raw", nr * lr)
  nr <- length(r) %/% lr
  idx <- outer(off_g + 0:7, (seq_len(nr) - 1) * lr, "+") + 1
  readBin(r[idx], "double", n = nr, size = 8, endian = "little")
}
comprobar_las <- function(ruta, epsg_ficha) {
  suppressPackageStartupMessages(library(lidR))
  r <- tryCatch({
    h <- readLASheader(ruta)
    list(h = h, phb = h@PHB, epsg = sf::st_crs(h)$epsg)
  }, error = function(e) e)
  if (inherits(r, "error")) { inf_error("LAS: lectura", sprintf("%s (%s)", ruta, conditionMessage(r))); return(list(rango = NULL, epsg = NULL)) }
  phb <- r$phb; fmt <- phb[["Point Data Format ID"]]
  if (!fmt %in% c(1, 3:10)) { inf_error("LAS: gps_time", sprintf("el formato de punto %d no tiene gps_time", fmt)); return(list(rango = NULL, epsg = NULL)) }
  laz <- grepl("\\.laz$", ruta, ignore.case = TRUE)
  gt <- NULL
  if (!laz) {
    gt <- tryCatch(gps_primeros(ruta, phb, N_MUESTRA_LAS), error = function(e) e)
    if (inherits(gt, "error")) { inf_error("LAS: lectura", sprintf("%s (%s)", ruta, conditionMessage(gt))); return(list(rango = NULL, epsg = NULL)) }
  }
  epsg <- r$epsg; if (!is.null(epsg) && is.na(epsg)) epsg <- NULL
  inf_ok("LAS", sprintf("%.0f puntos, formato %d, versión %d.%d, gps_time presente", as.numeric(phb[["Number of point records"]]), fmt,
                        phb[["Version Major"]], phb[["Version Minor"]]))
  if (!is.null(epsg_ficha)) {
    if (!is.null(epsg) && as.integer(epsg) != as.integer(epsg_ficha))
      inf_aviso("LAS: CRS", sprintf("cabecera EPSG:%s, ficha EPSG:%s (manda la ficha; comprueba que es el correcto)", epsg, epsg_ficha))
    else inf_ok("LAS: CRS", sprintf("EPSG:%s (ficha)", epsg_ficha))
  } else if (!is.null(epsg)) inf_ok("LAS: CRS", sprintf("EPSG:%s (cabecera)", epsg))
  else inf_error("LAS: CRS", "la cabecera no trae CRS: indica crs.epsg en la ficha")
  e_trab <- epsg_ficha %||% epsg
  if (laz) {
    inf_aviso("LAS: gps_time", "LAZ comprimido: en R no se lee una muestra sin descomprimir la nube; rango de tiempo del LAS sin comprobar")
    return(list(rango = NULL, epsg = e_trab))
  }
  if (!length(gt)) { inf_error("LAS: puntos", "la nube está vacía"); return(list(rango = NULL, epsg = epsg)) }
  if (all(gt == 0)) { inf_error("LAS: gps_time", "gps_time vale 0 en todos los puntos leídos"); return(list(rango = NULL, epsg = epsg)) }
  list(rango = c(min(gt), max(gt)), epsg = e_trab)
}

# ---- LAS-1.EP y rtk.gga -----------------------------------------------------------------
comprobar_ep <- function(ruta) {
  el <- tryCatch(leer_elevaciones_ep(ruta), error = function(e) e)
  if (inherits(el, "error")) { inf_error("LAS-1.EP", sprintf("no se puede leer la tabla de elevaciones (%s)", conditionMessage(el))); return(invisible()) }
  if (!(all(is.finite(el)) && -30 < min(el) && min(el) < max(el) && max(el) < 30)) {
    inf_error("LAS-1.EP", sprintf("elevaciones fuera de rango: %.2f .. %.2f°", min(el), max(el))); return(invisible())
  }
  inf_ok("LAS-1.EP", sprintf("32 elevaciones de %.2f a %.2f°", min(el), max(el)))
}
comprobar_rtk <- function(ruta) {
  l <- readLines(ruta, warn = FALSE, encoding = "latin1")
  p <- strsplit(l, ",", fixed = TRUE)
  campo <- function(i) vapply(p, function(z) if (length(z) >= i) z[i] else "", "")
  ok <- lengths(p) > 6 & endsWith(campo(1), "GGA") & grepl("^[0-9]+$", campo(7))
  n <- sum(ok); n_fijo <- sum(ok & campo(7) == "4")
  if (n == 0) inf_aviso("rtk.gga", "no hay sentencias GGA (el 01_qc no podrá dar el % de RTK fijo)")
  else inf_ok("rtk.gga", sprintf("%d épocas GGA, %.1f %% con RTK fijo", n, 100 * n_fijo / n))
}

# a, b = c(t0, t1). Fracción de 'a' que cae dentro de 'b'.
solape <- function(nombre, a, b, minimo = 0.5) {
  if (is.null(a) || is.null(b)) return(invisible())
  s <- max(0, min(a[2], b[2]) - max(a[1], b[1]))
  frac <- s / max(a[2] - a[1], 1e-9)
  texto <- sprintf("%.0f %% dentro (%.1f-%.1f frente a %.1f-%.1f)", 100 * frac, a[1], a[2], b[1], b[2])
  if (frac >= minimo) inf_ok(nombre, texto)
  else if (frac > 0) inf_aviso(nombre, texto)
  else inf_error(nombre, paste0("no se solapan en el tiempo: ", texto, ". ¿Son de la misma captura?"))
}

# Espacio libre (GB) en la unidad de 'ruta' (R no lo trae de serie: PowerShell en Windows, df en el resto). NA si no se puede medir.
espacio_libre_gb <- function(ruta) {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      raiz <- substr(normalizePath(ruta, winslash = "\\"), 1, 3)
      v <- system2("powershell", c("-NoProfile", "-Command", shQuote(sprintf("[System.IO.DriveInfo]::new('%s').AvailableFreeSpace", raiz))),
                   stdout = TRUE, stderr = FALSE)
      as.numeric(trimws(v[length(v)])) / 2^30
    } else {
      v <- system2("df", c("-Pk", shQuote(ruta)), stdout = TRUE)
      as.numeric(strsplit(trimws(v[length(v)]), "[[:space:]]+")[[1]][4]) * 1024 / 2^30
    }
  }, error = function(e) NA_real_, warning = function(w) NA_real_)
}
dentro_de <- function(a, b) {             # ¿la ruta a está dentro de b (o es b)?
  a <- tolower(gsub("[\\\\/]+$", "", normalizePath(a, winslash = "/", mustWork = FALSE)))
  b <- tolower(gsub("[\\\\/]+$", "", normalizePath(b, winslash = "/", mustWork = FALSE)))
  a == b || startsWith(a, paste0(b, "/"))
}
rel_windows <- function(carpeta, fichero) {
  r <- file.path(carpeta, basename(fichero))
  if (.Platform$OS.type == "windows") gsub("/", "\\\\", r) else r
}

# ---- principal -------------------------------------------------------------------------
leer_args <- function(args) {
  out <- list(ensayo = if (nzchar(Sys.getenv(VAR_ENTORNO_ENSAYO))) Sys.getenv(VAR_ENTORNO_ENSAYO) else NULL, parametros = NULL)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (grepl("^--(ensayo|parametros)=", a)) { k <- sub("^--([a-z]+)=.*$", "\\1", a); out[[k]] <- sub("^--[a-z]+=", "", a) }
    else if (a %in% c("--ensayo", "--parametros") && i < length(args)) { out[[sub("^--", "", a)]] <- args[i + 1]; i <- i + 1 }
    else if (a %in% c("-h", "--help")) { cat("Uso: Rscript comprobar_entradas.R --ensayo ficha.yaml [--parametros parametros.yaml]\n"); quit(status = 0) }
    i <- i + 1
  }
  out
}

main_comprobar <- function(argv = commandArgs(TRUE)) {
  a <- leer_args(argv)
  cat("Comprobación de entradas del pipeline LLI\n", strrep("-", 78), "\n", sep = "")
  if (is.null(a$ensayo) || !file.exists(a$ensayo)) { inf_error("ficha del ensayo", sprintf("no existe: %s (usa --ensayo)", a$ensayo %||% "None")); return(1L) }
  ruta_f <- normalizePath(a$ensayo)
  base <- dirname(ruta_f)
  e <- tryCatch(yaml::read_yaml(ruta_f) %||% list(), error = function(ex) ex)
  if (inherits(e, "error")) { inf_error("ficha del ensayo", sprintf("YAML no válido (%s)", conditionMessage(e))); return(1L) }
  faltan <- setdiff(c("ensayo", "entradas", "croquis"), names(e))
  if (length(faltan)) { inf_error("ficha: secciones", sprintf("faltan %s (ver config/ensayo_plantilla.yaml)", paste(faltan, collapse = ", "))); return(1L) }
  inf_ok("ficha", ruta_f)
  par <- tryCatch({
    p_ <- cargar_parametros(ruta_abs(base, e$parametros) %||% a$parametros)
    fusionar(p_, e$ajustes %||% list())
  }, error = function(ex) ex)
  if (inherits(par, "error")) { inf_error("parametros.yaml", conditionMessage(par)); return(1L) }
  inf_ok("parametros.yaml + ajustes de la ficha", sprintf("modo de calibración: %s", par$pcap$modo_calibracion %||% "None"))
  ent <- e$entradas

  # ---- nube LAS
  ruta_las <- ruta_abs(base, ent$las)
  rango_las <- NULL; epsg <- NULL
  if (is.null(ruta_las) || !file.exists(ruta_las)) inf_error("LAS", sprintf("no existe: %s", ruta_las %||% "None"))
  else { r <- comprobar_las(ruta_las, e$crs$epsg); rango_las <- r$rango; epsg <- r$epsg }

  # ---- ficheros de la carpeta de brutos
  brutos <- ruta_abs(base, ent$brutos)
  elegidos <- ent$brutos_archivos %||% list()
  if (!is.null(brutos) && !dir.exists(brutos)) { inf_error("carpeta de brutos", sprintf("no existe: %s", brutos)); brutos <- NULL }
  rutas <- list()
  cg <- par$pcap$calibracion_global %||% list()
  necesita_ep <- identical(par$pcap$modo_calibracion, "global_bloque") && identical(cg$elevaciones, "EP")
  for (nombre in names(ARCHIVOS_BRUTOS)) {
    def <- ARCHIVOS_BRUTOS[[nombre]]
    oblig <- isTRUE(def$obligatorio) || (nombre == "calib_ep" && necesita_ep)
    r <- ruta_abs(base, elegidos[[nombre]])
    if (!is.null(r)) {
      if (file.exists(r)) { rutas[[nombre]] <- r; inf_ok(paste0("brutos: ", nombre), paste0(r, " (indicado en la ficha)")) }
      else inf_error(paste0("brutos: ", nombre), sprintf("indicado en la ficha pero no existe: %s", r))
      next
    }
    if (is.null(brutos)) { (if (oblig) inf_error else inf_aviso)(paste0("brutos: ", nombre), "sin carpeta de brutos ni ruta en brutos_archivos"); next }
    patron_txt <- file.path(def$carpeta, def$patron)
    c_ <- sort(list.files(file.path(brutos, def$carpeta), pattern = def$patron, full.names = TRUE, ignore.case = TRUE))
    if (length(c_) == 1) { rutas[[nombre]] <- normalizePath(c_); inf_ok(paste0("brutos: ", nombre), rel_windows(def$carpeta, c_)) }
    else if (!length(c_)) (if (oblig) inf_error else inf_aviso)(paste0("brutos: ", nombre), sprintf("no se encuentra %s%s", patron_txt, if (oblig) "" else " (opcional)"))
    else inf_error(paste0("brutos: ", nombre), sprintf("%d candidatos para %s: %s. Indica uno en entradas.brutos_archivos.%s",
                                                       length(c_), patron_txt, paste(basename(c_), collapse = ", "), nombre))
  }
  t_posi <- if (!is.null(rutas$trayectoria_posi)) comprobar_texto_tiempo("PosI", rutas$trayectoria_posi, 2, 7) else NULL
  t_hf <- if (!is.null(rutas$trayectoria_hf)) comprobar_texto_tiempo("high_freq", rutas$trayectoria_hf, 0, 8) else NULL
  t_pcap <- if (!is.null(rutas$pcap)) comprobar_pcap(rutas$pcap) else NULL
  if (!is.null(rutas$rtk_gga)) comprobar_rtk(rutas$rtk_gga)
  if (!is.null(rutas$calib_ep)) comprobar_ep(rutas$calib_ep)
  else if (!necesita_ep) inf_ok("LAS-1.EP", "no hace falta con esta configuración")

  # ---- tiempos: todo debe ser de la misma captura
  solape("tiempo: LAS (primer trozo) en high_freq", rango_las, t_hf, minimo = 0.99)
  solape("tiempo: high_freq en PosI", t_hf, t_posi)
  solape("tiempo: high_freq en pcap", t_hf, t_pcap)

  # ---- archivo vectorial de hileras
  cq <- fusionar(list(capa = NULL, campo_id = "id_fila", atributos = list(), filtro = list(), normalizar_texto = TRUE,
                      reemplazos = list(), id_compuesto = FALSE), e$croquis)
  ruta_v <- ruta_abs(base, cq$ruta)
  if (is.null(ruta_v) || !file.exists(ruta_v)) inf_error("archivo vectorial", sprintf("no existe: %s", ruta_v %||% "None"))
  else if (is.null(epsg)) inf_error("archivo vectorial", "sin CRS de trabajo (LAS sin CRS y ficha sin crs.epsg): no se puede comprobar")
  else {
    cq$ruta <- ruta_v
    ej <- e$ejecucion_ensayo %||% list()
    cfg <- par
    cfg$croquis <- cq
    cfg$rutas <- list(filas_gpkg = ruta_v); cfg$rutas["filas_capa"] <- list(cq$capa)
    cfg$crs$epsg <- as.integer(epsg)
    cfg$recorte$modo <- ej$modo %||% "prueba"
    cfg$recorte$filas_prueba_id <- as.character(unlist(ej$hileras_prueba %||% list()))
    tryCatch({
      g <- cargar_filas(cfg, todas = TRUE)
      ids <- as.character(g$id_fila)
      inf_ok("archivo vectorial", sprintf("%d hileras tras el filtro %s; ids %s", nrow(g),
             if (length(cq$filtro)) as.character(jsonlite::toJSON(cq$filtro, auto_unbox = TRUE)) else "{}",
             paste0(paste(head(ids, 15), collapse = ", "), if (length(ids) > 15) " ..." else "")))
      largo <- as.numeric(sf::st_length(g))
      if (any(largo < 2)) inf_aviso("archivo vectorial: longitud", sprintf("hileras de menos de 2 m: %s", paste(ids[largo < 2], collapse = ", ")))
      if (cfg$recorte$modo == "prueba") {
        pid <- cfg$recorte$filas_prueba_id
        no <- setdiff(pid, ids)
        if (!length(pid)) inf_error("hileras de prueba", "modo prueba sin ejecucion_ensayo.hileras_prueba")
        else if (length(no)) inf_error("hileras de prueba", sprintf("no están en el archivo vectorial filtrado: %s", paste(no, collapse = ", ")))
        else inf_ok("hileras de prueba", paste(pid, collapse = ", "))
      }
      if (!is.null(ruta_las) && file.exists(ruta_las)) {
        ph <- readLASheader(ruta_las)@PHB
        x0 <- ph[["Min X"]]; y0 <- ph[["Min Y"]]; x1 <- ph[["Max X"]]; y1 <- ph[["Max Y"]]
        bb <- as.numeric(sf::st_bbox(g))            # xmin, ymin, xmax, ymax
        if (bb[1] >= x0 - 1 && bb[3] <= x1 + 1 && bb[2] >= y0 - 1 && bb[4] <= y1 + 1)
          inf_ok("hileras dentro de la nube", "la extensión de las hileras cae dentro de la del LAS")
        else if (bb[3] < x0 || bb[1] > x1 || bb[4] < y0 || bb[2] > y1)
          inf_error("hileras dentro de la nube", "las hileras no tocan la nube: ¿CRS distinto o captura equivocada?")
        else inf_aviso("hileras dentro de la nube", "parte de las hileras queda fuera de la extensión del LAS")
      }
    }, error = function(ex) inf_error("archivo vectorial", conditionMessage(ex)))
  }

  # ---- salida (la de R: <salidas.carpeta>_R, como cargar_config de utils.R)
  sal <- e$salidas %||% list()
  ruta_s <- paste0(sub("[\\\\/]+$", "", ruta_abs(base, sal$carpeta %||% "salidas")), "_R")
  logs_fijos <- ruta_abs(base, sal$logs)
  ruta_logs <- logs_fijos %||% file.path(ruta_s, "logs")
  if (!is.null(brutos) && dentro_de(ruta_s, brutos)) {
    inf_error("carpeta de salida", sprintf("está dentro de la carpeta de brutos: %s", ruta_s)); ruta_logs <- NULL
  } else {
    ok_ <- tryCatch({ dir.create(ruta_logs, recursive = TRUE, showWarnings = FALSE); dir.exists(ruta_logs) }, error = function(ex) FALSE)
    if (!ok_) { inf_error("carpeta de salida", sprintf("no se puede crear %s", ruta_s)); ruta_logs <- NULL }
    else {
      libre <- espacio_libre_gb(ruta_s)
      if (is.na(libre)) inf_aviso("carpeta de salida", sprintf("%s (no se pudo medir el espacio libre)", ruta_s))
      else (if (libre > 20) inf_ok else inf_aviso)("carpeta de salida", sprintf("%s (%.0f GB libres%s)", ruta_s, libre,
                                                   if (libre > 20) "" else "; un bloque completo ocupa decenas de GB"))
    }
  }
  cat(strrep("-", 78), "\n", sep = "")
  n_err <- n_nivel("ERROR"); n_av <- n_nivel("AVISO")
  cat(sprintf("RESULTADO: %d errores, %d avisos. %s\n", n_err, n_av,
              if (n_err == 0) "Se puede lanzar el pipeline." else "Corrige los ERROR antes de lanzar el pipeline."))
  if (!is.null(ruta_logs) && dir.exists(ruta_logs)) {
    nom <- if (is.null(logs_fijos)) "comprobar_entradas.json" else "R_comprobar_entradas.json"   # logs compartidos con Python: no pisar el suyo
    guardar_json(list(ficha = ruta_f, errores = n_err, avisos = n_av, comprobaciones = INF$filas), file.path(ruta_logs, nom))
    cat("informe:", file.path(ruta_logs, nom), "\n")
  }
  if (n_err) 1L else 0L
}

if (sys.nframe() == 0) quit(status = main_comprobar())
