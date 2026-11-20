# run_pipeline.R: lanza el flujo completo o desde un paso (réplica de run_pipeline.py).
#   Rscript run_pipeline.R --ensayo ficha.yaml                          # 01 .. 11 en el modo de la ficha
#   Rscript run_pipeline.R --ensayo ficha.yaml --desde 08 --hasta 11
#   Rscript run_pipeline.R --ensayo ficha.yaml --solo 09
#   Rscript run_pipeline.R --ensayo ficha.yaml --modo completo --bloque 5          # solo la hilera 5 entera
#   Rscript run_pipeline.R --ensayo ficha.yaml --modo completo --bloque 1,2,3 --procesos 2   # varias hileras en paralelo
#   --calibracion propia|python   (por defecto propia: el ajuste de la extrinseca se hace en R; python = lee las calibraciones calle_*.json de Python)
#   --parametros otro.yaml        (por defecto <repo>/config/parametros.yaml)
#   --recalibrar si               (ignora la caché de calibraciones propias)
# La ficha del ensayo y parametros.yaml son los MISMOS que usa Python (ver config/ensayo_plantilla.yaml).
# Las salidas van a <salidas.carpeta de la ficha>_R (NO se pisa la de Python) y los logs a <salida_R>/logs/R_*.log.
# "bloque" = unidad de trabajo del pipeline (en modo completo, una hilera), no el bloque del ensayo.
a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
DIR_SCRIPTS <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
source(file.path(DIR_SCRIPTS, "utils.R"))

PASOS <- c("01_qc", "02_recorte", "03_filtrado_sor", "04_suelo", "05_dtm", "06_normalizacion", "07_chm", "08_segmentacion", "09_metricas", "10_lli", "11_exportar")
GLOBALES <- c("01_qc", "02_recorte")
FUNC <- c("01_qc" = "main_01", "02_recorte" = "main_02", "03_filtrado_sor" = "main_03", "04_suelo" = "main_04", "05_dtm" = "main_05",
          "06_normalizacion" = "main_06", "07_chm" = "main_07", "08_segmentacion" = "main_08", "09_metricas" = "main_09", "10_lli" = "main_10", "11_exportar" = "main_11")

parse_args <- function(x) {
  o <- list(ensayo = NULL, parametros = NULL, desde = "01", hasta = "11", solo = NULL, modo = NULL, bloque = NULL, procesos = 1L, calibracion = "propia", recalibrar = "no")
  i <- 1
  while (i <= length(x)) {
    k <- sub("^--", "", x[i]); v <- if (i < length(x)) x[i + 1] else NA
    if (k %in% names(o)) { o[[k]] <- v; i <- i + 2 } else i <- i + 1
  }
  o$procesos <- as.integer(o$procesos)
  o
}

cargar_pasos <- function(pasos) for (n in unique(sub(":.*$", "", pasos))) source(file.path(DIR_SCRIPTS, paste0(n, ".R")), chdir = FALSE)

# Ejecuta los pasos por bloque (03..11) de una hilera. Se llama en el proceso principal o en un worker.
trabajo_bloque <- function(cfg, b, pasos) {
  cargar_pasos(pasos)
  cb <- cfg_bloque(cfg, b)
  tiempos <- list(); t0 <- Sys.time()
  for (nom in pasos) {
    t1 <- Sys.time()
    cat(strrep("=", 70), "\nPASO", nom, "| bloque", b %||% "prueba", "\n"); flush.console()
    mod <- sub(":.*$", "", nom); fase <- if (grepl(":", nom)) sub("^.*:", "", nom) else NULL    # "09_metricas:calibrar"
    ok <- tryCatch({ if (is.null(fase)) get(FUNC[[mod]])(cb) else get(FUNC[[mod]])(cb, fase = fase); TRUE }, error = function(e) { cat("ERROR en", nom, "hilera", b, ":", conditionMessage(e), "
"); FALSE })
    if (!ok) return(list(bloque = b, tiempo_s = NA, pasos = tiempos, error = nom))
    tiempos[[nom]] <- round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1)
    cat(sprintf("  -> %s [%s]: %.1f s\n", nom, b, tiempos[[nom]]))
    invisible(gc(FALSE))
  }
  list(bloque = b, tiempo_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), pasos = tiempos)
}

consolidar <- function(cfg, bloques) {
  out <- file.path(dir_salida(cfg, "11_exportar"), "consolidado")
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  for (nom in c("lli_hileras.csv", "lli_tramos.csv", "lli_rebanadas.csv", "perfil_anchura.csv")) {
    l <- lapply(bloques, function(b) { r <- file.path(dir_salida(cfg_bloque(cfg, b), "11_exportar"), nom); if (file.exists(r)) fread(r, colClasses = list(character = "id_fila")) else NULL })
    l <- Filter(Negate(is.null), l)
    if (length(l)) fwrite(rbindlist(l, fill = TRUE), file.path(out, nom))
  }
  cat("consolidado en", out, "\n")
}

resumen_final <- function(cfg, bloques) {
  l <- lapply(bloques, function(b) { r <- file.path(dir_salida(cfg_bloque(cfg, b), "11_exportar"), "lli_hileras.csv")
    if (file.exists(r)) fread(r, colClasses = list(character = "id_fila")) else data.table(id_fila = b, estado_P = "NO PROCESADA (falta lli_hileras.csv)") })
  df <- rbindlist(l, fill = TRUE)
  cat(strrep("=", 70), "\nRESUMEN FINAL POR HILERA\n")
  print(as.data.frame(df[, intersect(c("id_fila", atributos_croquis(cfg), "n_ventanas_validas", "n_rebanadas_con_P", "n_rebanadas", "P_media", "LLI_estricta_media", "LLI_percentil_media", "estado_P"), names(df)), with = FALSE]))
  fwrite(df, file.path(file.path(dir_salida(cfg, "11_exportar"), "consolidado"), "resumen_hileras.csv"))
  df
}

main <- function() {
  a <- parse_args(commandArgs(TRUE))
  cfg <- cargar_config(ensayo = a$ensayo, parametros = a$parametros, modo = a$modo, calibracion = a$calibracion)
  Sys.setenv(LLI_ENSAYO = cfg$.ruta_ensayo)               # la heredan los procesos hijos
  cfg$recalibrar <- identical(a$recalibrar, "si")
  if (!is.null(a$bloque)) cfg$.bloques_sel <- trimws(strsplit(a$bloque, ",")[[1]])
  sel <- if (!is.null(a$solo)) PASOS[startsWith(PASOS, a$solo)] else PASOS[PASOS >= paste0(a$desde, "_") & PASOS < paste0(sprintf("%02d", as.integer(a$hasta) + 1L), "_")]
  bloques <- lista_bloques(cfg)
  t0 <- Sys.time()
  cargar_pasos(intersect(sel, GLOBALES))
  for (nom in intersect(sel, GLOBALES)) { cat(strrep("=", 70), "\nPASO", nom, "\n"); get(FUNC[[nom]])(cfg) }
  por_bloque <- setdiff(sel, GLOBALES)
  resumen <- list()
  ejecutar <- function(pasos, etiqueta = "") {
    if (!length(pasos)) return(invisible(NULL))
    if (a$procesos > 1L && length(bloques) > 1L) {
      cat(sprintf("Ejecución en paralelo%s: %d procesos, %d bloques\n", etiqueta, a$procesos, length(bloques)))
      cl <- parallel::makeCluster(a$procesos)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("DIR_SCRIPTS", "PASOS", "FUNC", "cargar_pasos", "trabajo_bloque"), envir = environment(trabajo_bloque) %||% globalenv())
      parallel::clusterEvalQ(cl, source(file.path(DIR_SCRIPTS, "utils.R")))
      res <- parallel::parLapplyLB(cl, bloques, function(b, cfg, pasos) trabajo_bloque(cfg, b, pasos), cfg = cfg, pasos = pasos)
      for (r in res) resumen[[paste0(as.character(r$bloque), etiqueta)]] <<- r
    } else {
      for (b in bloques) { r <- trabajo_bloque(cfg, b, pasos); resumen[[paste0(as.character(r$bloque %||% "prueba"), etiqueta)]] <<- r }
    }
    errs <- Filter(function(r) !is.null(r$error), resumen)
    if (length(errs)) stop(sprintf("fallaron %d bloques (%s)", length(errs), paste(names(errs), collapse = ", ")))
  }
  if ("09_metricas" %in% por_bloque && identical(cfg$pcap$modo_calibracion %||% "por_calle", "global_bloque")) {
    # Calibración global del bloque: (1) pasos previos + 09 en fase "calibrar" en todas las hileras; (2) ajuste conjunto de un
    # único montaje con todas las calles (calibracion_global.R, "09b"); (3) 09 (P con ese montaje) .. 11.
    ejecutar(c(por_bloque[por_bloque < "09"], "09_metricas:calibrar"), " [fase 1]")
    cat(strrep("=", 70), "\nPASO 09b_calibracion_global (montaje único del bloque)\n")
    t1 <- Sys.time()
    source(file.path(DIR_SCRIPTS, "calibracion_global.R"))
    main_09b(cfg)
    resumen[["09b_calibracion_global"]] <- list(tiempo_s = round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1))
    ejecutar(por_bloque[por_bloque >= "09"], " [fase 2]")
  } else ejecutar(por_bloque)
  if (!identical(bloques, list(NULL)) && "11_exportar" %in% sel) { consolidar(cfg, bloques); resumen_final(cfg, bloques) }
  total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  guardar_json(list(tiempo_total_s = round(total, 1), n_procesos = a$procesos, bloques = resumen, modo = cfg$recorte$modo, pasos = sel),
               file.path(cfg$rutas$logs, sprintf("R_tiempos_%s-%s%s.json", substr(sel[1], 1, 2), substr(sel[length(sel)], 1, 2), if (is.null(a$bloque)) "" else paste0("_bloque_", gsub(",", "-", a$bloque)))))
  cat(sprintf("TOTAL %.1f s\n", total))
}
main()
