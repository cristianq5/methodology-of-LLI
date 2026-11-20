# 01_qc.R: control de calidad previo (réplica de 01_qc.py; no recorta ni modifica nada).
# Cabecera LAS, trayectoria (hf -> UTM, ajuste vs PosI), % RTK fijo, cobertura bilateral por hilera, pasos del recorte.
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R")); source(file.path(DIR_SCRIPTS, "trayectoria.R")); source(file.path(DIR_SCRIPTS, "geometria.R"))

pct_rtk_fijo <- function(ruta) {
  l <- readLines(ruta, warn = FALSE, encoding = "latin1")
  l <- l[grepl("^\\$..GGA", l)]
  q <- suppressWarnings(as.integer(vapply(strsplit(l, ",", fixed = TRUE), function(p) if (length(p) > 6) p[7] else NA_character_, "")))
  q <- q[!is.na(q)]
  list(pct = if (length(q)) 100 * mean(q == 4) else NaN, n = length(q))
}

main_01 <- function(cfg = cargar_config()) {
  out <- dir_salida(cfg, "01_qc")
  con_paso(cfg, "01_qc", entrada = cfg$rutas$las_entrada, parametros = cfg$qc, fn = function(p) {
    suppressPackageStartupMessages(library(lidR))
    res <- list()
    h <- readLASheader(cfg$rutas$las_entrada)
    phb <- h@PHB
    res$las <- list(n_puntos = phb$`Number of point records`, mins = c(phb$`Min X`, phb$`Min Y`, phb$`Min Z`),
                    maxs = c(phb$`Max X`, phb$`Max Y`, phb$`Max Z`), version = paste0(phb$`Version Major`, ".", phb$`Version Minor`),
                    formato = phb$`Point Data Format ID`)
    p$n_antes <- res$las$n_puntos
    p$info("LAS: n=%s formato=%s", format(res$las$n_puntos, scientific = FALSE), res$las$formato)
    tray <- cargar_trayectoria(cfg, cache_dir = out)
    dt <- diff(tray$t)
    res$trayectoria <- list(fuente = cfg$trayectoria_modelo$fuente, n = length(tray$t), t0 = tray$t[1], t1 = tail(tray$t, 1),
                            dt_mediana = median(dt), huecos_mayores_tol = sum(dt > cfg$tiempo$tolerancia_hueco_trayectoria_s),
                            ajuste_vs_posi = tray$info)
    if (!is.null(cfg$rutas$rtk_gga)) {
      rt <- pct_rtk_fijo(cfg$rutas$rtk_gga)
      res$rtk <- list(pct_fijo = rt$pct, n_epocas = rt$n, umbral = cfg$qc$min_pct_rtk_fijo, aviso = rt$pct < cfg$qc$min_pct_rtk_fijo)
    } else {
      rt <- list(pct = NA_real_, n = 0L)
      res$rtk <- list(pct_fijo = NULL, aviso = TRUE, motivo = "no hay rtk.gga en los brutos")
    }
    p$info("trayectoria: n=%d dt_med=%.4f rms_vs_posi=%.3f giro_z=%.3f", length(tray$t), median(dt), tray$info$rms_vs_posi_m, tray$info$giro_z_deg)
    p$info("RTK fijo: %.2f %% (%d épocas)", rt$pct, rt$n)
    # cobertura bilateral (todas las hileras del archivo vectorial filtrado)
    filas <- cargar_filas(cfg, todas = TRUE)
    rng <- cfg$qc$distancia_tray_eje_esperada_m
    lst <- list()
    for (i in seq_len(nrow(filas))) {
      m <- marco_de_linea(sf::st_geometry(filas)[[i]]); st <- marco_st(m, tray$p[, 1], tray$p[, 2])
      fila <- c(list(id_fila = filas$id_fila[i]), valores_atributos(cfg, filas, i), list(largo_m = round(m$largo, 2)))
      for (nm in c("izq", "der")) {
        sg <- if (nm == "izq") 1 else -1
        ok <- abs(st$t) >= rng[1] & abs(st$t) <= rng[2] + 1 & sign(st$t) == sg & st$s > 0 & st$s < m$largo
        fila[[paste0("cobertura_", nm, "_pct")]] <- round(100 * length(unique(floor(st$s[ok]))) / ceiling(m$largo), 1)
        fila[[paste0("t_med_", nm)]] <- if (any(ok)) round(median(st$t[ok]), 2) else NA
      }
      fila$bilateral <- fila$cobertura_izq_pct > 50 && fila$cobertura_der_pct > 50
      lst[[i]] <- fila
    }
    df <- rbindlist(lst)
    fwrite(df, file.path(out, "qc_cobertura_bilateral.csv"))
    p$info("cobertura bilateral: %d/%d hileras bilaterales", sum(df$bilateral), nrow(df))
    res$n_filas <- nrow(df); res$filas_bilaterales <- sum(df$bilateral)
    fp <- cargar_filas(cfg)
    vt <- ventanas_tiempo(cfg, tray, fp)
    res$recorte_ventanas_tiempo <- vt$ventanas; res$recorte_pasos <- vt$pasos
    p$info("ventanas de tiempo del recorte: %s", paste(vapply(vt$ventanas, function(v) sprintf("[%.1f, %.1f]", v[1], v[2]), ""), collapse = " "))
    guardar_json(res, file.path(out, "qc_resumen.json"))
    png(file.path(out, "mapa_trayectoria_filas.png"), width = 800, height = 800, res = 100)
    bb <- sf::st_bbox(filas); mg <- cfg$qc$margen_mapa_m %||% 10       # ventana = extensión del archivo vectorial + margen
    mm <- tray$p[, 1] > bb[["xmin"]] - mg & tray$p[, 1] < bb[["xmax"]] + mg & tray$p[, 2] > bb[["ymin"]] - mg & tray$p[, 2] < bb[["ymax"]] + mg
    plot(tray$p[mm, 1], tray$p[mm, 2], pch = ".", asp = 1, col = "blue", xlab = "x", ylab = "y", main = sprintf("QC %s: trayectoria (hf), hileras y recorte", cfg$proyecto$nombre))
    for (i in seq_len(nrow(filas))) { cc <- sf::st_coordinates(sf::st_geometry(filas)[[i]]); lines(cc[, 1], cc[, 2], col = "darkgreen") }
    # recorte = unión de los rectángulos de cada hilera (como poligono_recorte de geometria.py)
    pol <- sf::st_union(sf::st_sfc(lapply(rect_recorte(cfg, filas), function(r) {
      xy <- rbind(marco_xy(r$marco, r$s0, -r$ancho), marco_xy(r$marco, r$s1, -r$ancho), marco_xy(r$marco, r$s1, r$ancho), marco_xy(r$marco, r$s0, r$ancho))
      sf::st_polygon(list(rbind(xy, xy[1, ])))
    })))
    for (g in sf::st_geometry(sf::st_cast(pol, "POLYGON"))) { cc <- sf::st_coordinates(g); lines(cc[, 1], cc[, 2], col = "red", lwd = 1.2) }
    legend("topright", c("trayectoria", "hileras", "recorte"), col = c("blue", "darkgreen", "red"), lty = c(NA, 1, 1), pch = c(16, NA, NA), pt.cex = 0.5, bg = "white")
    dev.off()
    p$n_despues <- p$n_antes
    p$salida <- out
    p$extra$filas_bilaterales <- sprintf("%d/%d", res$filas_bilaterales, res$n_filas)
    res
  })
}
if (sys.nframe() == 0) main_01()
