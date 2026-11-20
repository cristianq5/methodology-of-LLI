# 08_segmentacion.R: puntos de cada hilera = buffer del GPKG (hilera$buffer_m) sobre el tramo (réplica de 08_segmentacion.py).
# 1. Eje inicial = línea del GPKG; eje reajustado por PCA con los puntos de copa (HAG >= reajuste_hag_min_m).
#    Longitud: con hilera$extension = "nube" (modo completo) inicio y final salen de la copa (extension_copa en utils.R);
#    la línea del GPKG solo identifica la hilera y sirve de referencia.
# 2. Coordenadas por punto: s, t (+ izquierda), lado = sign(t). 3. Comprobación bilateral.
# Salida: hileras.rds (un data.table por hilera; equivale a hileras.npz), marcos.json y segmentacion.las (dims extra).
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R")); source(file.path(DIR_SCRIPTS, "trayectoria.R")); source(file.path(DIR_SCRIPTS, "geometria.R"))

main_08 <- function(cfg = cargar_config()) {
  suppressPackageStartupMessages(library(lidR))
  out <- dir_salida(cfg, "08_segmentacion")
  ent <- file.path(dir_salida(cfg, "06_normalizacion"), "recorte_hag.las")
  h <- cfg$hilera; u <- cfg$normalizacion$umbral_altura_m
  filas <- cargar_filas(cfg)
  # código numérico de la dimensión LAS 'id_hilera': el id si es entero; si no, la posición en el archivo vectorial filtrado
  todas_ids <- as.character(cargar_filas(cfg, todas = TRUE)$id_fila)
  numericos <- all(grepl("^[0-9]+$", todas_ids)) && all(suppressWarnings(as.numeric(todas_ids)) < 32767)
  codigo <- stats::setNames(if (numericos) as.integer(todas_ids) else seq_along(todas_ids), todas_ids)
  con_paso(cfg, "08_segmentacion", entrada = ent, parametros = h, fn = function(p) {
    las <- readLAS(ent, select = "*")
    d <- las@data
    x <- d$X; y <- d$Y; z <- d$Z; hag <- d$HeightAboveGround; cls <- d$Classification
    p$n_antes <- length(x)
    tray <- cargar_trayectoria(cfg, cache_dir = dir_salida(cfg, "01_qc"))
    marcos <- list(); datos <- list(); resumen <- list()
    idx_tot <- c(); id_h <- c(); s_tot <- c(); t_tot <- c(); lado_tot <- c()
    for (i in seq_len(nrow(filas))) {
      m0 <- marco_de_linea(sf::st_geometry(filas)[[i]])
      sa <- marco_st(m0, x, y)
      lim <- tramo_s(cfg, m0)
      en_buf <- abs(sa$t) <= h$buffer_m & sa$s >= lim[1] & sa$s <= lim[2] & cls != 7L
      m <- m0; giro <- 0; desp <- 0
      if (isTRUE(h$reajustar_eje)) {
        # el eje se ajusta solo con la copa a lo largo de la línea de referencia (no con objetos de la cabecera)
        cop <- en_buf & hag >= h$reajuste_hag_min_m & sa$s >= max(lim[1], 0) & sa$s <= min(lim[2], m0$largo)
        if (h$reajuste_metodo == "pca" && sum(cop) > h$pca_min_puntos) {
          r <- reajustar_eje_pca(m0, x[cop], y[cop]); m <- r$marco; giro <- r$giro; desp <- r$desp
        }
      }
      sl <- marco_st(m, marco_xy(m0, lim[1], 0)[1], marco_xy(m0, lim[1], 0)[2])$s
      sh <- marco_st(m, marco_xy(m0, lim[2], 0)[1], marco_xy(m0, lim[2], 0)[2])$s
      st <- marco_st(m, x, y)
      ext_info <- list(metodo = "gpkg")
      if (extension_desde_nube(cfg)) {
        # longitud real de la hilera = extensión de la copa en la nube; la línea del GPKG solo es referencia
        g <- c(marco_st(m, marco_xy(m0, 0, 0)[1], marco_xy(m0, 0, 0)[2])$s, marco_st(m, marco_xy(m0, m0$largo, 0)[1], marco_xy(m0, m0$largo, 0)[2])$s)
        cop_e <- abs(st$t) <= h$buffer_m & st$s >= sl & st$s <= sh & cls != 7L & hag >= h$extension_hag_min_m
        e <- extension_copa(cfg, st$s[cop_e], c(sl, sh), g)
        sl <- e$s0; sh <- e$s1; ext_info <- e$info
        if (isTRUE(ext_info$toca_limite_busqueda))
          p$info("AVISO hilera %s: la copa llega al límite de la zona de búsqueda (margen_busqueda_m)", filas$id_fila[i])
      }
      sel <- abs(st$t) <= h$buffer_m & st$s >= sl & st$s <= sh & cls != 7L
      ii <- which(sel)
      lado <- as.integer(sign(st$t[ii])); lado[lado == 0L] <- 1L
      clave <- filas$id_fila[i]
      marcos[[clave]] <- list(marco = marco_a_lista(m), marco_gpkg = marco_a_lista(m0), s_rango = c(sl, sh), extension = ext_info,
                              atributos = valores_atributos(cfg, filas, i), codigo_las = codigo[[clave]],
                              giro_deg = giro, desplazamiento_m = desp)
      datos[[clave]] <- data.table(idx = ii, s = st$s[ii], t = st$t[ii], hag = hag[ii], x = x[ii], y = y[ii], z = z[ii], gps_time = d$gpstime[ii], lado = lado)
      cop_i <- hag[ii] >= u
      n_izq <- sum(lado > 0 & cop_i); n_der <- sum(lado < 0 & cop_i)
      pasos <- pasos_sensor(cfg, tray, m, sl, sh)
      lados_tray <- sort(unique(vapply(pasos, function(q) q$lado, 0L)))
      nreb <- as.integer((sh - sl) / 0.1) + 1L
      reb <- floor((st$s[ii][cop_i] - sl) / cfg$metricas$rebanada_m)
      reb <- reb[reb >= 0 & reb < nreb]
      cuenta <- function(mk) { r <- floor((st$s[ii][cop_i][mk] - sl) / cfg$metricas$rebanada_m); tabulate(r[r >= 0 & r < nreb] + 1L, nreb) }
      izq <- cuenta(lado[cop_i] > 0); der <- cuenta(lado[cop_i] < 0)
      ambos <- mean(izq >= cfg$metricas$min_puntos_rebanada / 2 & der >= cfg$metricas$min_puntos_rebanada / 2)
      resumen[[length(resumen) + 1]] <- c(list(id_fila = clave), valores_atributos(cfg, filas, i), list(n_puntos = length(ii), n_copa = sum(cop_i),
                                              n_copa_izq = n_izq, n_copa_der = n_der, rebanadas_con_ambos_lados_pct = round(100 * ambos, 1),
                                              lados_trayectoria = lados_tray, giro_eje_deg = round(giro, 2), desp_eje_m = round(desp, 3),
                                              s_rango = round(c(sl, sh), 2), largo_m = round(sh - sl, 2),
                                              largo_gpkg_m = round(m0$largo, 2), extension = ext_info$metodo,
                                              dif_inicio_m = ext_info$dif_inicio_m, dif_final_m = ext_info$dif_final_m))
      p$info("hilera %s: n=%d copa=%d (izq %d / der %d) ambos_lados=%.1f%% giro=%.2f desp=%.3f", clave, length(ii), sum(cop_i), n_izq, n_der, 100 * ambos, giro, desp)
      idx_tot <- c(idx_tot, ii); id_h <- c(id_h, rep(codigo[[clave]], length(ii))); s_tot <- c(s_tot, st$s[ii]); t_tot <- c(t_tot, st$t[ii]); lado_tot <- c(lado_tot, lado)
    }
    saveRDS(datos, file.path(out, "hileras.rds"))
    guardar_json(list(marcos = marcos, resumen = resumen), file.path(out, "marcos.json"))
    sub <- lidR::LAS(d[idx_tot], las@header, crs = sf::st_crs(las), check = FALSE)
    sub <- add_lasattribute(sub, as.integer(id_h), "id_hilera", "id de hilera")
    sub <- add_lasattribute(sub, as.numeric(s_tot), "s_hilera", "s a lo largo de la hilera (m)")
    sub <- add_lasattribute(sub, as.numeric(t_tot), "t_hilera", "t perpendicular (m)")
    sub <- add_lasattribute(sub, as.integer(lado_tot), "lado", "lado (+1 izq, -1 der)")
    writeLAS(sub, file.path(out, "segmentacion.las"))
    png(file.path(out, "segmentacion_planta.png"), width = 800, height = 800, res = 100)
    plot(NA, xlim = range(unlist(lapply(datos, function(q) q$x))), ylim = range(unlist(lapply(datos, function(q) q$y))), asp = 1, xlab = "x", ylab = "y", main = "08: copa (HAG>=0,3) y ejes")
    for (k in names(datos)) { q <- datos[[k]]; q <- q[hag >= u]; q <- q[seq(1, .N, length.out = min(.N, 30000))]; points(q$x, q$y, pch = ".", col = ifelse(q$t > 0, "red", "blue")) }
    dev.off()
    p$n_despues <- length(idx_tot)
    p$salida <- out
    p$extra$resumen <- vapply(resumen, function(r) sprintf("%s:%d", r$id_fila, r$n_puntos), "")
  })
}
if (sys.nframe() == 0) main_08()
