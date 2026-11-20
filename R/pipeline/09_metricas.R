# 09_metricas.R: MH, MW y P por rebanada (0,1 m) y por lado de cada hilera (réplica de 09_metricas.py).
#   MH, MW : estrictos (máximo) y con percentil. P: vía pcap (disparos decodificados del .pcap en R con readBin) con la pose de la
#   trayectoria y la extrínseca; vía rayos-LAS como contraste. Bi = rayos que ATRAVIESAN la copa; P = Bi/Bt; L = 1 - P.
#
# CALIBRACIÓN: por defecto "propia" (cfg$calibracion, opción --calibracion de run_pipeline.R): el ajuste de la extrínseca se hace en R
# (extrinseca.R::calibrar), sin leer nada de Python. Con "python" se usan los parámetros ajustados por Python (descrito a continuación).
# (Texto de la versión anterior) DECLARACIÓN: con calibracion = "python" el ajuste NO se replica en R: el AJUSTE de la extrínseca (ICP/RANSAC/mínimos cuadrados no lineales con nudos)
# no se ha reimplementado: sus parámetros se leen del intermedio de Python (salidas\09_metricas\_calles\calle_*.json, o
# salidas\09_metricas\extrinseca.json en el recorte de prueba) y se copian a salidas_R\...\_calibracion_python\ con fecha.
# Lo demás (decodificación del pcap, giro del cabezal, rayos, Bt/Bi, P) y la VALIDACIÓN de la extrínseca (hold-out, QA-22) se calculan en R.
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
for (f in c("utils.R", "trayectoria.R", "geometria.R", "pcap_hesai.R", "extrinseca.R", "lli_core.R", "porosidad.R", "calibracion_global.R")) source(file.path(DIR_SCRIPTS, f))

# Calibración de Python por ventana: devuelve lista de informes (json) emparejados por tiempo con 'lims'
buscar_calibraciones <- function(cfg, lims, destino) {
  base <- file.path(cfg$rutas$salida_python, "09_metricas")
  dir.create(destino, showWarnings = FALSE, recursive = TRUE)
  infos <- list(); fuentes <- character(0)
  if (cfg$recorte$modo == "prueba") {
    e <- fromJSON(file.path(base, "extrinseca.json"), simplifyVector = FALSE)$ventanas
    tv <- t(vapply(e, function(v) as.numeric(unlist(v$ventana_t)), c(0, 0)))
    ficheros <- file.path(base, "extrinseca.json")
    for (w in lims) {
      j <- which.min(abs(tv[, 1] - w[1]) + abs(tv[, 2] - w[2]))
      stopifnot(abs(tv[j, 1] - w[1]) + abs(tv[j, 2] - w[2]) < 2)
      infos[[length(infos) + 1]] <- e[[j]]
    }
    file.copy(ficheros, file.path(destino, "extrinseca.json"), overwrite = TRUE)
    fuentes <- sprintf("%s (modificado %s)", ficheros, format(file.mtime(ficheros), "%Y-%m-%d %H:%M"))
  } else {
    fs <- list.files(file.path(base, "_calles"), pattern = "^calle_[0-9]+_[0-9]+\\.json$", full.names = TRUE)
    tt <- do.call(rbind, lapply(strsplit(sub("\\.json$", "", sub("^calle_", "", basename(fs))), "_"), as.numeric))
    for (w in lims) {
      j <- which.min(abs(tt[, 1] - w[1]) + abs(tt[, 2] - w[2]))
      stopifnot(abs(tt[j, 1] - w[1]) + abs(tt[j, 2] - w[2]) < 4)
      file.copy(fs[j], file.path(destino, basename(fs[j])), overwrite = TRUE)
      infos[[length(infos) + 1]] <- fromJSON(file.path(destino, basename(fs[j])), simplifyVector = FALSE)
      fuentes <- c(fuentes, sprintf("%s (modificado %s)", basename(fs[j]), format(file.mtime(fs[j]), "%Y-%m-%d %H:%M")))
    }
  }
  list(infos = infos, fuentes = fuentes)
}
# Calibracion PROPIA de R de una calle (ventana): cache por calle en salidas_R, carpeta 09_metricas, subcarpeta _calles, fichero calle_<t0>_<t1>.json y un directorio
# .lock para que dos procesos no calibren la misma calle. Arranque: extrínsecas (Er, omega) de las calles ya calibradas por R y VÁLIDAS
# (hold-out); nunca de Python. Devuelve el informe (lista) con holdout_*, ajuste_* y valida (= hold-out).
calibracion_propia <- function(cfg, tray, w, lim, las_t, las_xyz, dir_calles, log) {
  dir.create(dir_calles, showWarnings = FALSE, recursive = TRUE)
  ruta <- file.path(dir_calles, sprintf("calle_%d_%d.json", round(lim[1]), round(lim[2])))
  lock <- paste0(ruta, ".lock")
  leer <- function() if (file.exists(ruta)) tryCatch(fromJSON(ruta, simplifyVector = FALSE), error = function(e) NULL) else NULL
  inf <- if (isTRUE(cfg$recalibrar)) NULL else leer()
  if (!is.null(inf) && is.null(inf$k_ajuste)) { log("la caché no registra los grupos del ajuste (k_ajuste): se recalibra"); inf <- NULL }
  if (!is.null(inf)) { log(sprintf("calibración propia reutilizada (%s)", basename(ruta))); return(inf) }
  if (!dir.create(lock, showWarnings = FALSE)) {
    log("otra hilera calibra esta calle; se espera")
    t_w <- Sys.time()
    while (dir.exists(lock) && as.numeric(difftime(Sys.time(), t_w, units = "secs")) < 4 * 3600) { Sys.sleep(20); inf <- leer(); if (!is.null(inf)) return(inf) }
    inf <- leer(); if (!is.null(inf)) return(inf)
    dir.create(lock, showWarnings = FALSE)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  sem <- list()
  for (f in setdiff(list.files(dir_calles, pattern = "^calle_[0-9]+_[0-9]+[.]json$", full.names = TRUE), ruta)) {
    j <- tryCatch(fromJSON(f, simplifyVector = TRUE), error = function(e) NULL)
    if (!is.null(j) && isTRUE(j$valida) && isTRUE(j$origen == "ajuste propio en R (extrinseca.R)")) sem[[length(sem) + 1]] <- list(rotvec_e = as.numeric(j$rotvec_e_deg) * pi / 180, omega = mean(j$omega))
  }
  log(sprintf("calibrando la calle [%.1f, %.1f] con %d semillas de calles propias válidas", lim[1], lim[2], length(sem)))
  inf <- calibrar_calle(cfg, tray, w, lim, las_t, las_xyz, semillas = sem, log = function(s) log(s))
  guardar_json(inf, ruta)
  inf
}
# Distancia entre hileras (m) para el TRV (réplica de distancia_hileras de 09_metricas.py). metricas$distancia_entre_hileras_m si está
# en la configuración; si es null, se mide en el archivo vectorial: distancia del punto medio de la línea de esta hilera (marco_gpkg) a la
# línea más próxima de OTRA hilera (todas las líneas del fichero, sin el filtro del ensayo, para no saltarse vecinas de otro tratamiento).
# Devuelve list(dist, origen). Las líneas se leen una sola vez por fichero/capa/EPSG (.CACHE_LINEAS).
.CACHE_LINEAS <- new.env()
distancia_hileras <- function(cfg, d, p) {
  v <- cfg$metricas$distancia_entre_hileras_m
  if (!is.null(v) && isTRUE(as.numeric(v) != 0)) return(list(dist = as.numeric(v), origen = "config"))
  suppressPackageStartupMessages(library(sf))
  clave <- paste(cfg$rutas$filas_gpkg, cfg$rutas$filas_capa %||% "", cfg$crs$epsg, sep = "|")
  if (is.null(.CACHE_LINEAS[[clave]])) {
    g <- if (!is.null(cfg$rutas$filas_capa)) st_read(cfg$rutas$filas_gpkg, layer = cfg$rutas$filas_capa, quiet = TRUE) else st_read(cfg$rutas$filas_gpkg, quiet = TRUE)
    ep <- st_crs(g)$epsg
    if (!is.na(st_crs(g)) && !is.null(ep) && !is.na(ep) && ep != cfg$crs$epsg) g <- st_transform(g, cfg$crs$epsg)
    assign(clave, st_geometry(g), envir = .CACHE_LINEAS)
  }
  geo <- .CACHE_LINEAS[[clave]]
  m0 <- marco_de_lista(list(origen = d$marco_gpkg$origen, e = d$marco_gpkg$e, largo = d$marco_gpkg$largo))
  pt <- st_sfc(st_point(as.numeric(marco_xy(m0, 0.5 * m0$largo, 0))), crs = st_crs(geo))
  dist <- as.numeric(st_distance(geo, pt))
  dist <- dist[is.finite(dist) & dist > cfg$metricas$distancia_min_otra_hilera_m]       # fuera la propia línea
  if (!length(dist)) {
    p$info("AVISO: no hay otra hilera en el archivo vectorial: TRV sin calcular (indica metricas.distancia_entre_hileras_m)")
    return(list(dist = NA_real_, origen = "sin_hilera_vecina"))
  }
  list(dist = min(dist), origen = "vectorial")
}

.simp <- function(inf) fromJSON(toJSON(inf, auto_unbox = TRUE, digits = NA), simplifyVector = TRUE)

# fase = "todo" (MH, MW, P) | "calibrar" (solo calibración por calle y grupos; la usa run_pipeline.R en el modo global_bloque
# antes del ajuste conjunto del bloque, calibracion_global.R)
main_09 <- function(cfg = cargar_config(), fase = "todo") {
  suppressPackageStartupMessages(library(lidR))
  out <- dir_salida(cfg, "09_metricas")
  seg_dir <- dir_salida(cfg, "08_segmentacion")
  m_ <- cfg$metricas; pc <- cfg$porosidad; rc <- pc$rayos
  con_paso(cfg, "09_metricas", entrada = seg_dir, parametros = list(metricas = m_, porosidad = pc, pcap_extrinseca = cfg$pcap$extrinseca), fn = function(p) {
    t_ini <- Sys.time()
    datos <- readRDS(file.path(seg_dir, "hileras.rds"))
    info_seg <- cargar_json(file.path(seg_dir, "marcos.json"))
    meta <- fromJSON(file.path(dir_salida(cfg, "02_recorte"), "recorte_meta.json"), simplifyVector = FALSE)
    tray <- cargar_trayectoria(cfg, cache_dir = dir_salida(cfg, "01_qc"))
    dtm <- leer_dtm(file.path(dir_salida(cfg, "05_dtm"), "dtm.tif"))
    u <- cfg$normalizacion$umbral_altura_m
    u06 <- pc$umbral_P_sensibilidad_m
    j06 <- as.integer(round((u06 - u) / m_$franja_m))
    ids <- names(datos)
    hmax <- max(vapply(datos, function(d) max(d$hag), 0))
    nf <- as.integer(ceiling((hmax - u) / m_$franja_m)) + 1L
    # ---------------- MH, MW y envolvente por hilera
    hileras <- list()
    for (k in ids) {
      d <- info_seg$marcos[[k]]
      marco <- marco_de_lista(list(origen = d$marco$origen, e = d$marco$e, largo = d$marco$largo))
      sl <- d$s_rango[1]; sh <- d$s_rango[2]
      n_reb <- as.integer(floor((sh - sl) / m_$rebanada_m + 1e-9))
      dd <- datos[[k]]
      perfil <- perfiles_rebanada(dd$s, dd$t, dd$hag, sl, n_reb, m_$rebanada_m, m_$franja_m, u, m_$min_puntos_franja, m_$min_puntos_rebanada,
                                  m_$mh_percentil, m_$mw_percentil_t, n_franjas = nf, env_percentil_t = pc$envolvente_percentil_t)
      env <- envolvente(perfil, if (pc$envolvente_estadistico == "max") "max" else "percentil")
      dh <- distancia_hileras(cfg, d, p)
      hileras[[k]] <- list(id = k, marco = marco, s_ini = sl, s_fin = sh, perfil = perfil, env = env,
                           dtm_fn = dtm_fn_hilera(dtm, marco), atributos = d$atributos %||% list(),
                           copa = metricas_copa(perfil, dh$dist), dist_hileras = dh$dist, origen_dist = dh$origen)
      p$info("hilera %s: %d rebanadas, MH_max medio %.3f, MW_max medio %.3f; distancia entre hileras %.2f m (%s)", k, n_reb,
             nanmean(perfil$MH_max), nanmean(perfil$MW_max), dh$dist, dh$origen)
    }
    # ---------------- LAS normalizado sin ruido (solo gps_time y xyz)
    las <- readLAS(file.path(dir_salida(cfg, "06_normalizacion"), "recorte_hag.las"), select = "*")
    ok <- las@data$Classification != 7L
    tl <- las@data$gpstime[ok]; xyz <- cbind(las@data$X[ok], las@data$Y[ok], las@data$Z[ok])
    rm(las); invisible(gc(FALSE))
    p$n_antes <- length(tl)
    # ---------------- pcap + extrinseca por ventana: calibración PROPIA de R (por defecto) o parámetros de Python (cfg$calibracion)
    lims <- lapply(meta$ventanas_tiempo, function(v) as.numeric(unlist(v)))
    modo_cal <- cfg$calibracion %||% "propia"
    p$info("calibración de la extrínseca: %s", modo_cal)
    calles_dir <- if (cfg$recorte$modo == "prueba") file.path(out, "_calles") else file.path(dirname(out), "_calles")
    if (modo_cal == "python") {
      cal <- buscar_calibraciones(cfg, lims, file.path(out, "_calibracion_python"))
      for (f in cal$fuentes) p$info("calibración (intermedio de Python, NO ajustada en R): %s", f)
      pars <- lapply(cal$infos, function(inf) par_desde_informe(.simp(inf)))
    } else pars <- list()
    ventanas <- list(); resid <- list()
    ce <- cfg$pcap$extrinseca
    inf_ext <- list()
    glob_ <- modo_global(cfg)
    if (glob_ && modo_cal != "propia") stop("pcap.modo_calibracion = global_bloque exige --calibracion propia")
    for (i in seq_along(lims)) {
      ruta_c <- file.path(calles_dir, sprintf("calle_%d_%d.json", round(lims[[i]][1]), round(lims[[i]][2])))
      if (fase == "calibrar" && file.exists(ruta_c) && !isTRUE(cfg$recalibrar)) {      # calle ya calibrada (y con grupos): nada que hacer
        inf0 <- tryCatch(fromJSON(ruta_c, simplifyVector = TRUE), error = function(e) NULL)
        if (!is.null(inf0) && (!glob_ || grupos_validos(ruta_c, inf0))) { p$info("ventana %d: calle ya calibrada (%s)", i, basename(ruta_c)); next }
      }
      t0 <- Sys.time()
      w <- leer_disparos(cfg$pcap$ruta, lims[[i]][1], lims[[i]][2])
      ventanas[[i]] <- w
      p$info("pcap decodificado [%.1f, %.1f] s: %d bloques en %.0f s", lims[[i]][1], lims[[i]][2], length(w$t), as.numeric(difftime(Sys.time(), t0, units = "secs")))
      mw <- tl >= lims[[i]][1] & tl <= lims[[i]][2]
      if (modo_cal == "propia") {
        inf <- .simp(calibracion_propia(cfg, tray, w, lims[[i]], tl[mw], xyz[mw, , drop = FALSE], calles_dir, function(s) p$info("ventana %d: %s", i, s)))
        if (glob_ && !grupos_validos(ruta_c, inf) &&
            guardar_grupos(cfg, tray, ruta_c, inf, w, lims[[i]], tl[mw], xyz[mw, , drop = FALSE], cfg[[".bloque"]]))
          p$info("ventana %d: grupos de ajuste y validación guardados para la calibración global", i)
        if (fase == "calibrar") { ventanas[i] <- list(NULL); next }
        pars[[i]] <- par_desde_informe(inf)
        inf_ext[[i]] <- list(ventana_t = lims[[i]], calibracion = "propia (R)", residuo_holdout_n = inf$holdout_n, residuo_holdout_mediana_m = inf$holdout_mediana_m,
                             residuo_holdout_p90_m = inf$holdout_p90_m, valida_holdout_R = isTRUE(inf$valida),
                             ajuste_mediana_m = inf$ajuste_mediana_m, ajuste_p90_m = inf$ajuste_p90_m, tiempo_calibracion_s = inf$tiempo_s)
        p$info("ventana %d [%.1f, %.1f]: ajuste propio mediana %.3f p90 %.3f | hold-out mediana %.3f p90 %.3f -> %s", i, lims[[i]][1], lims[[i]][2],
               inf$ajuste_mediana_m, inf$ajuste_p90_m, inf$holdout_mediana_m, inf$holdout_p90_m, if (isTRUE(inf$valida)) "VALIDA" else "NO VALIDA: excluida de P")
        next
      }
      # calibración de Python: validación hold-out (QA-22) en R con puntos LAS de la ventana
      k_aj <- unlist(cal$infos[[i]]$k_ajuste)
      if (is.null(k_aj)) stop("La calibración de Python no registra sus grupos (k_ajuste): recalibra en Python ",
                              "(--limpiar-cache) para validar con datos no usados en el ajuste")
      r <- residuos_holdout(cfg, tray, pars[[i]], w, tl[mw], xyz[mw, , drop = FALSE], semilla = 20260520L + i, excluir = k_aj)
      med <- if (length(r)) median(r) else NA_real_; p90 <- if (length(r)) percentil(r, 90) else NA_real_
      val_aj <- isTRUE(cal$infos[[i]]$valida)
      inf_ext[[i]] <- list(ventana_t = lims[[i]], calibracion = "python", residuo_holdout_n = length(r), residuo_holdout_mediana_m = med, residuo_holdout_p90_m = p90,
                           frac_holdout_menor_0_15 = mean(r < 0.15),
                           valida_holdout_R = isTRUE(med < ce$residuo_mediana_max_m && p90 < ce$residuo_p90_max_m),
                           python_residuo_ajuste_mediana_m = cal$infos[[i]]$residuo_mediana_m, python_residuo_ajuste_p90_m = cal$infos[[i]]$residuo_p90_m,
                           python_valida_ajuste = val_aj)
      p$info("ventana %d [%.1f, %.1f]: hold-out R mediana %.3f m, p90 %.3f m, frac<0,15 %.2f -> %s | (Python, sobre el ajuste: mediana %.3f p90 %.3f -> %s)",
             i, lims[[i]][1], lims[[i]][2], med, p90, mean(r < 0.15), if (inf_ext[[i]]$valida_holdout_R) "VALIDA" else "NO VALIDA: excluida de P",
             cal$infos[[i]]$residuo_mediana_m, cal$infos[[i]]$residuo_p90_m, if (val_aj) "valida" else "no valida")
    }
    if (fase == "calibrar") { p$info("fase 'calibrar' terminada: %d ventanas", length(lims)); return(invisible(NULL)) }
    if (glob_) {                                   # P con el montaje global del bloque y la lectura de orientación asociada
      mont <- cargar_montaje(cfg)
      noms <- vapply(mont$calles, function(z) z$calle, "")
      for (i in seq_along(lims)) {
        nom <- sprintf("calle_%d_%d.json", round(lims[[i]][1]), round(lims[[i]][2]))
        k <- match(nom, noms)
        if (is.na(k)) stop(sprintf("La ventana %d (%s) no está en el montaje global del bloque: repite la calibración global", i, nom))
        c_ <- mont$calles[[k]]
        pars[[i]] <- par_ventana(mont, k, isTRUE(c_$valida))
        inf_ext[[i]] <- list(ventana_t = lims[[i]], calibracion = "global_bloque (R)", calle = nom, residuo_holdout_n = c_$holdout_n,
                             residuo_holdout_mediana_m = c_$holdout_mediana_m, residuo_holdout_p90_m = c_$holdout_p90_m, valida_holdout_R = isTRUE(c_$valida),
                             por_calle_mediana_m = c_$por_calle_mediana_m, por_calle_p90_m = c_$por_calle_p90_m, por_calle_valida = isTRUE(c_$por_calle_valida))
        p$info("ventana %d: %s con el montaje global (hold-out %.3f / %.3f) -> %s", i, nom, c_$holdout_mediana_m, c_$holdout_p90_m,
               if (isTRUE(c_$valida)) "VALIDA" else "NO VALIDA: excluida de P")
      }
      tray <- trayectoria_global(cfg)
    }
    valida_vent <- vapply(inf_ext, function(z) z$valida_holdout_R, TRUE)
    n_val <- sum(valida_vent)
    p$info("criterio de validez (hold-out): mediana < %.2f m Y p90 < %.2f m; ventanas válidas %d/%d", ce$residuo_mediana_max_m, ce$residuo_p90_max_m, n_val, length(pars))
    p$extra$ventanas_validas <- sprintf("%d/%d", n_val, length(pars))
    # ---------------- recuento de rayos por hilera
    vias <- list(); tabs <- list(); resumen_P <- list()
    metodo <- pc$metodo
    filas <- list(); P_ids <- list()
    for (h in hileras) {
      n_reb <- h$perfil$n_reb
      acc_v <- nuevo_acc(n_reb, nf); acc_t <- nuevo_acc(n_reb, nf); acc_r <- nuevo_acc(n_reb, nf)
      if (metodo %in% c("pcap", "ambos")) {
        t0 <- Sys.time()
        rayos_pcap(cfg, tray, pars, ventanas, valida_vent, h, acc_v, acc_t, log = function(s) p$info("hilera %s: %s", h$id, s))
        p$info("vía pcap lista en %.0f s: %s", as.numeric(difftime(Sys.time(), t0, units = "secs")), toJSON(acc_v$est, auto_unbox = TRUE))
      }
      if (metodo %in% c("rayos", "ambos")) {
        t0 <- Sys.time()
        mw <- rep(FALSE, length(tl)); for (l in lims) mw <- mw | (tl >= l[1] & tl <= l[2])
        rayos_las(cfg, tray, lims, tl[mw], xyz[mw, , drop = FALSE], h, acc_r)
        p$info("vía rayos-LAS (solo contraste) lista en %.0f s: %s", as.numeric(difftime(Sys.time(), t0, units = "secs")), toJSON(acc_r$est, auto_unbox = TRUE))
      }
      vias_h <- list(pcap = acc_v, pcap_todas = acc_t, rayos = acc_r)
      P_ <- h$perfil
      sc <- h$s_ini + (seq_len(n_reb) - 0.5) * m_$rebanada_m
      xy <- marco_xy(h$marco, sc, rep(0, n_reb))
      df <- data.table(id_fila = h$id, rebanada = seq_len(n_reb) - 1L, s_ini = h$s_ini + (seq_len(n_reb) - 1) * m_$rebanada_m,
                       s_centro = sc, x_eje = xy[, 1], y_eje = xy[, 2], n_pts = P_$n_pts, MH_max = P_$MH_max, MH_p = P_$MH_p, MW_max = P_$MW_max, MW_p = P_$MW_p)
      for (c_ in names(h$copa)) set(df, j = c_, value = h$copa[[c_]])               # métricas de copa (lli_core.R::metricas_copa)
      df[, dist_hileras_m := h$dist_hileras]; df[, origen_dist_hileras := h$origen_dist]
      for (a_ in names(h$atributos)) set(df, j = a_, value = h$atributos[[a_]])     # atributos del archivo vectorial (variedad, bloque...)
      setcolorder(df, c("id_fila", names(h$atributos)))
      for (via in names(vias_h)) {
        if (via == "rayos" && !(metodo %in% c("rayos", "ambos"))) next
        if (via != "rayos" && !(metodo %in% c("pcap", "ambos"))) next
        cn <- vias_h[[via]]$cont
        for (mm in list(c("suma_rayos", ""), c("media_franjas", "_mf"))) {
          r <- p_por_rebanada(cn$Bt, cn$Bi, rc$min_rayos_celda, pc$min_franjas_fiables, mm[1])
          df[[paste0("P_izq_", via, mm[2])]] <- r$Pl[1, ]; df[[paste0("P_der_", via, mm[2])]] <- r$Pl[2, ]; df[[paste0("P_", via, mm[2])]] <- r$Pm
        }
        r6 <- p_por_rebanada(cn$Bt, cn$Bi, rc$min_rayos_celda, pc$min_franjas_fiables, "suma_rayos", j_min = j06)
        df[[paste0("P_izq_", via, "_h06")]] <- r6$Pl[1, ]; df[[paste0("P_der_", via, "_h06")]] <- r6$Pl[2, ]; df[[paste0("P_", via, "_h06")]] <- r6$Pm
        cv <- cobertura_P(cn$Bt, P_$MH_max, u, m_$franja_m, rc$min_rayos_celda)
        df[[paste0("h_max_P_izq_", via)]] <- cv$h_max[1, ]; df[[paste0("h_max_P_der_", via)]] <- cv$h_max[2, ]
        df[[paste0("cob_MH_izq_", via)]] <- cv$frac[1, ]; df[[paste0("cob_MH_der_", via)]] <- cv$frac[2, ]; df[[paste0("cob_MH_", via)]] <- cv$frac_min
        df[[paste0("Bt_izq_", via)]] <- rowSums(cn$Bt[1, , ]); df[[paste0("Bt_der_", via)]] <- rowSums(cn$Bt[2, , ])
        df[[paste0("Bi_izq_", via)]] <- rowSums(cn$Bi[1, , ]); df[[paste0("Bi_der_", via)]] <- rowSums(cn$Bi[2, , ])
        df[[paste0("franjas_fiables_izq_", via)]] <- rowSums(cn$Bt[1, , ] >= rc$min_rayos_celda); df[[paste0("franjas_fiables_der_", via)]] <- rowSums(cn$Bt[2, , ] >= rc$min_rayos_celda)
      }
      filas[[h$id]] <- df
      resumen_P[[h$id]] <- lapply(names(vias_h), function(via) list(via = via, P_media = nanmean(df[[paste0("P_", via)]]), n_reb_con_P = sum(!is.na(df[[paste0("P_", via)]])), n_reb = n_reb, rayos = vias_h[[via]]$est))
      if (metodo %in% c("pcap", "ambos")) saveRDS(lapply(vias_h, function(z) z$cont), file.path(out, sprintf("conteos_rayos_hilera%s.rds", h$id)))
    }
    tabla <- rbindlist(filas, fill = TRUE)
    via_ppal <- if (metodo %in% c("pcap", "ambos")) "pcap" else "rayos"
    tabla[, P := get(paste0("P_", via_ppal))]; tabla[, P_h06 := get(paste0("P_", via_ppal, "_h06"))]
    for (c_ in c("h_max_P_izq", "h_max_P_der", "cob_MH_izq", "cob_MH_der", "cob_MH")) tabla[[c_]] <- tabla[[paste0(c_, "_", via_ppal)]]
    tabla[, via_P := via_ppal]
    pm <- tabla$P; vm <- ifelse(!is.na(pm), "pcap", "")
    for (alt in list(c("pcap_todas", "pcap_ventana_no_valida"), c("rayos", "rayos_LAS"))) {
      cn <- paste0("P_", alt[1])
      if (cn %in% names(tabla)) { falta <- is.na(pm) & !is.na(tabla[[cn]]); pm[falta] <- tabla[[cn]][falta]; vm[falta] <- alt[2] }
    }
    vm[is.na(pm)] <- "sin_P"          # NaN rotulado (como en Python tras el QA)
    tabla[, P_mejor := pm]; tabla[, via_P_mejor := vm]; tabla[, n_ventanas_validas := n_val]
    mot <- rep("", nrow(tabla))
    ppal <- is.na(tabla$P); sin_mh <- is.na(tabla$MH_max)
    if (via_ppal == "pcap" && n_val == 0) {
      mot[ppal] <- sprintf("extrinseca no valida: 0/%d ventanas cumplen mediana<%.2f y p90<%.2f", length(pars), ce$residuo_mediana_max_m, ce$residuo_p90_max_m)
    } else {
      izq_na <- is.na(tabla[[paste0("P_izq_", via_ppal)]]); der_na <- is.na(tabla[[paste0("P_der_", via_ppal)]])
      mot[ppal & izq_na & !der_na] <- sprintf("lado izquierdo sin franjas fiables (<%d franjas con >=%d rayos)", pc$min_franjas_fiables, rc$min_rayos_celda)
      mot[ppal & der_na & !izq_na] <- sprintf("lado derecho sin franjas fiables (<%d franjas con >=%d rayos)", pc$min_franjas_fiables, rc$min_rayos_celda)
      mot[ppal & izq_na & der_na] <- "ambos lados sin franjas fiables (sin rayos suficientes)"
      if (n_val < length(pars)) mot[ppal] <- paste0(ifelse(nzchar(mot[ppal]), paste0(mot[ppal], " | "), ""), sprintf("%d/%d ventanas excluidas por extrinseca", length(pars) - n_val, length(pars)))
    }
    mot[sin_mh & ppal] <- "rebanada sin puntos de copa (MH/MW NaN)"
    tabla[, P_motivo_nan := mot]
    fwrite(tabla, file.path(out, "metricas_rebanadas.csv"))
    # ---------------- perfil de anchura por altura (hilera entera y tramos de hilera$tramo_m)
    perf <- list()
    for (h in hileras) {
      P_ <- h$perfil
      tramo <- as.integer(floor((seq_len(P_$n_reb) - 1) * m_$rebanada_m / cfg$hilera$tramo_m + 1e-9))
      ph_ <- NULL
      for (nivel in c("hilera", "tramo")) {
        f_ <- perfil_anchura(P_, if (nivel == "tramo") tramo else NULL)
        if (!nrow(f_)) next
        f_[, tramo := if (nivel == "tramo") as.integer(grupo) else NA_integer_]
        f_[, grupo := NULL]; f_[, nivel := nivel]; f_[, id_fila := h$id]
        for (a_ in names(h$atributos)) set(f_, j = a_, value = h$atributos[[a_]])
        setcolorder(f_, c("id_fila", names(h$atributos), "nivel", "tramo"))
        perf[[length(perf) + 1L]] <- f_
        if (nivel == "hilera") ph_ <- f_
      }
      if (!is.null(ph_)) {
        hc <- 0.5 * (ph_$h_ini_m + ph_$h_fin_m)
        png_lineas(list(list(ph_$anchura_media_m, hc, "anchura estricta (franjas ocupadas)"),
                        list(ph_$anchura_media_p_m, hc, "anchura robusta (p2,5-p97,5)"),
                        list(ph_$anchura_media_con_huecos_m, hc, "anchura contando huecos")),
                   file.path(out, sprintf("perfil_anchura_hilera%s.png", h$id)), sprintf("Hilera %s: anchura por altura", h$id),
                   "anchura (m)", "altura sobre el suelo (m)")
      }
    }
    fwrite(rbindlist(perf, fill = TRUE), file.path(out, "perfil_anchura.csv"))
    guardar_json(list(via_principal = via_ppal, por_hilera = resumen_P, extrinseca = inf_ext, calibracion_python = if (modo_cal == "python") cal$fuentes else "no usada (calibración propia de R)",
                      umbral_oficial_m = u, umbral_sensibilidad_m = u06, validacion = "hold-out en R (QA-22)"), file.path(out, "comparacion_P.json"))
    guardar_json(list(ventanas = inf_ext, n_ventanas_validas = n_val), file.path(out, "extrinseca.json"))
    for (k in ids) {
      sub <- tabla[id_fila == k]
      series <- list(list(sub$s_centro, sub$P_pcap, "P pcap"), list(sub$s_centro, sub$P_rayos, "P rayos"), list(sub$s_centro, sub$MH_max / 4, "MH/4"), list(sub$s_centro, sub$MW_max / 4, "MW/4"))
      png_lineas(series, file.path(out, sprintf("perfil_P_hilera%s.png", k)), sprintf("Hilera %s: P por rebanada", k), "s (m)", "P (0-1)")
    }
    p$n_despues <- nrow(tabla)
    p$salida <- file.path(out, "metricas_rebanadas.csv")
    p$extra$tiempo_09_s <- round(as.numeric(difftime(Sys.time(), t_ini, units = "secs")))
    tabla
  })
}
if (sys.nframe() == 0) main_09()
