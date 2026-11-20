# 02_recorte.R: extrae del LAS original (solo lectura) los bloques de trabajo (réplica de 02_recorte.py).
# Modo "prueba": un bloque con las hileras filas_prueba_id (tramo tramo_prueba_s_m, buffer buffer_recorte_m).
# Modo "completo": un bloque por hilera (zona de búsqueda de tramo_s: con hilera$extension = "nube", la línea ampliada
# margen_busqueda_m por cada extremo; la longitud real la fija 08 con la copa). NO se submuestrea.
# Diferencia con Python: Python hace una sola pasada y reparte los 12 bloques; aquí cada hilera lee el LAS con el filtro
# de caja (-keep_xy de LASlib) y se procesa hilera a hilera (RAM acotada).
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R")); source(file.path(DIR_SCRIPTS, "trayectoria.R")); source(file.path(DIR_SCRIPTS, "geometria.R"))

# subconjunto de un LAS por máscara lógica (actualiza la cabecera)
sub_las <- function(las, mask) {
  suppressPackageStartupMessages(library(lidR))
  d <- las@data[mask]
  lidR::LAS(d, las@header, crs = sf::st_crs(las), check = FALSE)
}

# Suelo aproximado: mínimo de z por celda de 1 m; consulta = celda del punto (vecino más próximo si la celda no tiene datos)
suelo_aproximado <- function(x, y, z, celda = 1.0) {
  ix <- floor(x / celda); iy <- floor(y / celda)
  dt <- data.table(ix, iy, z)[, .(zmin = min(z)), by = .(ix, iy)]
  function(xq, yq) {
    jx <- floor(xq / celda); jy <- floor(yq / celda)
    k <- paste(jx, jy); kd <- paste(dt$ix, dt$iy)
    v <- dt$zmin[match(k, kd)]
    f <- which(is.na(v))
    if (length(f)) {
      cx <- (dt$ix + .5) * celda; cy <- (dt$iy + .5) * celda
      for (u in unique(k[f])) {
        i <- which(k == u)[1]
        j <- which.min((cx - xq[i])^2 + (cy - yq[i])^2)
        v[k == u] <- dt$zmin[j]
      }
    }
    v
  }
}

main_02 <- function(cfg = cargar_config()) {
  bloques <- lista_bloques(cfg)
  tray <- cargar_trayectoria(cfg, cache_dir = dir_salida(cfg, "01_qc"))
  res_total <- 0
  con_paso(cfg, "02_recorte", entrada = cfg$rutas$las_entrada,
           parametros = c(cfg$recorte, list(bloques = vapply(bloques, function(b) as.character(b %||% "prueba"), ""))), fn = function(p) {
    phb <- lidR::readLASheader(cfg$rutas$las_entrada)@PHB          # nº de puntos del LAS de entrada (cabecera)
    p$n_antes <- max(c(phb[["Number of point records"]], phb[["Number of points records"]], phb[["Extended number of point records"]], 0), na.rm = TRUE)
    for (b in bloques) {
      cb <- cfg_bloque(cfg, b)
      f <- cargar_filas(cb)
      rects <- rect_recorte(cb, f)
      bb <- bbox_rects(rects)
      out <- dir_salida(cb, "02_recorte")
      las <- leer_las(cfg$rutas$las_entrada, filtro = sprintf("-keep_xy %.4f %.4f %.4f %.4f", bb[1], bb[2], bb[3], bb[4]))
      dd <- las@data
      en_poli <- dentro_rects(rects, dd$X, dd$Y)
      las <- sub_las(las, en_poli)
      dd <- las@data
      n_poly <- nrow(dd)
      hz <- dd$Z - suelo_aproximado(dd$X, dd$Y, dd$Z)(dd$X, dd$Y)
      ok <- hz >= cfg$recorte$margen_z_m[1] & hz <= cfg$recorte$margen_z_m[2]
      las <- sub_las(las, ok)
      dd <- las@data
      suppressPackageStartupMessages(library(lidR))
      writeLAS(las, file.path(out, "recorte.las"))
      vt <- ventanas_tiempo(cb, tray, f)
      meta <- list(bloque = cb[[".bloque"]], n_en_poligono = n_poly, n_final = nrow(dd),
                   bbox = c(min(dd$X), min(dd$Y), max(dd$X), max(dd$Y)), gps_time = range(dd$gpstime),
                   ventanas_tiempo = vt$ventanas, pasos = vt$pasos,
                   filas = lapply(seq_len(nrow(f)), function(i) c(list(id_fila = f$id_fila[i]), valores_atributos(cb, f, i), list(largo = marco_de_linea(sf::st_geometry(f)[[i]])$largo))))
      guardar_json(meta, file.path(out, "recorte_meta.json"))
      res_total <- res_total + nrow(dd)
      p$info("bloque %s: %d puntos en polígono, %d tras margen z, %d ventanas de tiempo", cb[[".bloque"]] %||% "prueba", n_poly, nrow(dd), length(vt$ventanas))
      rm(las, dd); gc(FALSE)
    }
    p$n_despues <- res_total
    p$salida <- dir_salida(cfg, "02_recorte")
  })
}
if (sys.nframe() == 0) main_02()
