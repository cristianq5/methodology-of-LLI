# 11_exportar.R: GPKG (hileras, tramos, rebanadas), CSV (incluido el perfil de anchura por altura), GeoTIFF (DTM, CHM, LLI por rebanada), nube con dimensiones extra y copia de la
# configuración (réplica de 11_exportar.py).
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R"))

main_11 <- function(cfg = cargar_config()) {
  suppressPackageStartupMessages({library(sf); library(terra)})
  out <- dir_salida(cfg, "11_exportar")
  e <- cfg$exportar; epsg <- cfg$crs$epsg
  d10 <- dir_salida(cfg, "10_lli"); d09 <- dir_salida(cfg, "09_metricas"); d08 <- dir_salida(cfg, "08_segmentacion")
  con_paso(cfg, "11_exportar", entrada = d10, parametros = e, fn = function(p) {
    H <- fread(file.path(d10, "lli_hileras.csv"), colClasses = list(character = "id_fila")); T <- fread(file.path(d10, "lli_tramos.csv"), colClasses = list(character = "id_fila")); Rb <- fread(file.path(d10, "lli_rebanadas.csv"), colClasses = list(character = "id_fila"))
    marcos <- cargar_json(file.path(d08, "marcos.json"))$marcos
    p$n_antes <- nrow(Rb)
    mk <- function(k) { d <- marcos[[as.character(k)]]; marco_de_lista(list(origen = d$marco$origen, e = d$marco$e, largo = d$marco$largo)) }
    linea <- function(a, b) st_linestring(rbind(a, b))
    gh <- lapply(seq_len(nrow(H)), function(i) { m <- mk(H$id_fila[i]); s <- marcos[[as.character(H$id_fila[i])]]$s_rango; linea(marco_xy(m, s[1], 0), marco_xy(m, s[2], 0)) })
    GH <- st_sf(as.data.frame(H)[, setdiff(names(H), "id_fila")], id_fila = as.character(H$id_fila), geometry = st_sfc(gh, crs = epsg))
    gt <- lapply(seq_len(nrow(T)), function(i) { m <- mk(T$id_fila[i]); linea(marco_xy(m, T$s_ini[i], 0), marco_xy(m, T$s_fin[i], 0)) })
    GT <- st_sf(as.data.frame(T), geometry = st_sfc(gt, crs = epsg))
    gr <- lapply(seq_len(nrow(Rb)), function(i) {
      m <- mk(Rb$id_fila[i]); w <- if (is.finite(Rb$MW_max[i])) Rb$MW_max[i] / 2 else 0.05
      xy <- rbind(marco_xy(m, Rb$s_ini[i], -w), marco_xy(m, Rb$s_ini[i] + 0.1, -w), marco_xy(m, Rb$s_ini[i] + 0.1, w), marco_xy(m, Rb$s_ini[i], w))
      st_polygon(list(rbind(xy, xy[1, ])))
    })
    GR <- st_sf(as.data.frame(Rb), geometry = st_sfc(gr, crs = epsg))
    ruta_gpkg <- file.path(out, "lli_resultados.gpkg")
    if (file.exists(ruta_gpkg)) file.remove(ruta_gpkg)
    if ("gpkg" %in% e$formatos) {
      st_write(GH, ruta_gpkg, layer = "hileras", quiet = TRUE)
      st_write(GT, ruta_gpkg, layer = "tramos_1_5m", quiet = TRUE)
      st_write(GR, ruta_gpkg, layer = "rebanadas", quiet = TRUE)
    }
    if ("csv" %in% e$formatos) {
      for (nm in c("lli_hileras.csv", "lli_tramos.csv", "lli_rebanadas.csv", "lli_resumen.json")) file.copy(file.path(d10, nm), out, overwrite = TRUE)
      for (nm in c("extrinseca.json", "comparacion_P.json", "metricas_rebanadas.csv", "perfil_anchura.csv")) file.copy(file.path(d09, nm), out, overwrite = TRUE)
    }
    if ("geotiff" %in% e$formatos) {
      file.copy(file.path(dir_salida(cfg, "05_dtm"), "dtm.tif"), file.path(out, "dtm.tif"), overwrite = TRUE)
      file.copy(file.path(dir_salida(cfg, "07_chm"), "chm.tif"), file.path(out, "chm.tif"), overwrite = TRUE)
      chm <- rast(file.path(out, "chm.tif"))
      r <- rast(xmin = xmin(chm), xmax = xmax(chm), ymin = ymin(chm), ymax = ymax(chm), resolution = 0.1, crs = paste0("EPSG:", epsg))
      ok <- is.finite(Rb$LLI_estricta)
      if (any(ok)) {
        v <- vect(GR[ok, ]); ras <- rasterize(v, r, field = "LLI_estricta", background = -9999)
        writeRaster(ras, file.path(out, "lli_estricta_rebanadas.tif"), overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
      }
    }
    if (isTRUE(e$las_con_dimensiones_extra)) file.copy(file.path(d08, "segmentacion.las"), file.path(out, "nube_hileras.las"), overwrite = TRUE)
    if (isTRUE(e$guardar_copia_config)) copiar_config(cfg, out)
    p$n_despues <- nrow(GR)
    p$extra$capas_gpkg <- sprintf("hileras(%d), tramos_1_5m(%d), rebanadas(%d)", nrow(GH), nrow(GT), nrow(GR))
    p$salida <- out
  })
}
if (sys.nframe() == 0) main_11()
