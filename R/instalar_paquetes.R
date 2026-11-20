# instalar_paquetes.R: instala desde CRAN los paquetes que usan los scripts de R/ (probado con R 4.4.2).
# No reinstala los ya presentes. Versiones de referencia en requirements_R.txt.
#   Rscript R/instalar_paquetes.R
paquetes <- c(
  data.table = "1.18.2.1",  # tablas (todo el flujo)
  yaml       = "2.3.12",    # lectura de parametros.yaml y de la ficha del ensayo
  jsonlite   = "2.0.0",     # JSON de salidas y logs
  lidR       = "4.3.2",     # lectura/escritura LAS y nubes (pasos 01-09)
  RCSF       = "1.0.2",     # algoritmo CSF que usa lidR::csf() en 04_suelo
  sf         = "1.1.0",     # hileras GPKG, trayectoria, exportación
  terra      = "1.9.11",    # rásteres DTM/CHM (05, 07, 11)
  RANN       = "2.6.2",     # kd-tree en extrinseca.R
  minpack.lm = "1.2.4"      # ajuste de la extrínseca (extrinseca.R)
)
repos <- getOption("repos"); if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) repos <- c(CRAN = "https://cloud.r-project.org")
instalados <- rownames(installed.packages())
faltan <- setdiff(names(paquetes), instalados)
if (length(faltan)) { cat("Instalando:", paste(faltan, collapse = ", "), "\n"); install.packages(faltan, repos = repos) } else cat("Todos los paquetes ya están instalados.\n")
for (p in names(paquetes)) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  aviso <- if (is.na(v)) "  <-- NO INSTALADO" else if (package_version(v) < package_version(paquetes[[p]])) "  <-- versión inferior a la de referencia" else ""
  cat(sprintf("%-11s %-10s (ref. %s)%s\n", p, v, paquetes[[p]], aviso))
}
