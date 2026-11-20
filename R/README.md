# Pipeline LLI en R (réplica)

Réplica en R del pipeline de Python, para contrastar resultados con otra implementación. Lee **los mismos** ficheros de configuración: `config/parametros.yaml` y la ficha del ensayo. Escribe en `<salidas.carpeta>_R`, así que no pisa las salidas de Python.

`pipeline/` tiene los pasos `01_qc.R` … `11_exportar.R`, el lanzador `run_pipeline.R` y los módulos `utils.R`, `trayectoria.R`, `geometria.R`, `pcap_hesai.R`, `extrinseca.R`, `lli_core.R` y `porosidad.R`. Cada uno se llama igual que su equivalente de Python y hace lo mismo (ver [`../docs/PIPELINE.md`](../docs/PIPELINE.md)). `LLI_scripts_R.Rproj` es el proyecto de RStudio.

## Instalación

R 4.4.2. Desde la raíz del repositorio:

```bash
Rscript R/instalar_paquetes.R      # instala desde CRAN solo los que faltan
```

## Ejecución

```bash
Rscript R/pipeline/run_pipeline.R --ensayo /ruta/mi_ensayo.yaml
Rscript R/pipeline/run_pipeline.R --ensayo /ruta/mi_ensayo.yaml --modo completo --procesos 2
Rscript R/pipeline/run_pipeline.R --ensayo /ruta/mi_ensayo.yaml --desde 08 --hasta 11
#   también: --solo 09, --bloque 1,2,3, --parametros otro.yaml, --recalibrar si
#   --calibracion propia|python : propia = R calibra la extrínseca; python = reutiliza las calle_*.json de Python
```

## Dependencias y por qué estas

| Paquete | Versión probada | Para qué se usa | Por qué este y no otro |
|---|---|---|---|
| **lidR** | 4.3.2 | lectura de LAS con filtros, SOR (`classify_noise(sor())`), suelo (`classify_ground(csf())`), escritura de LAS con atributos extra | estándar en R para LiDAR forestal y agrícola; es el equivalente de laspy + PDAL |
| **RCSF** | 1.0.2 | el algoritmo CSF que llama `lidR::csf()` | es la misma librería C++ de Zhang et al. (2016) que usa PDAL, así que el suelo es comparable con Python; sin ella `csf()` no funciona |
| **sf** | 1.1.0 | archivo vectorial de hileras (GPKG, SHP...), filtro, reproyección, transformación del PosI al CRS del ensayo, GPKG de salida | estándar de vectores en R, sobre GDAL, PROJ y GEOS (las mismas librerías que geopandas) |
| **terra** | 1.9.11 | GeoTIFF del DTM, el CHM y el LLI | sucesor de `raster`, más rápido y con CRS explícito |
| **data.table** | 1.18.2.1 | tablas por rebanada, hilera y tramo; lectura del pcap por bloques; CSV | rápido y con poca memoria en tablas de millones de filas (rayos) |
| **yaml** | 2.3.12 | lectura de `parametros.yaml` y de la ficha | mismo formato que Python, así que hay una sola configuración para las dos implementaciones |
| **jsonlite** | 2.0.0 | JSON de resultados, de las calibraciones y de los logs | lee también las `calle_*.json` de Python (`--calibracion python`) |
| **RANN** | 2.6.2 | kd-tree (`nn2`) de las correspondencias en la calibración | equivalente directo de `scipy.spatial.cKDTree` |
| **minpack.lm** | 1.2.4 | ajuste por mínimos cuadrados de la extrínseca (`nls.lm`) | Levenberg-Marquardt robusto; R base no tiene un equivalente de `least_squares` de scipy |
| parallel | (base R) | hileras en paralelo (`--procesos`) | viene con R; funciona en Windows (clúster de sockets) |

**Diferencias deliberadas con Python:**
- **Optimizador:** Levenberg-Marquardt en R frente a `trf` en Python.
- **Generador aleatorio:** el de R, así que los grupos de tiempo elegidos al azar no son los mismos que en Python.
- **Envolvente rellena:** solo está en Python. En R `porosidad.envolvente_relleno` se ignora; el valor por defecto en las dos implementaciones es `false`.
- **Suelo SMRF:** solo está en Python (PDAL `filters.smrf`). lidR no implementa SMRF, así que en R `suelo.metodo: "smrf"` para el paso con un mensaje; usa `"csf"`.

Por estas diferencias, R y Python dan resultados parecidos, pero no idénticos bit a bit.
