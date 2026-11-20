# Pipeline LLI en Python (implementación de referencia)

`pipeline/` contiene los pasos `01_qc.py` … `11_exportar.py`, el lanzador `run_pipeline.py`, los módulos comunes y `tests/` (tests sintéticos con resultados calculados a mano). La descripción del método y de cada paso está en [`../docs/PIPELINE.md`](../docs/PIPELINE.md).

| Módulo | Qué hace |
|---|---|
| `utils.py` | Carga la configuración (parámetros + ficha del ensayo), localiza los brutos del RS30, lee el archivo vectorial de hileras, gestiona los logs por paso, el marco de hilera (s, t) y la escritura de LAS. **Hay que importarlo el primero**: prepara el entorno de las DLL de conda en Windows |
| `trayectoria.py` | Trayectoria post-SLAM llevada al CRS del ensayo (ajuste rígido frente al PosI) e interpolación de pose (lineal + slerp) |
| `geometria.py` | Polígono de recorte por tramo de hilera y pasos del sensor por cada calle (un paso = una calle en un sentido) |
| `pcap_hesai.py` | Decodificador del pcap del Hesai XT32: disparos con su tiempo, azimut, canal y distancia (0 = sin eco) |
| `extrinseca.py` | Calibración escáner ↔ trayectoria por calle, con el giro del cabezal del RS30 |
| `validacion.py` | Validación de la calibración con datos no usados en el ajuste (mediana < 0,15 m y p90 < 0,30 m) |
| `porosidad.py` | Recuento de rayos Bt y Bi por hilera: vía pcap (oficial) y vía rayos-LAS (contraste) |
| `lli_core.py` | Núcleo numérico, sin lectura ni escritura de archivos: MH, MW, envolvente de copa, clasificación de rayos, P y LLI |

## Instalación

```bash
conda env create -f python/environment.yml
conda activate lidar-agrolab
```

**Por qué conda.** PDAL y GDAL son librerías en C++ con muchas dependencias del sistema (PROJ, GEOS, libtiff...). En conda-forge vienen compiladas y con versiones compatibles entre sí. Con pip en Windows no hay ruedas oficiales de PDAL. `requirements.txt` repite la lista para pip, solo como referencia en Linux.

## Ejecución

```bash
cd python/pipeline
python run_pipeline.py --ensayo /ruta/mi_ensayo.yaml                    # modo de la ficha (empieza por "prueba")
python run_pipeline.py --ensayo /ruta/mi_ensayo.yaml --modo completo    # todas las hileras
python run_pipeline.py --ensayo /ruta/mi_ensayo.yaml --desde 08 --hasta 10
#   también: --solo 09, --bloque 2,3 (hileras, en modo completo), --procesos N, --limpiar-cache, --parametros otro.yaml
python tests/test_lli_sintetico.py       # tests sintéticos (sin datos); los demás igual
```

Para lanzar un paso suelto (`python 05_dtm.py`), define antes la variable de entorno `LLI_ENSAYO` con la ruta de la ficha.

## Dependencias y por qué estas

Versiones probadas: Python 3.11.16 en Windows 11 (entorno `lidar-agrolab`). Las versiones mínimas del `environment.yml` son las que se han usado; no se ha probado con versiones anteriores.

| Paquete | Versión probada | Para qué se usa | Por qué este y no otro |
|---|---|---|---|
| **numpy** | 2.4.6 | todo el cálculo vectorial (rebanadas, rayos, residuos) | base de todo lo demás; sin él no hay nada |
| **scipy** | 1.17.1 | `cKDTree` (correspondencias de la calibración), `least_squares` (ajuste robusto de la extrínseca), `Rotation`/`Slerp` (orientación de la trayectoria), `distance_transform_edt` (relleno del DTM) | reúne en un paquete estable el kd-tree, el optimizador robusto (soft-L1) y la interpolación de rotaciones; evita sumar Open3D o PyTorch solo para esto |
| **pandas** | 3.0.5 | tablas por rebanada, hilera y tramo; CSV; lectura del PosI | estándar para tablas; geopandas ya depende de él |
| **pyyaml** | 6.0.3 | lectura de `parametros.yaml` y de la ficha; copia de la configuración efectiva | YAML admite comentarios, así que cada parámetro lleva anotado su origen (`[FICHA]`, `[DATO]`, `[SUP]`); R lo lee igual con `yaml` |
| **laspy** + **lazrs** | 2.7.0 | lectura por trozos de LAS de más de 100 M puntos, escritura de LAS con dimensiones extra; `lazrs` solo para escribir `.laz` | lectura por trozos (memoria acotada) y acceso directo a `gps_time`, imprescindible para emparejar con el pcap; `lazrs` es el backend LAZ en Rust, sin LAStools |
| **pdal** + **python-pdal** | 3.5.3 | SOR (`filters.outlier`), suelo (`filters.csf` / `filters.smrf`), HAG (`filters.hag_dem`), DTM/CHM (`writers.gdal`) | implementaciones de referencia y probadas de CSF, SMRF y SOR; R usa el mismo CSF (librería C++ de Zhang et al. 2016), lo que permite comparar Python y R |
| **gdal** | 3.12.3 | lo usan por debajo PDAL (`writers.gdal`), rasterio y pyogrio | no se importa directamente; se fija en el entorno para que las tres librerías usen la misma |
| **rasterio** | 1.4.4 | lectura y escritura de GeoTIFF (DTM, CHM, LLI), `fillnodata`, `rasterize` | API sencilla sobre GDAL, con transformaciones afines y CRS explícitos |
| **geopandas** + **pyogrio** | 1.1.4 / 0.12.1 | lectura del archivo vectorial (GPKG, SHP...), filtro por atributos, reproyección; escritura de `lli_resultados.gpkg` | acepta cualquier formato vectorial de GDAL; `pyogrio` es su motor de lectura y escritura, más rápido que fiona |
| **shapely** | 2.1.2 | polígonos de recorte, unión de tramos, líneas multiparte del archivo vectorial, geometrías de salida | dependencia de geopandas; la versión 2 es vectorizada |
| **pyproj** | 3.7.2 | lat/lon/h del PosI al CRS proyectado del ensayo | la librería PROJ oficial; admite cualquier EPSG |
| **matplotlib** | 3.11.1 | PNG de control (mapa de QC, planta de la segmentación, perfiles de P y LLI) | backend Agg sin pantalla; los PNG sirven para revisar el modo prueba |
| **pillow** | 12.3.0 | reserva para los PNG de matriz si matplotlib falla | ya viene con matplotlib |
| **psutil** | 7.2.2 | RAM pico del proceso y sus hijos (`logs/ram_pico.json`) | multiplataforma; permite dimensionar `ejecucion.n_procesos` |

**Lo que no se usa:**
- **Open3D**: el SOR y el kd-tree ya los dan PDAL y scipy.
- **PyTorch o GPU**: el trazado de rayos es vectorial en numpy y se hace por lotes.
- **Lectores propios del `.EP` del RS30**: es un formato propietario y no documentado, así que la calibración se estima con los datos.
