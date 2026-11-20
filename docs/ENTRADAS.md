# Datos de entrada de un ensayo

Cada **captura** del CHCNAV RS30 es un ensayo para el pipeline: una nube LAS, la carpeta de brutos de esa misma captura, un
archivo vectorial con las hileras y una ficha (`config/ensayo_plantilla.yaml`). Si un ensayo se capturó en varias pasadas
(por ejemplo, un bloque por captura), cada una lleva su ficha y su carpeta de salida.

Los datos de entrada se abren en **solo lectura**: el pipeline nunca los modifica, renombra ni mueve.

## Organización recomendada

```
mi_ensayo/
├── mi_ensayo.yaml                     ficha del ensayo (copia de config/ensayo_plantilla.yaml)
├── nube.las                           nube exportada del RS30 (LAS o LAZ)
├── hileras.gpkg                       líneas de las hileras
└── brutos/                            carpeta de la captura TAL CUAL la deja el RS30 (sin renombrar ni reorganizar)
    ├── SLAM_DATA/
    │   └── <captura>.PosI             obligatorio
    ├── POST_SLAM_DATA/
    │   └── result_trajectory3d_high_freq.txt   obligatorio
    ├── GPS/Rtk/
    │   └── rtk.gga                    opcional (solo control de calidad)
    └── SCAN/Scanner1/
        ├── <fecha_hora>.pcap          obligatorio
        └── LAS-1.EP                   obligatorio con la calibración global y elevaciones "EP" (opción por defecto)
```

Las rutas de la ficha pueden ser absolutas o relativas a la carpeta de la ficha. La carpeta de brutos puede contener todo
lo demás que genera el RS30 (`GPS/Rover`, `INFO`, `TMP_*.bin`...): el pipeline solo busca los ficheros de la tabla.

## Qué se usa de cada fichero

| Fichero | Dónde lo busca (dentro de `entradas.brutos`) | Obligatorio | Formato esperado | Para qué |
|---|---|---|---|---|
| Nube | `entradas.las` (fuera de los brutos) | sí | LAS/LAZ con `gps_time` (formato de punto con tiempo GPS), CRS proyectado en metros en la cabecera o `crs.epsg` en la ficha | MH, MW, envolvente de copa y correspondencias con el pcap para calibrar |
| Trayectoria PosI | `SLAM_DATA/*.PosI` | sí | texto; 2 líneas de cabecera y columnas `t, lat, lon, h, roll, pitch, yaw` | georreferenciar la trayectoria (ajuste rígido local → CRS del ensayo) |
| Trayectoria post-SLAM | `POST_SLAM_DATA/result_trajectory3d_high_freq.txt` | sí | texto; columnas `t, x, y, z` (marco local) y 4 de orientación | posición y orientación del sensor a ~200 Hz |
| RTK | `GPS/Rtk/rtk.gga` | no | NMEA GGA | % de épocas con RTK fijo (aviso de deriva) |
| pcap del escáner | `SCAN/Scanner1/*.pcap` | sí | captura UDP del Hesai XT32 (registros de 862 B, protocolo v6.1) | todos los disparos, también los sin eco: Bt y Bi |
| Calibración del escáner | `SCAN/Scanner1/*.EP` | con `pcap.calibracion_global.elevaciones: "EP"` | `LAS-1.EP` (Protocol Buffers) | tabla de elevaciones de los 32 canales |

Notas:
- Hay un `result_trajectory3d_high_freq.txt` en `SLAM_DATA` y otro en `POST_SLAM_DATA`. Se usa el de **`POST_SLAM_DATA`**.
- El tiempo de todos los ficheros debe estar en la misma base (`gps_time` del LAS en segundos Unix UTC, igual que la
  trayectoria y la fecha UTC del pcap). El 01_qc comprueba que el `gps_time` del LAS cae dentro de la trayectoria.
- Si un patrón encuentra **más de un fichero** (p. ej. dos `.pcap`), el pipeline se para y pide que se indique cuál en
  `entradas.brutos_archivos` de la ficha. Ahí también se puede dar la ruta de cualquier fichero que esté en otro sitio.

## La nube LAS

- Debe ser la exportación del RS30 **de la misma captura** que los brutos, con las coordenadas y el `gps_time` originales:
  la calibración empareja cada punto del LAS con su eco en el pcap por el sello de tiempo. Una nube recortada, reproyectada
  con otra herramienta, fusionada con otras capturas o submuestreada rompe esas correspondencias.
- No hace falta que esté ordenada por tiempo ni clasificada: el pipeline clasifica el suelo y el ruido.
- Cuantos más ecos conserve, mejor calibra: en IFAPA Cabra, los tramos exportados a menor densidad de puntos por segundo
  fueron los que no pasaron la validación de la calibración.

## El archivo vectorial de hileras

- Una capa de **líneas** (GPKG, SHP o cualquier formato que lea GeoPandas/sf), **una línea por hilera**, trazada sobre su
  eje. Basta con una precisión de unos decímetros: el eje se reajusta con la copa y la longitud se toma de la nube.
- Un campo identificador (`croquis.campo_id`). Si se repite, `croquis.id_compuesto: true` lo combina con los atributos.
- Los atributos que se quieran arrastrar a los resultados (`croquis.atributos`: variedad, bloque, tratamiento...).
- `croquis.filtro` elige las hileras de esta captura (p. ej. `{bloque: 1, manejo: riego}`).
- Si está en otro CRS que la nube, se reproyecta en memoria.

## Comprobar las entradas antes de lanzar

```
python python/pipeline/comprobar_entradas.py --ensayo mi_ensayo.yaml
Rscript R/pipeline/comprobar_entradas.R --ensayo mi_ensayo.yaml
```

Tarda segundos y solo lee (cabeceras y unas pocas líneas o paquetes de cada fichero). Cada comprobación sale como **OK**,
**AVISO** (se puede lanzar, pero conviene mirarlo) o **ERROR** (el pipeline fallaría):

| Qué | Cómo |
|---|---|
| Ficha y parámetros | secciones obligatorias; `parametros.yaml` + `ajustes` cargan |
| Nube | se abre; formato de punto con `gps_time` y no todo a 0; CRS en la cabecera o en la ficha (aviso si no coinciden) |
| Ficheros de los brutos | cada uno se encuentra una sola vez (o está en `brutos_archivos`); el `LAS-1.EP` es obligatorio si la configuración lo usa |
| PosI y `high_freq` | número de columnas numéricas y rango de tiempo (primera y última línea) |
| pcap | cabecera de pcap; tramas del Hesai XT32 (862 B, payload EE FF) al principio y al final; rango de tiempo |
| `rtk.gga` | épocas GGA y % con RTK fijo |
| `LAS-1.EP` | se lee la tabla de 32 elevaciones y están en un rango razonable |
| Tiempos | el `gps_time` del LAS (primer millón de puntos) cae en la trayectoria, y la trayectoria en el PosI y en el pcap: si no se solapan, los ficheros no son de la misma captura |
| Archivo vectorial | capa de líneas, campo id único, atributos y campos del filtro presentes, el filtro deja hileras, hileras de prueba existentes, hileras dentro de la extensión de la nube |
| Salida | se puede crear, no está dentro de los brutos, espacio libre |

El informe queda en `<salida>/logs/comprobar_entradas.json` y el script termina con código 1 si hay algún ERROR.

Además, al arrancar el pipeline, `cargar_config` vuelve a comprobar que existen las entradas obligatorias y que hay CRS.
