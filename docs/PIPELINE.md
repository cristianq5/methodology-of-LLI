# Descripción del pipeline LLI

El pipeline calcula el **Leafiness-LiDAR Index** de cada hilera de un seto (olivar, viña, frutal en seto) a partir de una captura del LiDAR de mano **CHCNAV RS30** (escáner Hesai XT32 + SLAM):

```
LLI = MH × MW × (1 − P)        P = Bi / Bt
```

| Símbolo | Qué es | Cómo se obtiene |
|---|---|---|
| MH | altura del seto (m) | máximo (estricta) o percentil 99 (robusta) de la altura sobre el suelo en cada rebanada |
| MW | anchura del seto (m) | media, por franjas horizontales de 0,1 m, de la anchura `t_max − t_min` (estricta) o entre percentiles 2,5–97,5 (robusta) |
| Bt | haces que entran en la copa | disparos reales del escáner (pcap) que entran en la envolvente de copa por su lado |
| Bi | haces que atraviesan la copa | de esos, los que salen por el otro lado sin que su eco caiga dentro |
| P | porosidad (0–1) | `Bi / Bt` por celda (rebanada × franja), por lado y media de los dos lados |

El cálculo se hace por **rebanadas de 0,1 m** a lo largo de la hilera. Después se agrega por hilera y por tramos de 1,5 m, tras quitar las rebanadas a más de ±2,5 desviaciones típicas.

**Por qué el pcap.** El LAS solo guarda los ecos. Los haces que atraviesan la copa sin devolver nada (los que hacen falta para Bi) solo están en el `.pcap` del escáner, como disparos de distancia 0. Por eso el pipeline decodifica el pcap y reconstruye cada rayo en el mundo con la trayectoria y la calibración escáner–trayectoria (extrínseca) de cada calle.

---

## Entradas de un ensayo

Todo se describe en una **ficha de ensayo** (`config/ensayo_plantilla.yaml`). La organización de la carpeta y el formato de cada fichero están en [`ENTRADAS.md`](ENTRADAS.md).

1. **Nube LAS/LAZ** exportada por el RS30, con `gps_time`.
2. **Carpeta de brutos** de la captura, tal como la deja el RS30. El pipeline busca dentro:

   | Fichero | Ruta dentro de los brutos | Uso |
   |---|---|---|
   | trayectoria PosI | `SLAM_DATA/*.PosI` | posición absoluta (lat/lon/h) para llevar la trayectoria al CRS del ensayo |
   | trayectoria post-SLAM | `POST_SLAM_DATA/result_trajectory3d_high_freq.txt` | posición y orientación a 200 Hz |
   | RTK | `GPS/Rtk/rtk.gga` | solo control de calidad (opcional) |
   | pcap del escáner | `SCAN/Scanner1/*.pcap` | disparos con y sin eco (Bt, Bi) |

   Si un fichero no está o hay más de un candidato, el pipeline se para y pide que se indique en `entradas.brutos_archivos`.
3. **Archivo vectorial de hileras**: una capa de **líneas**, una por hilera, trazada sobre el eje. Se indica qué campo es el identificador (`campo_id`), qué atributos se copian a los resultados (variedad, bloque, tratamiento...) y qué hileras se procesan (`filtro`). El eje se reajusta después con los puntos de copa, así que basta con un archivo vectorial aproximado (unos decímetros). La **longitud** de la hilera tampoco sale de la línea: con `hilera.extension: nube` (por defecto) el inicio y el final se toman de la copa en la nube, y la línea solo identifica la hilera y sirve de referencia (ver 08).

Los datos de entrada se abren en **solo lectura**: el pipeline nunca los modifica. Las correcciones (minúsculas, erratas, reproyección del archivo vectorial) se hacen en memoria.

---

## Pasos

Cada paso escribe en `<salida>/NN_nombre/` (en modo completo, en `<salida>/NN_nombre/hilera_<id>/`) y deja su log en `<salida>/logs/NN_nombre[_hilera_<id>].log`. El log recoge la entrada, los parámetros, los puntos antes y después, la salida y el tiempo.

| Paso | Qué hace | Salidas principales |
|---|---|---|
| **01_qc** | Control de calidad sin tocar nada. Lee la cabecera del LAS y lleva la trayectoria post-SLAM al CRS del ensayo con un ajuste rígido frente al PosI. Calcula el % de RTK fijo, comprueba la cobertura bilateral de cada hilera (que se haya recorrido por las dos calles) y busca los pasos del sensor | `qc_resumen.json`, `qc_cobertura_bilateral.csv`, `trayectoria_utm.npz`, `mapa_trayectoria_filas.png` |
| **02_recorte** | Extrae del LAS, por trozos, la zona de trabajo: en modo *prueba*, un tramo de las hileras de prueba; en modo *completo*, un bloque por hilera: la línea del archivo vectorial ampliada `hilera.margen_busqueda_m` (1 m: entre bloques hay poca distancia) por cada extremo, como zona de búsqueda. A lo ancho toma ±`recorte.buffer_recorte_m` (2,4 m) del eje, más que el buffer de la hilera (1,2 m), para conservar el suelo de la calle y el «fondo» de los rayos. **No submuestrea** | `recorte.las`, `recorte_meta.json` |
| **03_filtrado_sor** | Detecta puntos aislados (SOR de PDAL, k = 10, ×2,5). Se **marcan** como clase 7, no se borran | `recorte_sor.las` |
| **04_suelo** | Clasifica el suelo con CSF (o SMRF) sin contar el ruido | `recorte_suelo.las` |
| **05_dtm** | Crea el DTM (IDW, 0,25 m) solo con el suelo de las calles, porque bajo el seto suele haber hierba, y rellena los huecos. Reclasifica como no-suelo los puntos que quedan muy por encima. Comprueba la altura del sensor sobre el DTM | `dtm.tif`, `qa_suelo.json`, `qa_altura_sensor_sobre_dtm.json` |
| **06_normalizacion** | Calcula la altura sobre el suelo (HAG) con el DTM. Conserva todos los puntos: el umbral de 0,3 m se aplica después | `recorte_hag.las` |
| **07_chm** | Modelo de altura de copa (máximo de HAG en celdas de 0,1 m), solo como control visual | `chm.tif`, `chm.png` |
| **08_segmentacion** | Asigna los puntos de cada hilera (buffer de 1,2 m a cada lado del eje, 2,4 m de ancho), reajusta el eje por PCA y calcula las coordenadas de hilera: `s` a lo largo, `t` perpendicular (+ a la izquierda) y lado. Fija la **longitud real de la hilera con la nube**: un tramo de 0,1 m es copa si tiene puntos de HAG ≥ 0,5 m, al menos 20 y al menos el 15 % de un tramo típico de la hilera. Dentro de la línea del archivo vectorial, la hilera va del primer al último tramo con copa (un olivo que falte no la corta). Si la copa llega al extremo de la línea, se hace un **barrido fuera en pasos de 0,1 m hasta 1 m** (`extension_paso_m`, `margen_busqueda_m`) y el límite es el último paso con copa pegada a la hilera; así no entra la copa del bloque vecino ni la vegetación suelta de la cabecera. `hilera.recorte_extremos_m` (0 por defecto) quitaría metros en cada cabecera después. `marcos.json` guarda la longitud de la nube y la de la línea | `hileras.npz`, `marcos.json`, `segmentacion.las`, `segmentacion_planta.png` |
| **09_metricas** | MH, MW y envolvente de copa por rebanada y franja. Decodifica el pcap en las ventanas de cada calle, calibra la extrínseca de cada calle y la **valida con datos no usados en el ajuste**. Después clasifica los rayos (Bt, Bi) y calcula P por rebanada y por lado. Con la misma malla rebanada × franja calcula las **métricas de copa** (ver más abajo) y el **perfil de anchura por altura** | `metricas_rebanadas.csv`, `perfil_anchura.csv`, `extrinseca.json`, `_calles/calle_*.json`, `comparacion_P.json`, `perfil_P_hilera*.png`, `perfil_anchura_hilera*.png` |
| **10_lli** | LLI por rebanada (variantes estricta y percentil), filtro de ±2,5 DE y medias por hilera y por tramo de 1,5 m. Agrega también las métricas de copa (medias y volumen total). Comprueba el signo de L = 1 − P | `lli_rebanadas.csv`, `lli_hileras.csv`, `lli_tramos.csv`, `lli_resumen.json` |
| **11_exportar** | Exporta un GPKG (capas hileras, tramos y rebanadas), CSV (también el perfil de anchura), GeoTIFF (DTM, CHM, LLI), la nube con dimensiones extra (HAG, id_hilera, s, t, lado) y la configuración usada | `lli_resultados.gpkg`, CSV, TIF, `nube_hileras.las`, `ficha_ensayo_usada.yaml`, `config_efectiva.yaml` |

En modo completo, `run_pipeline` reparte las hileras entre `ejecucion.n_procesos` procesos y junta los CSV en `11_exportar/consolidado/`. También escribe `resumen_hileras.csv`, que dice para cada hilera si está completa o no y por qué, y anota la RAM pico en `logs/ram_pico.json`.

---

## Métricas de copa (09 y 10)

Además de MH, MW y P, el paso 09 calcula otras métricas habituales en LiDAR de cultivos leñosos. Salen de la misma malla que MH y MW: rebanadas de 0,1 m a lo largo de la hilera y franjas horizontales de 0,1 m desde 0,3 m. Por eso no hace falta leer otra vez la nube. Solo cuenta la copa por encima de 0,3 m.

| Columna (por rebanada) | Qué es |
|---|---|
| `area_seccion_m2` (`_p_m2`) | Área de la sección transversal del seto: suma de la anchura de cada franja ocupada × 0,1 m. Con la anchura estricta (t_max − t_min) o la robusta (percentiles 2,5–97,5). A diferencia de MH × MW, no supone que la sección sea un rectángulo |
| `vol_copa_m3` (`_p_m3`) | Volumen de copa de la rebanada = área × 0,1 m |
| `h_base_copa_m` | Altura de inicio de copa (falda): borde inferior de la franja ocupada más baja. Si la copa baja de 0,3 m, vale 0,3 (no se mide por debajo) |
| `prof_copa_m` | Profundidad de copa = MH − altura de inicio de copa |
| `anchura_max_m`, `h_anchura_max_m` | Anchura de la franja más ancha y su altura (centro de la franja) |
| `TRV_m3_ha` (`TRV_p_m3_ha`) | Tree Row Volume = MH × MW × 10 000 / distancia entre hileras (m³/ha) |
| `dist_hileras_m`, `origen_dist_hileras` | Distancia entre hileras usada en el TRV: `metricas.distancia_entre_hileras_m` si se indica; si no, la distancia del punto medio de la línea de la hilera a la línea más próxima de otra hilera del archivo vectorial (sin el filtro del ensayo) |

En `lli_hileras.csv` y `lli_tramos.csv` se agregan como `<columna>_media`. Los volúmenes se agregan como `vol_copa_m3_total` (m³ de la hilera o del tramo) y `vol_copa_m3_por_m` (m³ por metro de hilera). Estas métricas **no pasan el filtro de ±2,5 DE** del LLI: las rebanadas sin copa ya son NaN.

**Perfil de anchura por altura** (`perfil_anchura.csv`, una fila por franja). Se calcula para la hilera entera (`nivel = hilera`) y para cada tramo de 1,5 m (`nivel = tramo`), solo con las rebanadas válidas:
- `n_ocupadas`, `frac_ocupada`: en cuántas rebanadas hay copa a esa altura;
- `anchura_media_m` (`_p_m`): anchura media donde hay copa;
- `anchura_media_con_huecos_m`: anchura media contando 0 donde no hay copa (= anchura media × fracción ocupada).

`perfil_anchura_hilera<id>.png` lo dibuja para cada hilera.

---

## Calibración escáner–trayectoria (09)

El XT32 del RS30 gira lentamente sobre su eje (unas 30 rpm), así que la relación entre escáner y trayectoria cambia con el tiempo. El fabricante no documenta el fichero que la describe (`.EP`). Por eso `extrinseca.py` la estima para cada **calle**, es decir, cada paso continuo del sensor por una calle en un sentido:

- **Modelo.** Giro del cabezal θ(t), base de giro → cuerpo y corrección rígida del marco local al mundo, con deriva por tramos.
- **Ajuste.** Se hace con correspondencias exactas entre ecos del pcap y puntos del LAS que comparten el sello de tiempo.
- **Validación (hold-out).** Se usan otros grupos de tiempo: la calibración guarda los sellos de todos los grupos que usó (`k_ajuste`) y la validación los excluye del sorteo. Una calibración sin `k_ajuste` (versión anterior) no se da por válida y hay que recalibrarla. Una calle es válida si el residuo cumple **mediana < 0,15 m y p90 < 0,30 m** (`pcap.extrinseca.residuo_*_max_m`). Este umbral no se cambia.
- **Calles no válidas.** Sus rayos no entran en la P oficial. Se registra el motivo y se ofrece, aparte, `P_mejor` / `via_P_mejor`, con prioridad: pcap válido, luego pcap no válido, luego rayos desde la trayectoria hasta los puntos LAS.

Las calibraciones se guardan en `09_metricas/_calles/` con una huella (parámetros, ventana, pcap, código y semilla). Se reutilizan solo si la huella coincide.

### Calibración global por bloque (`pcap.modo_calibracion: "global_bloque"`, por defecto)

Con la lectura original de la orientación de la trayectoria (columnas 5–8 del `high_freq` como qx qy qz qw), cada calle necesitaba su propio montaje y muchas no pasaban la validación. Con otra lectura de esas columnas (`pcap.calibracion_global.lectura_cuaternion`, `ultimo_yzx_-+-`) y la tabla de elevaciones del fichero `LAS-1.EP` del RS30, **un único montaje del escáner vale para todas las calles del bloque**. En IFAPA Cabra (experimento del 2026-11-17/18) quedaron válidas 39 de 45 calles, frente a 12 de 45 calle a calle.

- **Común a todas las calles:** giro escáner → cuerpo, eje del motor, palanca, escala, desfases de azimut y elevación, y retardo de reloj.
- **Propio de cada calle:** estado del motor (ω, φ y nudos de fase) y la corrección lenta entre la trayectoria `high_freq` y la que usó el fabricante para el LAS (dW, dV, dR y nudos de dW).
- **Flujo (`run_pipeline`):**
  1. Pasos 03–08 y el paso 09 en su fase `calibrar` en todas las hileras: calibración por calle como antes, más los grupos de ajuste y de validación (`_calles/calle_*_grupos.pkl` o `.rds`).
  2. **09b** (`calibracion_global.py` / `.R`): ajuste conjunto del montaje. Arranca de la calle del bloque con menor p90 de validación (entre las válidas, si las hay).
  3. Paso 09 completo (P con el montaje global), 10 y 11.
- **Validación:** el mismo hold-out por calle y el mismo umbral. Las calles no válidas con el montaje global se excluyen de P.
- **Salidas:** `09_metricas/_global/montaje_bloque.json`, con huella, la validez de cada calle con el montaje global y su validez calle a calle para comparar.
- **R:** mismo modelo. El ajuste conjunto se hace con `minpack.lm::nls.lm` y un jacobiano por diferencias finitas agrupadas, equivalente al jacobiano disperso de `scipy`.

---

## Configuración

| Fichero | Contenido |
|---|---|
| `config/parametros.yaml` | Parámetros del método, comunes a todos los ensayos. Cada valor indica su origen: `[FICHA]` método original, `[DATO]` comprobado con datos del RS30, `[SUP]` supuesto que conviene revisar con un análisis de sensibilidad |
| ficha del ensayo | Entradas, archivo vectorial, CRS, modo (prueba / completo), hileras de prueba, salida y `ajustes:` (cualquier valor de `parametros.yaml` que cambie en ese ensayo) |

La configuración efectiva (parámetros + ficha + ajustes) se guarda con los resultados (`config_efectiva.yaml`), de modo que cada resultado queda reproducible.

**CRS.** Se trabaja en el CRS proyectado (en metros) de `crs.epsg` o, si la ficha no lo indica, en el de la cabecera del LAS. Si el archivo vectorial está en otro CRS, se reproyecta en memoria.

---

## Cómo se lanza un ensayo nuevo

1. Copia `config/ensayo_plantilla.yaml` junto a los datos del ensayo y rellena `entradas`, `croquis`, `ejecucion_ensayo` y `salidas` (qué debe haber en cada carpeta: [`ENTRADAS.md`](ENTRADAS.md)).
2. **Comprueba las entradas** (unos segundos, solo lee): `python python/pipeline/comprobar_entradas.py --ensayo mi_ensayo.yaml` (o `Rscript R/pipeline/comprobar_entradas.R --ensayo mi_ensayo.yaml`). Revisa que la ficha y los parámetros cargan, la nube (gps_time, CRS), cada fichero de los brutos, que los tiempos del LAS, la trayectoria y el pcap se solapan (misma captura), el archivo vectorial (campos, filtro, hileras de prueba, que las hileras caigan dentro de la nube) y la carpeta de salida. Deja el informe en `<salida>/logs/comprobar_entradas.json` y termina con código 1 si hay algún ERROR.
3. Lanza el **modo prueba** (un tramo de 10–15 m de 1–3 hileras vecinas) y revisa los PNG de `01_qc` y `08_segmentacion`, la validez de las calles en el log de `09_metricas` y `lli_hileras.csv`.
4. Si hace falta, ajusta parámetros en `ajustes:` (por ejemplo `hilera.buffer_m` si las calles son más estrechas, `qc.distancia_tray_eje_esperada_m` o `metricas.distancia_entre_hileras_m` para el TRV).
5. Lanza el **modo completo** (`--modo completo`).

---

## Limitaciones

- **El pipeline está hecho para el RS30.** Depende de la estructura de carpetas de sus brutos, del formato de su trayectoria y del pcap del Hesai XT32 (paquetes de 32 canales, unidad de 5 mm). Para otro sensor habría que adaptar `trayectoria.py`, `pcap_hesai.py` y `extrinseca.py`.
- **Ángulos verticales por canal.** Con la calibración global se usa la tabla del fichero `LAS-1.EP` (campo 7 de la tabla 7 del registro 0). Su significado no lo documenta el fabricante, pero mejora el residuo frente al LAS. Con `elevaciones: "xt32_estandar"` (o calle a calle) se usan los estándar del XT32 (de +15° a −16°), no verificados con el manual.
- **Calles que no calibran.** Si una calle no cumple el umbral de validación, sus rebanadas quedan sin P oficial. El motivo se escribe en `P_motivo_nan`.
- **Supuestos por revisar.** Los parámetros marcados `[SUP]` son supuestos iniciales y conviene revisarlos con un análisis de sensibilidad en cada tipo de plantación.
