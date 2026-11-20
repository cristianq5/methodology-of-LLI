# Methodology of LLI (Leafiness-LiDAR Index)

Pipeline para calcular el **Leafiness-LiDAR Index (LLI)** de cada hilera de un seto (olivar, viña, frutal) a partir de una captura con el LiDAR de mano **CHCNAV RS30** (escáner Hesai XT32 + SLAM). Sirve para cualquier ensayo capturado con ese equipo: se le da la nube, la carpeta de brutos y un archivo vectorial con las líneas de las hileras, y calcula.

```
LLI = MH × MW × (1 − P),   P = Bi / Bt
```
- **MH**: altura del seto.
- **MW**: anchura del seto.
- **P**: porosidad = haces que atraviesan la copa (Bi) / haces que entran en ella (Bt). Bi y Bt se cuentan por trazado de rayos con los disparos reales del escáner (el `.pcap`), porque el LAS no guarda los haces sin eco.

El cálculo se hace por rebanadas de 0,1 m  y se agrega por hilera y por tramos de 1,5 m.

El método está implementado dos veces, con la misma configuración: en **Python** y en **R**.

**Descripción completa del pipeline (entradas, pasos, salidas, calibración, limitaciones): [`docs/PIPELINE.md`](docs/PIPELINE.md).**
**Qué datos hacen falta y cómo organizarlos (nube, carpeta de brutos del RS30, hileras): [`docs/ENTRADAS.md`](docs/ENTRADAS.md).**

## Estructura

```
├── config/
│   ├── parametros.yaml        Parámetros del método (comunes a todos los ensayos)
│   └── ensayo_plantilla.yaml  Plantilla de la ficha de un ensayo
├── python/
│   ├── pipeline/              01_qc … 11_exportar, run_pipeline.py, módulos y tests/
│   ├── environment.yml        Entorno conda (recomendado)
│   └── requirements.txt
├── R/
│   ├── pipeline/              Réplica 01 … 11, run_pipeline.R y módulos
│   ├── instalar_paquetes.R
│   └── requirements_R.txt
└── docs/
    ├── PIPELINE.md            Descripción del pipeline
    └── ENTRADAS.md            Datos de entrada y organización de la carpeta
```

Las dependencias de cada lenguaje, y por qué se usan esas, están en [`python/README.md`](python/README.md) y [`R/README.md`](R/README.md).

## Uso rápido

1. **Instalar.**
   ```bash
   conda env create -f python/environment.yml && conda activate lidar-agrolab
   Rscript R/instalar_paquetes.R          # solo si se va a usar la versión de R
   ```
2. **Preparar la ficha del ensayo.** Copia `config/ensayo_plantilla.yaml` junto a los datos y rellena la nube, la carpeta de brutos del RS30, el archivo vectorial (campo id, atributos y filtro de hileras), las hileras de prueba y la carpeta de salida.
3. **Comprobar que están todos los datos** (segundos; solo lee). Corrige los ERROR antes de seguir:
   ```bash
   python python/pipeline/comprobar_entradas.py --ensayo /ruta/mi_ensayo.yaml
   ```
4. **Probar con un tramo corto** (modo `prueba`, el de la plantilla) y revisar los PNG y los logs:
   ```bash
   python python/pipeline/run_pipeline.py --ensayo /ruta/mi_ensayo.yaml
   ```
5. **Calcular todas las hileras:**
   ```bash
   python python/pipeline/run_pipeline.py --ensayo /ruta/mi_ensayo.yaml --modo completo
   Rscript R/pipeline/run_pipeline.R --ensayo /ruta/mi_ensayo.yaml --modo completo    # réplica en R
   ```

**Resultados.**
- **Final:** `<salida>/11_exportar/` (GPKG con hileras, tramos y rebanadas; CSV; GeoTIFF).
- **Resumen por hilera:** `11_exportar/consolidado/resumen_hileras.csv`.
- **Logs:** `<salida>/logs/`, uno por paso.

## Datos

El repositorio **no contiene datos**: el `.gitignore` impide subir nubes, pcap, trayectorias, vectoriales, salidas y logs. Los datos de entrada son de **solo lectura**: ningún paso escribe sobre ellos.

## Reglas del pipeline
- Todos los parámetros están en `config/` (con su origen anotado) o en la ficha del ensayo. No hay valores fijos en el código.
- Cada resultado es reproducible: la configuración efectiva se guarda junto a él.
- Validez de la calibración del escáner: **mediana < 0,15 m y p90 < 0,30 m**, con datos no usados en el ajuste. Si una calle no cumple, sus rebanadas quedan sin P oficial y el motivo se escribe en el resultado.
- Documentación y comentarios en español.
