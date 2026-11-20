# Configuración

Python y R leen **la misma** configuración, repartida en dos ficheros. Los scripts no llevan parámetros fijos.

| Fichero | Qué contiene | Quién lo cambia |
|---|---|---|
| `parametros.yaml` | Parámetros del **método**, comunes a todos los ensayos: filtrado, suelo, DTM, buffer de hilera, rebanadas, porosidad, calibración, LLI | casi nunca; un ensayo concreto lo ajusta en su ficha |
| `ensayo_plantilla.yaml` | Plantilla de la **ficha de ensayo**: LAS, carpeta de brutos del RS30, archivo vectorial de hileras (campo id, atributos, filtro), CRS, modo prueba/completo, salida y `ajustes:` | se copia y se rellena una por ensayo |

La configuración efectiva es `parametros.yaml` + la ficha + sus `ajustes:`. Se guarda con los resultados (`11_exportar/config_efectiva.yaml` y `ficha_ensayo_usada.yaml`).

Cada valor de `parametros.yaml` indica su origen:
- `[FICHA]`: método LLI original.
- `[DATO]`: comprobado con datos del RS30.
- `[SUP]`: supuesto que hay que justificar con un análisis de sensibilidad.

## Rutas
En la ficha, las rutas relativas se resuelven desde la carpeta donde está la ficha. Si se deja una ficha junto a los datos del ensayo, la misma ficha sirve en cualquier máquina.

## Reglas fijas
- **Umbral de validez de la calibración del escáner:** mediana < 0,15 m y p90 < 0,30 m, evaluado con datos no usados en el ajuste (`pcap.extrinseca.residuo_*_max_m`). No se cambia.
- **Datos originales en solo lectura:** las salidas se escriben siempre en la carpeta `salidas` de la ficha.
