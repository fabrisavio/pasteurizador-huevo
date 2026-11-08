# Pasteurizador de ovoproductos: modelo paramétrico de diseño

Una planta de huevo en polvo seca por spray tres productos líquidos: **huevo entero**, **yema** y **clara concentrada** por ultrafiltración. Antes del secador, los tres pasan por el mismo pasteurizador de placas, en campañas distintas. Cada uno empuja el diseño hacia un lado diferente:

| Producto | Caudal (kg/h) | Entra a (°C) | Tratamiento | Lo que complica |
|---|---|---|---|---|
| Huevo entero (HP) | 798 | 10–25 | 60,0 °C × 3,5 min | Más horas de operación al año: domina la energía |
| Yema (YP) | 1083 | 7 | 61,1 °C × 3,5 min | Decenas de veces más viscosa y pseudoplástica: domina el tubo de retención |
| Clara concentrada (CP) | 822 | 35 | 56,7 °C × 4,96 min, a pH 9,3 | Es la más sensible al calor. Tras el desazucarado (pH 7) se lleva a pH 9,3 con soda: así alcanza con 56,7 °C, pero empieza a espesar justo ahí y el agua no puede pasar de 59 °C |

**La pregunta:** ¿qué equipo único cumple con los tres, y cómo funciona cada uno adentro de ese equipo?

El modelo responde en dos pasadas. Primero calcula qué necesitaría cada producto por separado. Después arma un equipo donde, en cada sección, manda el producto que más pide. Al final verifica cada producto en ese equipo real, en todo su rango de temperatura de entrada.

```
crudo ─► [regeneración, lado frío] ─► [calentamiento] ─► tubo de retención ─┐
                                          ▲ agua caliente                   │
al secador ◄── [regeneración, lado caliente] ◄─────────────────────────────┘
```

Todos los datos de entrada viven en `config/base.yaml`, cada uno con su fuente. Si cambia un caudal o una temperatura, se cambia ahí y el diseño se recalcula solo.

## Resultados de la configuración base (v0.5.0)

| Elemento | Resultado | Lo gobierna |
|---|---|---|
| Regeneración | 10 placas de 0,15 m² (1,50 m²) | Clara: necesita ε ≈ 0,55 contra ~0,3 de huevo y yema |
| Calentamiento | 24 placas (3,60 m²) | Yema |
| Tubo de retención | DN40, 97,5 m enrollados en una espiral de 1 m de diámetro × 1,6 m de alto (31 vueltas), 111 L, Δp 0,38 bar | Clara: 4,96 min (la yema queda casi igual, con 3,5 min y más caudal) |
| Aislación del tubo | Camisa desmontable de 50 mm alrededor de toda la espiral: pierde 0,3 °C y 0,74 MWh/año | Contra coquilla a lo largo del tubo: 0,9 °C, 2,2 MWh/año y 27 placas de calentamiento en vez de 24 |
| Calor | 49,4 MWh/año | Huevo entero: 88 % del total |

Como la regeneración se dimensiona para la clara, huevo y yema regeneran de más y llegan al secador a unos 31–32 °C en lugar de 45 °C. Es aceptable y queda informado para el diseño del atomizador.

**El tubo no se achica cambiando el diámetro.** En laminar, el volumen lo fijan el caudal, los 3,5 min y el factor de la partícula más rápida; el DN solo lo reparte entre largo y Δp. Enrollarlo en espiral cambia el factor: en cada vuelta la fuerza centrífuga arma vórtices de Dean que sacan a la partícula rápida del centro del tubo. En una hélice, la partícula más rápida tarda al menos 0,586 veces el tiempo medio (Ruthven 1971; Trivedi y Vasudeva 1974; Gutierrez et al. 2010 midieron 0,70 o más), contra 0,5 en tubo recto. Con ese valor, el más conservador publicado, el tubo es un 15 % más corto que en recto. A cambio sube la fricción por vuelta: DN32 queda al borde del Δp máximo y se descarta con el factor de seguridad sobre la viscosidad de la yema.

**Aislar la espiral entera sale mejor que aislar el tubo.** Una camisa alrededor de la espiral tiene 9 m² de superficie; la coquilla a lo largo de 97 m de tubo, 43 m². Con la camisa el producto pierde un tercio, y eso no solo ahorra aislante y calor: como el agua caliente trabaja a 1–2 °C del tope de cada producto, cada décima que el producto no pierde en el tubo es área que no hace falta en el calentamiento. Son 3 placas menos.

**Validación.** El mismo código reproduce el pasteurizador de yema publicado por Gut et al. (2003), Tabla 7, con diferencias de 0,02 °C o menos en las cinco temperaturas que reporta el paper: la yema a la salida de la regeneración (T2), del calentamiento (T3), del tubo de retención (T4), de la regeneración del lado pasteurizado (T5) y del enfriamiento (T6). Tolerancia: 0,5 °C. `python correr.py --validar` muestra la tabla con un esquema de dónde está cada punto.

## Cómo correrlo

Requiere Python 3.10 o superior.

```bash
python -m venv .venv
source .venv/Scripts/activate     # Windows (Git Bash). En Mac/Linux: source .venv/bin/activate
pip install -e ".[dev]"

pytest                            # 60 tests, incluida la validación
python correr.py                  # diseña con config/base.yaml → resultados/reporte.md
python correr.py --validar        # compara con Gut et al. (2003)
```

Para probar otro caso (otro caudal, otro diámetro de espiral, tubo recto…), copiá la configuración base, editá la copia y corré esa:

```bash
cp config/base.yaml config/mi_caso.yaml
python correr.py config/mi_caso.yaml
```

El reporte separa **alertas** (algo no cumple: hay que actuar) de **avisos** (supuestos, extrapolaciones y datos que necesitan otros equipos).

## Estructura

```
config/
  base.yaml                 inputs del diseño, con fuentes
  validacion_gut2003.yaml   caso publicado para validar
src/pasteurizador/
  propiedades.py            Choi–Okos (térmicas), reología medida por producto, agua IAPWS-IF97
  cinetica.py               modelo D–z: logaritmos logrados
  holding.py                tubo de retención: Metzner–Reed, partícula más rápida, espiral, elección de DN
  aislacion.py              pérdida de calor del tubo: coquilla o camisa en bloque
  seccion_phe.py            una sección de placas en contracorriente (ε-NTU), diseño y verificación
  pasteurizador.py          las dos pasadas: diseño por producto, equipo único, verificación
  validacion.py             encadena las secciones del caso de Gut
  reporte.py                reporte en Markdown y tablas CSV
tests/                      balances, límites físicos, validación
docs/
  pseudocodigo.md           el diseño del modelo, paso a paso
  justificacion_clara.tex   por qué la clara va a pH 9,3 con soda, listo para el informe
  fuentes.bib               todas las fuentes, listas para LaTeX
correr.py                   punto de entrada
```

## Hoja de ruta

| Versión | Qué agrega | Estado |
|---|---|---|
| v0 | U fijo por sección, tres productos, diseño y verificación, tubo con Metzner–Reed, validación | Hecho (0.1.0) |
| v0.2 | Tubo enrollado en espiral: fricción por curvatura, Re crítico corrido | Hecho (0.2.0) |
| v0.3 | Partícula más rápida en la hélice con valor de literatura, recorrido de T del producto en consola y reporte, validación explicada | Hecho (0.3.0) |
| v0.4 | Pérdida del tubo calculada con la aislación (coquilla o camisa en bloque) y comparada, margen de control sobre la T legal, factor de seguridad en el Δp del tubo | Hecho (0.4.0) |
| v0.5 | Clara a pH 9,3 con soda y 56,7 °C × 4,96 min (FSIS, Tabla III.A); aviso de espesamiento por producto; pH en el reporte. Cierre de v0 | Hecho (0.5.0) |
| v1 | Coeficientes h y fricción con correlaciones de placa chevron (Kumar 1984; Martin 1996 de contraste), Δp en placas, chequeo de CIP | Próximo |
| v2 | Perfil axial en las placas y letalidad acumulada (F) de placas más tubo | — |
| Después | Sensibilidad (tornado → Sobol), optimización de costo anual, Monte Carlo | — |

## Supuestos y límites de v0

- **U fijo** de 756 W/m²K en regeneración y 864 W/m²K en calentamiento, tomado del caso de yema de Gut et al. (2003) con placas lisas. Se usa para los tres productos hasta que v1 calcule h con correlaciones.
- **Sin Δp en placas**: necesita el factor de fricción de la correlación (v1). Solo se calcula el Δp del tubo.
- **Tubo en transición** (2100 < Re < 4000 en recto; en espiral, la banda se corre hacia arriba): se diseña con el factor laminar de la partícula más rápida, que es el lado seguro.
- **Espiral:** las correcciones de fricción de White (1929) e Ito (1959) y el t_min/t_medio = 0,586 son resultados para newtonianos; se aplican también a yema y huevo entero, que son levemente pseudoplásticos y tienen un perfil más plano. El diámetro de la hélice (1 m) es un supuesto.
- **Aislación:** lana mineral de 50 mm, k = 0,040 W/m·K, sala a 20 °C. En la camisa en bloque se supone el aire interior a la T del producto, que es la cota superior de la pérdida.
- **Margen de control:** el control trabaja 0,5 °C sobre la T legal (supuesto: el doble de la histéresis mínima de la válvula de desvío).
- **Δp del tubo con factor 1,4:** cubre la mayor diferencia de K de la yema entre la ecuación y la tabla de Gut et al. (2003). Sin él, DN32 entra con 1 % de margen.
- **Letalidad solo en el tubo.** Las placas también matan; Gut encontró que suben el F de 3,5 a 5,5 min. Eso entra en v2.
- **Reologías extrapoladas** a la temperatura de retención en los tres productos. El reporte lo avisa.
- **Espesamiento:** huevo entero (desde 56 °C) y clara a pH 9 (desde 56,7 °C) pasan por el tubo arriba de la T donde empiezan a espesar. Las reologías medidas no llegan ahí: el Δp real del tubo puede ser algo mayor. El reporte lo avisa.
- **Clara a pH 9,3:** la forma de la viscosidad con T es la de Scalzo para clara a pH 7; se supone igual a pH 9,3. El sodio de la soda (~0,4 g/kg) no se agrega a la composición de Choi–Okos.
- **Desvíos del RBN declarados:** la clara se pasteuriza a 56,7 °C × 4,96 min y no a 60 °C × 3,5 min (15.2.11), y lleva soda, que no está en la lista de aditivos para huevo. La justificación está en `docs/justificacion_clara.tex`.

## Fuentes

Las claves entre corchetes en `config/base.yaml` remiten a `docs/fuentes.bib`. Las principales: FSIS (2026) y el Reglamento Bromatológico Nacional para los pares T–t; Choi y Okos (1986) para propiedades térmicas; Gut et al. (2003), León Munguía, Thapon y Scalzo et al. (1970) para reología; Metzner y Reed (1955) para el tubo; Ruthven (1971), Trivedi y Vasudeva (1974) y Gutierrez et al. (2010) para la partícula más rápida en la hélice; FSIS (2026, Tabla III.A) y el manual AEB para la clara a pH 9,3; Health Canada (2013), 21 CFR 184.1763 y el Codex GSFA para la soda como regulador de pH; CFIA para la operación del tubo de retención (bomba de tiempo, válvula de desvío, aislación desmontable).

---

Proyecto Industrial de Ingeniería Química, Facultad de Ingeniería, Udelar (2026): *Producción de huevo en polvo y valorización de la cáscara*.
