# Pasteurizador de ovoproductos: modelo paramétrico de diseño

Una planta de huevo en polvo seca por spray tres productos líquidos: **huevo entero**, **yema** y **clara concentrada** por ultrafiltración. Antes del secador, los tres pasan por el mismo pasteurizador de placas, en campañas distintas. Cada uno empuja el diseño hacia un lado diferente:

| Producto | Caudal (kg/h) | Entra a (°C) | Tratamiento | Lo que complica |
|---|---|---|---|---|
| Huevo entero (HP) | 798 | 10–25 | 60,0 °C × 3,5 min | Más horas de operación al año: domina la energía |
| Yema (YP) | 1083 | 7 | 61,1 °C × 3,5 min | Decenas de veces más viscosa y pseudoplástica: transfiere peor y es la que más Δp hace en las placas |
| Clara concentrada (CP) | 822 | 35 | 56,7 °C × 4,96 min, a pH 9,3 | Es la más sensible al calor. Tras el desazucarado (pH 7) se lleva a pH 9,3 con soda: así alcanza con 56,7 °C, pero empieza a espesar justo ahí y el agua no puede pasar de 59 °C |

**La pregunta:** ¿qué equipo único cumple con los tres, y cómo funciona cada uno adentro de ese equipo?

El modelo responde en dos pasadas. Primero calcula qué necesitaría cada producto por separado. Después arma un equipo donde, en cada sección, manda el producto que más pide. Al final verifica cada producto en ese equipo real, en todo su rango de temperatura de entrada.

```
crudo ─► [regeneración, lado frío] ─► [calentamiento] ─► tubo de retención ─┐
                                          ▲ agua caliente                   │
al secador ◄── [regeneración, lado caliente] ◄─────────────────────────────┘
```

Todos los datos de entrada viven en `config/base.yaml`, cada uno con su fuente. Si cambia un caudal o una temperatura, se cambia ahí y el diseño se recalcula solo.

## Resultados de la configuración base (v1.1.0)

| Elemento | Resultado | Lo gobierna |
|---|---|---|
| Regeneración | 15 placas de 0,15 m² (2,25 m²), 1 paso × 8 canales | Clara: necesita ε ≈ 0,55 contra ~0,3 de huevo y yema. Más pasos darían menos placas, pero la yema pasaría 1 bar |
| Calentamiento | 47 placas (7,05 m²), 4 pasos × 6 canales | Yema: el U más bajo (483 W/m²K) y la que más calor pide. El agua caliente pierde 0,96 bar de 1 admitido; con 3 pasos serían 53 placas y 0,36 bar |
| U en operación | Regeneración 339–504 W/m²K; calentamiento 483–660 W/m²K | Kumar (1984) para placa chevron de 45°; el ensuciamiento de diseño pone el 20–57 % de la resistencia |
| Tubo de retención | DN50 (fijado), 66 m en horquillas: 11 tramos de 5,8 m y 10 curvas de 180°, en dos planos. Sube 2 % continuo: 1,3 m de desnivel. Envolvente 6,1 × 0,33 × 1,5 m. 130 L, Δp 0,09 bar | Clara: 4,96 min (la yema queda casi igual, con 3,5 min y más caudal) |
| Aislación del tubo | Coquilla desmontable de 50 mm en tramos y curvas: pierde hasta 0,75 °C y 1,8 MWh/año | Con las curvas al aire: 1,0 °C, 2,3 MWh/año y 53 placas de calentamiento en vez de 47 |
| CIP | Placas: 2,6 m³/h para Re ≥ 2150; tubo: 10,6 m³/h a 1,5 m/s | Manda el tubo |
| Calor | 51,4 MWh/año | Huevo entero: 86 % del total |

**Qué cambió de v0.5 a v1.** En v0.5 el U era un dato: 756 y 864 W/m²K, del pasteurizador de yema de Gut et al. (2003). Con la correlación y los datos que publica el propio paper, ese U no se reconstruye: sale 2,6 a 3,8 veces menor (`python correr.py --validar`). En v1 cada U sale de los dos coeficientes de película, la pared y el ensuciamiento, para cada producto y cada sección. Sin ensuciamiento, v1 pide las mismas 34 placas que v0.5: el U de v0.5 equivalía a un equipo limpio. Con el ensuciamiento de diseño que usa Gut (el de leche de Lalande et al., 1979) pide 62. Ese margen es la decisión más grande que queda abierta. Para volver al cálculo de v0.5 alcanza con `modelo_U: fijo` en la config.

Como la regeneración se dimensiona para la clara, huevo y yema regeneran de más y llegan al secador a unos 31 °C (huevo entero) y 38 °C (yema) en lugar de 45 °C. Es aceptable y queda informado para el diseño del atomizador.

**El tubo en horquillas (v1.1).** Hasta v1.0 el tubo iba enrollado en espiral; a pedido de los tutores pasa a horquillas: tramos rectos unidos por curvas de 180°, que es la forma habitual. Tres cosas cambian. Primero, se pierde la ayuda de la hélice: en tubo recto la partícula más rápida va al doble de la media (contra 1,71 veces en la hélice), así que en el mismo DN40 el tubo se alarga de 97,5 a 114 m. Las curvas mezclan algo, pero no se cuenta: es el lado seguro. Segundo, CFIA pide que el tubo suba al menos 2 % continuo, curvas incluidas: el desnivel total no puede ser menor que 2 % del largo, se lo arme como se lo arme. Con 114 m serían 2,3 m. Para no gastar más altura que eso, los tramos alternan entre dos planos verticales y cada curva pasa de uno al otro: una hélice aplastada. (La espiral de v1.0 subía 1,6 %, menos de lo que pide CFIA.) Tercero, la aislación: con coquilla a lo largo, el tubo pierde más que con la camisa en bloque, y como el agua caliente trabaja cerca del tope de cada producto, eso cuesta placas de calentamiento. Con DN50 y las curvas también aisladas, la pérdida queda en 0,75 °C y el calentamiento vuelve a 47 placas; con las curvas al aire (5 % del largo, un cuarto de la pérdida) serían 53.

**El tubo no se achica cambiando el diámetro, pero sí se acorta.** En laminar, el volumen lo fijan el caudal, el tiempo y el factor de la partícula más rápida; el DN solo lo reparte entre largo y Δp. Por eso todos los DN dan 130 L. Pero en horquillas el largo es altura (2 % de él), así que el DN pesa. El más chico que cumple el Δp es DN40: 114 m, 19 tramos y 2,3 m de desnivel. Con DN50 quedan 66 m, 11 tramos y 1,3 m, con Δp 0,09 bar y menos pérdida de calor, a cambio de pedir 10,6 m³/h de CIP en vez de 6,1. Se eligió DN50 (`DN_fijo: 50`; vacío vuelve al criterio del DN más chico).

**Validación.** El mismo código reproduce el pasteurizador de yema publicado por Gut et al. (2003), Tabla 7, con diferencias de 0,02 °C o menos en las cinco temperaturas que reporta el paper: la yema a la salida de la regeneración (T2), del calentamiento (T3), del tubo de retención (T4), de la regeneración del lado pasteurizado (T5) y del enfriamiento (T6). Tolerancia: 0,5 °C. Eso valida balances, propiedades y encadenamiento. Las correlaciones de placa se controlan contra valores de referencia (constantes de Kumar y Martin), y el U de Gut no se usa porque no sale de sus propios datos. `python correr.py --validar` muestra las dos tablas.

## Cómo correrlo

Requiere Python 3.10 o superior.

```bash
python -m venv .venv
source .venv/Scripts/activate     # Windows (Git Bash). En Mac/Linux: source .venv/bin/activate
pip install -e ".[dev]"

pytest                            # 100 tests, incluida la validación
python correr.py                  # diseña con config/base.yaml → resultados/reporte.md
python correr.py --validar        # compara con Gut et al. (2003)
```

Para probar otro caso (otro caudal, otro largo de tramo, otro DN, el tubo en espiral…), copiá la configuración base, editá la copia y corré esa:

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
  holding.py                tubo de retención: Metzner–Reed, partícula más rápida, horquillas o espiral, elección de DN
  aislacion.py              pérdida de calor del tubo: coquilla (con o sin las curvas) o camisa en bloque
  dimensiones.py            paquete de placas y bastidor, para el plano
  seccion_phe.py            una sección de placas en contracorriente (ε-NTU), diseño y verificación
  placas.py                 placa chevron: h (Kumar, Martin), Δp de canales y puertos, U, caudal de CIP
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
| v1 | U por correlación de placa chevron (Kumar 1984; Martin 1996 de contraste) para cada producto y sección, arreglo de pasos elegido por Δp, Δp en placas, chequeo de CIP, chequeo del U de Gut | Hecho (1.0.0) |
| v1.1 | Tubo en horquillas con subida de 2 % (CFIA) y curvas por Hooper (1981); coquilla con o sin las curvas; dimensiones para el plano; DN fijable; comparación de un producto contra el agua | Hecho (1.1.0) |
| v2 | Perfil axial en las placas y letalidad acumulada (F) de placas más tubo | — |
| Después | Sensibilidad (tornado → Sobol), optimización de costo anual, Monte Carlo | — |

## Supuestos y límites

- **Placa típica, no de catálogo:** 0,15 m², 0,20 m de ancho, b = 3 mm, chevron de 45°, φ = 1,17, puertos de 50 mm. A 45° coinciden las dos convenciones de ángulo que circulan en la literatura.
- **Correlación de Kumar (1984)** para h y fricción, con las constantes que reproducen Kakaç et al. (2012). Martin (1996) se informa como contraste, pero está ajustada para Re ≥ 200 y los productos andan por Re 1–90.
- **Ensuciamiento de diseño:** 3·10⁻⁴ m²K/W por lado de producto en regeneración, 9·10⁻⁴ en calentamiento (leche, Lalande et al., 1979) y 3,4·10⁻⁵ en el agua caliente (Marriott, 1971), todos vía Gut et al. (2003). Pesa entre 20 y 57 % de la resistencia: sin él, el equipo pasa de 62 a 34 placas.
- **Yema y huevo entero en el canal:** viscosidad generalizada de Delplace y Leuliet (1995) con los parámetros de placas paralelas, como Gut. El canal chevron no es un par de placas paralelas: es una aproximación.
- **Arreglo:** mismos pasos y canales de los dos lados, en contracorriente, tratado como contracorriente pura (en Gut, ese modelo difiere ≤ 0,8 °C del distribuido). Δp máximo de 1 bar por sección y por corriente.
- **CIP de placas:** Re ≥ 2150 en los canales (el caudal mínimo efectivo de Galiński et al., 2025, trasladado como Re) con la solución a 70 °C.
- **Tubo en transición** (2100 < Re < 4000 en recto; en espiral, la banda se corre hacia arriba): se diseña con el factor laminar de la partícula más rápida, que es el lado seguro.
- **Horquillas:** partícula más rápida como en tubo recto (las curvas mezclan algo y no se cuenta). Curvas DIN 11852 (R = 70 mm en DN50, 60 mm en DN40); si los tubos aislados no entran a 2R entre planos, la curva lleva un niple. Δp de curvas por el método 2-K de Hooper (1981), pensado para newtonianos, con el Re de Metzner–Reed. Subida mínima 2 %, curvas incluidas (CFIA). Largo máximo de tramo 6 m (largo de fabricación de los tubos).
- **Espiral (alternativa):** las correcciones de fricción de White (1929) e Ito (1959) y el t_min/t_medio = 0,586 son resultados para newtonianos. El diámetro de la hélice (1 m) es un supuesto.
- **Aislación:** lana mineral de 50 mm, k = 0,040 W/m·K, sala a 20 °C. Las curvas al aire pierden como caño pelado (h = 10 W/m²K, sin la película interna: lado seguro), con 50 mm sin aislar a cada lado de la curva para la abrazadera (solo si las curvas no se aíslan). En la camisa en bloque (solo espiral) se supone el aire interior a la T del producto.
- **Dimensiones del bastidor:** las da el proveedor. Se toma de referencia un bastidor comercial con la misma placa (Alfa Laval M6: 0,92 × 0,32 m, 0,52–1,43 m de largo); el paquete de placas sí sale del modelo.
- **Margen de control:** el control trabaja 0,5 °C sobre la T legal (supuesto: el doble de la histéresis mínima de la válvula de desvío).
- **Δp del tubo con factor 1,4:** cubre la mayor diferencia de K de la yema entre la ecuación y la tabla de Gut et al. (2003).
- **Letalidad solo en el tubo.** Las placas también matan; Gut encontró que suben el F de 3,5 a 5,5 min. Eso entra en v2.
- **Reologías extrapoladas** a la temperatura de retención en los tres productos. El reporte lo avisa.
- **Espesamiento:** huevo entero (desde 56 °C) y clara a pH 9 (desde 56,7 °C) pasan por el tubo arriba de la T donde empiezan a espesar. Las reologías medidas no llegan ahí: el Δp real del tubo puede ser algo mayor. El reporte lo avisa.
- **Clara a pH 9,3:** la forma de la viscosidad con T es la de Scalzo para clara a pH 7; se supone igual a pH 9,3. El sodio de la soda (~0,4 g/kg) no se agrega a la composición de Choi–Okos.
- **Desvíos del RBN declarados:** la clara se pasteuriza a 56,7 °C × 4,96 min y no a 60 °C × 3,5 min (15.2.11), y lleva soda, que no está en la lista de aditivos para huevo. La justificación está en `docs/justificacion_clara.tex`.

## Fuentes

Las claves entre corchetes en `config/base.yaml` remiten a `docs/fuentes.bib`. Las principales: FSIS (2026) y el Reglamento Bromatológico Nacional para los pares T–t; Choi y Okos (1986) para propiedades térmicas; Kumar (1984) vía Kakaç et al. (2012) y Martin (1996) para las placas, con Delplace y Leuliet (1995) vía Gut et al. (2003) para los fluidos ley de potencia; Galiński et al. (2025) para el Δp de puertos y el CIP; Gut et al. (2003), León Munguía, Thapon y Scalzo et al. (1970) para reología; Metzner y Reed (1955) para el tubo; Ruthven (1971), Trivedi y Vasudeva (1974) y Gutierrez et al. (2010) para la partícula más rápida en la hélice; Hooper (1981) para las curvas y el catálogo DIN de Alfa Laval para sus radios; FSIS (2026, Tabla III.A) y el manual AEB para la clara a pH 9,3; Health Canada (2013), 21 CFR 184.1763 y el Codex GSFA para la soda como regulador de pH; CFIA para la operación del tubo de retención (bomba de tiempo, válvula de desvío, aislación desmontable, subida de 2 %).

---

Proyecto Industrial de Ingeniería Química, Facultad de Ingeniería, Udelar (2026): *Producción de huevo en polvo y valorización de la cáscara*.
