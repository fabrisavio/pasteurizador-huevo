# Pseudocódigo: modelo paramétrico del pasteurizador (v0.7: v0.5 implementada el 08/11, cierre de v0)

> **Cambios al implementar v0** (respecto de la v0.2 aprobada):
> - Agua con el paquete `iapws` (IAPWS-IF97, Python puro) en vez de CoolProp (IAPWS-95): en agua líquida difieren muy poco y se instala en cualquier compu sin compilar.
> - Tubo: se agrega `Re_turbulento: 4000`. Entre 2100 y 4000 (transición) se diseña con el factor laminar de la partícula más rápida, que es el lado seguro.
> - Regla de elección del DN: el más chico que cumple el Δp máximo (minimiza el caudal de CIP; en laminar el volumen no cambia). El reporte avisa si el margen es menor al 10 %.
> - La verificación corre cada producto en los dos extremos de su rango de T de entrada.
> - El reporte separa **alertas** (algo no cumple) de **avisos** (supuestos, extrapolaciones, datos para otros equipos).
> - v0 ya incluye los tres productos y Metzner–Reed en el tubo; quedan para v1 las correlaciones de h y f, el Δp de placas y el chequeo de CIP.
> - **v0.2 (31/10):** opción de tubo enrollado (`geometria: espiral`). Sube la fricción con White (1929) en laminar e Ito (1959) en turbulento, corre el Re crítico a 2·10⁴·(D/Dc)^0,32 y calcula vueltas y altura. La config base pasa a espiral de 1 m, y con eso el DN elegido pasa de 32 a 40.
> - **v0.3 (02/11):** la partícula más rápida en la hélice se toma de literatura: t_min/t_medio = 0,586, el extremo bajo de Ruthven (1971) y Trivedi y Vasudeva (1974), revisados por Gutierrez et al. (2010), que midieron ≥ 0,70. v_max/v_media = min(1/0,586; (3n+1)/(n+1)) si Dean > 11,6; si no, el de tubo recto. El tubo pasa de 106 m a 94 m. El reporte compara con el mismo DN en recto. Además: recorrido de T del producto en consola y reporte, puntos de la validación con nombre y esquema, la clara se documenta como estabilizada con láctico y aluminio (por eso soporta 60 °C), y mensaje claro si la config no existe.
> - **v0.4 (02/11):** la pérdida de T en el tubo se calcula: el tubo es un intercambiador contra el ambiente (Cr = 0), ΔT = (T_sale − T_amb)·(exp(UA/C) − 1). Dos formas de aislar: coquilla a lo largo del tubo, o camisa desmontable alrededor de toda la espiral (bloque). Se itera tubo ↔ pérdida por punto fijo. El reporte compara las dos, incluidas las placas de calentamiento que pide cada una (bloque 24, coquilla 27). Margen de control de 0,5 °C sobre la T legal: el tubo sale a T_legal + margen y el calentamiento entrega eso más la pérdida; los log se siguen contando a la T legal. Factor 1,4 sobre el Δp del tubo por la incertidumbre de K de la yema (sin él DN32 entraba con 1 % de margen).
> - **v0.5 (08/11, cierre de v0):** ruta de la soda para la clara. Al terminar el desazucarado se lleva a pH 9,3 con NaOH diluida y se pasteuriza a 56,7 °C × 4,96 min (FSIS, Tabla III.A, pH 9,3, 5,7 log). Cinética de CP a pH 9,3 por ajuste log-lineal de esa tabla: D₅₆,₇ = 0,854 min, z = 4,17 °C. Agua caliente de CP ≤ 59 °C (Cunningham, p. 296). Nuevos campos `pH` (informativo) y `T_espesa_C` por producto: si el tubo supera la T de espesamiento, el reporte da un aviso (la reología no lo ve). El reporte muestra pH y tratamiento legal por producto. Con esto CP pasa a gobernar regeneración y tubo. Justificación lista para el informe en `docs/justificacion_clara.tex`. Reemplaza la estabilización con alumbre documentada en v0.3.

Pasteurizador de placas (PHE) con regeneración y calentamiento con agua caliente, más tubo de retención (holding tube). Sin sección de enfriamiento, porque los tres productos salen a 45 °C directo al secador.

Unidades: el archivo de configuración usa unidades prácticas (kg/h, °C, min). El código convierte todo a SI (kg/s, K, Pa, W) al cargarlo y vuelve a unidades prácticas solo en el reporte.

---

## 0. El flujo completo en una línea

```
config.yaml → cargar → [por producto] diseñar → equipo único (máximo por sección) → [por producto] verificar en ese equipo → reporte
```

**Por qué en dos pasadas:** el mismo equipo procesa los tres productos y cada sección la gobierna uno distinto. Primero se calcula qué necesitaría cada producto por separado. Después se arma un equipo que cumpla con todos. Al final se verifica cómo funciona cada producto en ese equipo real.

---

## 1. Archivo de configuración (inputs, nada de números sueltos en el código)

```yaml
productos:
  HP:
    caudal_kg_h: 798
    T_entrada_C: {diseño: 10, rango: [10, 25]}
    T_holding_C: 60.0
    t_holding_min: 3.5
    T_salida_C: 45
    T_max_agua_C: 63
    composicion: {agua: 0.7651, proteina: 0.126, grasa: 0.10, carbohidratos: 0.003, cenizas: 0.0059}  # BMyBCrafa, Tabla 2
    reologia: leon_munguia            # K'(T), n' = 0.955
    cinetica: {D_ref_min: 1.25, T_ref_C: 60.0, z_C: 4.84}  # derivado de FSIS 2026, Ap. III
    horas_anio: 2496
  YP:
    caudal_kg_h: 1083
    T_entrada_C: {diseño: 7, rango: [7, 7.2]}   # tanque a 7 °C (FSIS: ≤ 7,2 °C, < 10 % de sal, espera ≤ 8 h)
    T_holding_C: 61.1
    t_holding_min: 3.5
    T_salida_C: 45
    T_max_agua_C: 64
    composicion: {agua: 0.562, proteina: 0.155, grasa: 0.256, carbohidratos: 0.0115, cenizas: 0.0155}  # yema comercial USDA
    reologia: gut2003                 # tablas K(T), n(T)
    cinetica: {D_ref_min: 1.68, T_ref_C: 61.1, z_C: 7.31}
    horas_anio: 168
  CP:
    caudal_kg_h: 822.4
    T_entrada_C: {diseño: 35, rango: [35, 35]}
    pH: 9.3                           # ruta de la soda (08/11): NaOH al final del desazucarado
    T_holding_C: 56.7                 # FSIS Tabla III.A, pH 9,3 (5,7 log)
    t_holding_min: 4.96
    T_salida_C: 45
    T_max_agua_C: 59                  # Cunningham (1995), p. 296
    T_espesa_C: 56.7                  # pp. 290–291: solo genera un aviso
    composicion: {agua: 0.77, proteina: 0.203, grasa: 0.0006, carbohidratos: 0.0148, cenizas: 0.0116}  # clara UF 23 %
    reologia: thapon                  # newtoniana, μ(T)
    cinetica: {D_ref_min: 0.854, T_ref_C: 56.7, z_C: 4.17}  # ajuste de FSIS III.A a pH 9,3
    horas_anio: 216

equipo:
  placa:                              # placa chevron "típica": cada dato dentro del rango de libro (PENDIENTE: fijar fuente)
    area_m2: 0.15                     # rango de fabricación 0,03–3,6 m²; para ~1 m³/h, placas chicas (Gut: 0,03 m² a 230 L/h)
    ancho_m: 0.20
    separacion_b_m: 0.003             # rango 1,5–5 mm
    angulo_chevron_grados: 45         # rango 25–65°: bajo = "blanda" (menos h y Δp), alto = "dura"
    factor_area: 1.17                 # rango 1,1–1,25
    espesor_m: 0.0005
    k_W_mK: 16.3                      # AISI 316L
  correlacion_placas: kumar1984       # base (constantes según ángulo y Re, cubre Re bajos); martin1996 como contraste
  ensuciamiento_m2K_W: {huevo_calentamiento: 9.0e-4, huevo_otras: 3.0e-4, agua: 3.4e-5}  # Lalande (vía Gut); se usan como rango
  holding:
    DN_diam_int_m: {25: 0.026, 32: 0.032, 40: 0.038, 50: 0.050, 65: 0.066}   # DIN 11850, serie 2
    perdida_T_C: 0.5
    dp_max_bar: 1.0                   # cota inferior del DN: bomba de desplazamiento positivo y presión de las juntas
    v_CIP_min_m_s: 1.5                # cota superior del DN: velocidad mínima de limpieza en cañería (regla práctica; citar)
servicios:
  agua_caliente: {DT_C: 4}            # variable de decisión v1
criterios:
  Re_CIP_min: 1500                    # Galiński et al. (2025)
  factor_turbulento_vmax: 1.25        # supuesto (conocimiento general)
  Re_critico: 2100
  logs_objetivo: 5
```

**Qué hace:** junta todos los inputs con su fuente (las claves coinciden con el registro).
**Por qué:** si mañana cambia un caudal o una temperatura, se cambia acá y el diseño se recalcula solo.

---

## 2. Propiedades (`propiedades.py`)

```
función propiedades(producto, T):
    # Choi–Okos: la misma fuente que el secador de Rafa
    para cada componente i en {agua, proteína, grasa, carbohidratos, cenizas}:
        cp_i = A + B·T + C·T²           # coeficientes = Tabla 11 de Rafa
        ρ_i  = polinomio(T)             # misma fuente (Rafa los confirma contra su libro)
        k_i  = polinomio(T)
    cp  = Σ x_i·cp_i                    # promedio másico
    ρ   = 1 / Σ (x_i/ρ_i)               # volúmenes aditivos
    φ_i = (x_i/ρ_i)·ρ                   # fracción volumétrica
    k   = Σ φ_i·k_i                     # modelo en paralelo (el de Choi–Okos)
    (K, n) = reologia(producto, T)
    devolver ρ, cp, k, K, n

función reologia(producto, T):
    HP: K = 0,1·(0,0689 − 0,0007·T) [Pa·sⁿ];  n = 0,955
        si T fuera de 10–40 °C → AVISO "extrapolación" (no ve la subida de viscosidad sobre 56 °C)
    YP: K, n = interpolar en las tablas de Gut (0,4–60,8 °C)
        (no la ec. 4 sobre 48 °C: subestima K en 23–38 %)
    CP: n = 1;  K = μ(T) = 9,3 mPa·s × [μ_Scalzo(T) / μ_Scalzo(30 °C)]
        # 9,3 = Thapon interpolado a 23 % de sólidos a 30 °C; la forma con T es de Scalzo (supuesto declarado)

función propiedades_agua(T):  usar iapws (IAPWS-IF97)
```

**Por qué así:** térmicas por composición, que es coherente con el secador y responde si cambia el % de sólidos. Reología medida, porque Choi–Okos no la da.

---

## 3. Cinética (`cinetica.py`) — verifica, no dimensiona

```
D(T)   = D_ref · 10^((T_ref − T)/z)
logs(t, T) = t / D(T)
```

**Por qué no dimensiona:** el diseño se basa en el par T–t de la norma (safe harbor o RBN). La cinética solo informa cuántos logaritmos se logran.

---

## 4. Tubo de retención (`holding.py`)

```
función holding(producto, DN, Q):
    A = π·Di²/4 ;  v_media = Q/A
    # Reynolds generalizado de Metzner–Reed (n = 1 → Re de siempre)
    μ_ef  = K·((3n+1)/(4n))ⁿ · (8·v_media/Di)^(n−1)
    Re_MR = ρ·v_media·Di / μ_ef
    si Re_MR < Re_critico:   factor = (3n+1)/(n+1)           # laminar: partícula más rápida (2 si n = 1)
                              f_Fanning = 16/Re_MR
    si no:                    factor = factor_turbulento_vmax
                              f_Fanning = Blasius o Colebrook
    v_max   = factor·v_media
    L       = v_max · t_holding                               # "every particle" (9 CFR 590.570)
    Δp      = 4·f_Fanning·(L/Di)·ρ·v_media²/2
    t_medio = L / v_media                                     # cuánto se sobreprocesa la partícula promedio
    devolver L, Δp, volumen, Re_MR, régimen, t_medio,
             logs_min = logs(t_holding, T_holding), logs_medio = logs(t_medio, T_holding)

función elegir_tubo(productos):
    para cada DN disponible:
        L_DN = máx sobre productos de holding(p, DN, Q_p).L   # un solo tubo físico; casi seguro lo gobierna YP
        Δp_p y t_medio_p de cada producto con ese L
        caudal_CIP_DN = v_CIP_min · área(DN)
    descartar DN con Δp_p > dp_max (DN chicos: Δp ∝ 1/D⁶ en laminar)
    entre los que quedan, elegir el que minimiza el largo sin disparar el caudal de CIP
    informar el sobreproceso de HP y CP (opción a evaluar: tramos de tubo seleccionables por producto)
```

**Qué sale:** diámetro, largo, Δp y cuánto tiempo pasa cada producto. Con el tubo dimensionado por el peor caso, los otros productos se quedan más de 3,5 min.

**Dato útil:** en laminar, el volumen del tubo no depende del diámetro (V = factor × caudal × t). El diámetro solo reparte ese volumen entre largo y Δp.

---

## 5. Sección de PHE genérica (`seccion_phe.py`)

Un solo bloque sirve para la regeneración (producto–producto), el calentamiento (agua–producto) y, para la validación, el enfriamiento.

```
función h_y_f(corriente, geometría, N_canales_por_paso, T_media):
    v  = m / (ρ·N·b·w) ;  De = 2b/φ
    Re = Reynolds generalizado para canal (Delplace & Leuliet 1995; con n = 1 queda el Re de siempre)
    Pr = cp·μ_ef/k
    Nu, f = correlación chevron: Kumar (1984), constantes según ángulo y rango de Re (reproducida en Kakaç & Liu)
            ; Martin (1996, VDI) como contraste en la sensibilidad
    h = Nu·k/De
    devolver h, f, Re

función diseñar_seccion(caliente, fria, geometría, Rf):
    # T de entrada y salida conocidas del lado producto; el otro lado sale por balance
    Q  = m_c·cp_c·(Tc_in − Tc_out)       # = m_f·cp_f·(Tf_out − Tf_in)
    ε  = Q / (C_min·(Tc_in − Tf_in)) ;  C_r = C_min/C_max
    NTU = inversa_contracorriente(ε, C_r)          # C_r ≈ 1 → NTU = ε/(1−ε)
    elegir N_canales_por_paso (empezar en 1):
        h_c, h_f ← h_y_f(...)
        U = 1 / (1/h_c + e/k_placa + 1/h_f + Rf_c + Rf_f)
        A = NTU·C_min/U ;  N_placas = techo(A/A_placa) + 1
        repetir hasta que N_placas y N_canales sean coherentes
    Δp = canales (f) + puertos
    devolver A, N_placas, U, h_c, h_f, Re, Δp, Q

función verificar_seccion(A, geometría, caliente(m, T_in), fria(m, T_in)):
    adivinar T de salida → propiedades a T media → U → NTU = U·A/C_min
    ε = ε_contracorriente(NTU, C_r) ;  Q = ε·C_min·(Tc_in − Tf_in) ;  salidas por balance
    repetir hasta que las T de salida no cambien (punto fijo o fsolve)
```

**Modelo de U:** v1 usa U constante por sección, con propiedades a T media. Es lo que recomienda Gut: casi no cambia las temperaturas. v2 pasa a perfil axial (ODE), que hace falta para calcular bien la letalidad acumulada en las placas.

---

## 6. Un producto, modo diseño (`pasteurizador.py`)

```
función diseñar_producto(p):
    T_salida_calent = T_holding + perdida_T          # el sensor y la válvula de desvío están a la SALIDA del tubo

    # Regeneración: producto pasteurizado (caliente) contra producto crudo (frío)
    caliente: T_holding → T_salida (45 °C)
    fria:     T_entrada → T_salida_regen (por balance; ≈ T_entrada + 15 si los cp son iguales)
    regen = diseñar_seccion(caliente, fria)

    # Calentamiento: agua caliente contra producto
    producto: T_salida_regen → T_salida_calent
    agua:     T_entrada = T_max_agua (garantiza T_pared ≤ T_max) ; T_salida = T_max_agua − DT_agua
    calent = diseñar_seccion(agua, producto)          # el caudal de agua sale por balance

    tubo = holding(p, DN, Q_p)

    servicios: Q_calent [kW], caudal de agua caliente, potencia de bombas = Q·Δp/η
    chequeos:  logs_min ≥ logs_objetivo (informa), approach > 0, Re dentro de la validez de las correlaciones
```

---

## 7. Un solo equipo para los tres (multi-producto)

```
para p en {HP, YP, CP}:  diseño[p] = diseñar_producto(p)

equipo.N_placas[regen]  = máx_p diseño[p].N_placas[regen]
equipo.N_placas[calent] = máx_p diseño[p].N_placas[calent]
equipo.tubo             = elegir_tubo({HP, YP, CP})

para p:  verif[p] = verificar_producto(equipo, p)
    # con más placas de las que necesita: ¿a qué T sale?, ¿qué T de agua precisa para llegar a T_holding?
    chequeos: T_holding alcanzada con T_agua ≤ T_max ; T_salida ≤ 45 °C (más frío está bien) ;
              viscosidad de la yema en el atomizador si sale muy fría ; t_residencia en el tubo ; Q_calent real
```

**Lo que ya se ve venir:** la regeneración la va a gobernar CP, porque entra a 35 °C y necesita ε ≈ 0,6 contra ~0,3 de HP e YP. Con esa regeneración, HP e YP regeneran de más y salen más fríos que 45 °C. **Decisión del 31/10: se acepta.** El único chequeo que queda es la viscosidad de la yema en el atomizador.

---

## 8. CIP (chequeo)

```
caudal_CIP = mínimo caudal de agua tal que Re_canal ≥ Re_CIP_min en todas las secciones del equipo final
```

---

## 9. Energía (para el peldaño 4)

```
para p:  E_calor[p] = Q_calent[p]·horas_anio[p] ;  E_bombas[p] = P_bombas[p]·horas_anio[p]
totales anuales ;  (después: costo anual = CAPEX(área, tubo, bombas) + OPEX(calor, electricidad))
```

---

## 10. Validación (peldaño 2): caso Gut et al. (2003)

```
config/validacion_gut2003.yaml:
    yema 230 L/h a 5 °C ; agua caliente 0,15 kg/s a 66 °C ; agua fría 0,15 kg/s a 2 °C
    placas planas: L 0,33 m ; w 0,09 m ; b 5 mm ; e 0,5 mm ; NC = 32 / 13 / 13
    Nu = 1,35·Re^0,36·Pr^0,33 (Sharifi) ; Rf 9e-4 (calentamiento), 3e-4 (resto) ; propiedades de Gut
correr verificar_seccion encadenando regeneración + calentamiento + enfriamiento (con enfriamiento, como Gut)
comparar T2…T6 y U con la Tabla 7 (columna de contracorriente pura)
criterio de aceptación: |ΔT| < 0,5 °C y |ΔU| < 5 %
```

**Por qué con enfriamiento:** el bloque de sección es genérico. Agregar la tercera sección es solo configuración, y así el código se valida contra un caso publicado.

---

## 11. Tests automáticos (`tests/`)

```
balance de cada sección: Q_caliente = Q_fria
ε-NTU: NTU → ∞ ⇒ ε → 1 ;  C_r = 1 ⇒ NTU = ε/(1−ε) ;  verificar(diseñar(x)) devuelve x
Choi–Okos: agua pura ≈ 998 kg/m³ y 4,18 kJ/(kg·K) a 20 °C
Metzner–Reed: con n = 1 da ρ·v·D/μ
holding: con n = 1 laminar el factor da 2
caso Gut dentro de tolerancia
```

---

## 12. Salidas (reporte)

- Por sección: placas, área, U, h, Re, Δp, Q.
- Tubo de retención: DN, L, volumen, Δp, régimen y residencias por producto.
- Por producto: T_salida real, T de agua caliente necesaria, caudal de agua, logaritmos logrados, chequeos (OK/aviso).
- Servicios: kW de calor, kW de bombas, kWh/año.
- Avisos: extrapolaciones de propiedades, correlaciones fuera de rango.

---

## 13. Versiones

| Versión | Qué agrega |
|---|---|
| v0 | U fijo de Gut (regen. 756 / calent. 864 W/m²K), un producto, tubo con factor 2: corre de punta a punta |
| v1 | Correlaciones de h y f, los tres productos, diseño + verificación, tubo con Metzner–Reed, CIP |
| v2 | Perfil axial en las placas (ODE) y letalidad acumulada (F), energía anual |
| después | Sensibilidad (tornado → Sobol), optimización, Monte Carlo |

---

## 14. Estructura del repositorio

```
pasteurizador/
├── README.md                    # abre con el problema de proceso
├── config/
│   ├── base.yaml
│   └── validacion_gut2003.yaml
├── src/pasteurizador/
│   ├── propiedades.py
│   ├── cinetica.py
│   ├── holding.py
│   ├── seccion_phe.py
│   ├── pasteurizador.py
│   └── reporte.py
├── tests/
├── docs/                        # este pseudocódigo + el registro de fuentes
└── correr.py                    # python correr.py config/base.yaml
```

Librerías: numpy, scipy (brentq/fsolve, interpolación, solve_ivp), pyyaml, iapws, pytest (pandas y matplotlib cuando hagan falta gráficos).
