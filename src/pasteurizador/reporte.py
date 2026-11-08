"""Salidas: reporte en Markdown y tablas CSV.

Markdown porque se lee bien en GitHub y en VS Code (Ctrl+Shift+V abre la
vista previa) y se pega fácil en el informe. Los CSV son para abrir en Excel
o levantar con pandas.
"""

from __future__ import annotations

import csv
from datetime import date
from pathlib import Path

from . import __version__
from .pasteurizador import Resultado


def _tabla(encabezado: list[str], filas: list[list]) -> str:
    """Arma una tabla Markdown. Los números chicos van en notación científica."""
    def fmt(x):
        if isinstance(x, float):
            return f"{x:.3g}" if abs(x) < 0.01 and x != 0 else f"{x:,.2f}".replace(",", " ")
        return str(x)
    lineas = ["| " + " | ".join(encabezado) + " |",
              "|" + "|".join("---" for _ in encabezado) + "|"]
    lineas += ["| " + " | ".join(fmt(x) for x in f) + " |" for f in filas]
    return "\n".join(lineas)


def describir_tubo(r: Resultado) -> str:
    """Una frase con la forma del tubo elegido (la usa también la consola)."""
    t, geo = r.equipo.tubo.elegido, r.equipo.tubo.geometria
    if not geo.es_espiral:
        return f"DN{t.DN} (Ø int {t.D*1000:.0f} mm), {t.L:.1f} m en tramos rectos"
    forma = (f"DN{t.DN} (Ø int {t.D*1000:.0f} mm), {t.L:.1f} m enrollados en espiral de "
             f"Ø {geo.D_espiral:.2f} m: {t.vueltas:.1f} vueltas")
    if t.alto is not None:
        forma += f", {t.alto:.2f} m de alto"
    return forma


NOMBRE_AISLACION = {"tubo": "Coquilla a lo largo del tubo", "bloque": "Camisa en bloque"}


def aislacion_md(r: Resultado) -> str:
    """Comparación de las dos formas de aislar el tubo y lo que ahorra la elegida."""
    cfg, eq = r.cfg, r.equipo
    a, ops = cfg.aislacion, eq.aislacion
    s = [f"\n### Aislación del tubo\n",
         f"Lana mineral de {a.espesor*1000:.0f} mm (k = {a.k} W/m·K), sala a {a.T_amb:.0f} °C. "
         f"Diseño con: **{NOMBRE_AISLACION[a.forma].lower()}**.\n"]
    claves = list(cfg.productos)
    placas = eq.placas_calent_por_aislacion
    encab = (["Forma", "Superficie exterior (m²)", "Aislante (m³)", "UA (W/K)"]
             + [f"Pérdida {k} (°C)" for k in claves]
             + ["Calor perdido (MWh/año)", "Placas de calentamiento"])
    filas = [[NOMBRE_AISLACION[o.forma], o.area_ext, o.volumen_aislante, o.UA]
             + [o.dT[k] for k in claves] + [o.E_MWh_anio, placas.get(o.forma, "—")]
             for o in ops.values()]
    s.append(_tabla(encab, filas))
    if "bloque" in ops and "tubo" in ops:
        b, t = ops["bloque"], ops["tubo"]
        peor = max(claves, key=lambda k: t.dT[k] - b.dT[k])
        s.append(f"\n**Lo que ahorra la camisa en bloque frente a la coquilla:** "
                 f"{t.area_ext - b.area_ext:.1f} m² menos de revestimiento ({1 - b.area_ext/t.area_ext:.0%}), "
                 f"{t.volumen_aislante - b.volumen_aislante:.2f} m³ menos de aislante "
                 f"({1 - b.volumen_aislante/t.volumen_aislante:.0%}) y "
                 f"{(t.E_MWh_anio - b.E_MWh_anio)*1000:.0f} kWh/año menos de calor perdido. "
                 f"El producto llega {t.dT[peor] - b.dT[peor]:.2f} °C menos caliente al tubo "
                 f"({peor}), y como el agua caliente trabaja cerca de su tope, eso ahorra "
                 f"{placas['tubo'] - placas['bloque']} placas de calentamiento "
                 f"({placas['tubo']} → {placas['bloque']}).\n")
    return "\n".join(s)


def markdown(r: Resultado, origen: str = "") -> str:
    cfg, eq = r.cfg, r.equipo
    claves = list(cfg.productos)
    t = eq.tubo.elegido
    s = [f"# Pasteurizador: resultados v{__version__}\n",
         f"Config: `{origen}` · generado el {date.today().isoformat()}\n"]

    # Equipo
    s.append("## Equipo\n")
    filas = []
    for nombre, n, A, U, gob, attr in [
        ("Regeneración", eq.placas_regen, eq.A_regen, cfg.U_regen, eq.gobierna_regen, "regen"),
        ("Calentamiento", eq.placas_calent, eq.A_calent, cfg.U_calent, eq.gobierna_calent, "calent"),
    ]:
        pedidas = []
        for k in claves:
            sec = getattr(eq.disenos[k], attr)  # getattr(obj, "regen") = obj.regen
            pedidas.append(f"{k} {sec.A:.2f}" if sec else f"{k} —")
        filas.append([nombre, n, A, U, gob, " · ".join(pedidas)])
    s.append(_tabla(["Sección", "Placas térmicas", "Área instalada (m²)", "U (W/m²K)",
                     "Gobierna", "Área que pide cada producto (m²)"], filas))
    s.append(f"\nPlaca de {cfg.area_placa} m². Bomba de agua caliente: "
             f"{eq.m_agua:.2f} kg/s ({eq.m_agua*3.6:.1f} m³/h aprox.).\n")

    # Tubo
    margen = 1 - t.dp_max / cfg.dp_max_holding
    geo = eq.tubo.geometria
    s.append("## Tubo de retención\n")
    s.append(f"**Elegido: {describir_tubo(r)}.** {t.volumen_L:.0f} L, "
             f"Δp {t.dp_max/1e5:.2f} bar (margen {margen:.0%} contra {cfg.dp_max_holding/1e5:.1f} bar), "
             f"CIP a {cfg.v_CIP} m/s: {t.Q_CIP_m3_h:.1f} m³/h. El largo lo fija {t.gobierna}.\n")
    if geo.es_espiral:
        s.append(f"Partícula más rápida en la hélice: t_min/t_medio = {geo.t_min_sobre_t_medio}, "
                 f"o sea v_max/v_media = {geo.factor_espiral:.2f} (Ruthven 1971; Trivedi y "
                 f"Vasudeva 1974; vía Gutierrez et al. 2010). Si un producto tiene Dean ≤ "
                 f"{geo.Dean_min} se usa el de tubo recto.\n")
    encab = ["DN", "D int (mm)", "Largo (m)", "Volumen (L)", "Δp máx (bar)", "Q CIP (m³/h)"]
    if geo.es_espiral:
        encab += ["Vueltas", "Alto (m)"]
    encab.append("Cumple Δp")
    filas = []
    for op in eq.tubo.opciones:
        # Δp con 3 decimales: si no, un 1,004 bar se ve "1.00" y parece que cumple
        fila = [op.DN, op.D * 1000, op.L, op.volumen_L, f"{op.dp_max / 1e5:.3f}", op.Q_CIP_m3_h]
        if geo.es_espiral:
            fila += [op.vueltas, op.alto if op.alto is not None else "—"]
        fila.append("sí" if op.factible else "no")
        filas.append(fila)
    s.append(_tabla(encab, filas))

    s.append(f"\nCada producto en el DN{t.DN}:\n")
    encab = ["Producto", "Régimen", "Re (M–R)"]
    if geo.es_espiral:
        encab += ["Dean", "f espiral / f recto"]
    encab += ["v_max/v_media", "t partícula rápida (min)", "t medio (min)",
              "log rápida", "log media", "Δp (bar)"]
    filas = []
    for k, x in t.productos.items():
        fila = [k, x.flujo.regimen.replace("transicion", "transición"), round(x.flujo.Re)]
        if geo.es_espiral:
            fila += [round(x.flujo.Dean), x.flujo.correccion_f]
        fila += [x.flujo.factor, x.t_rapida_min, x.t_media_min, x.logs_min, x.logs_medio,
                 x.dp / 1e5]
        filas.append(fila)
    s.append(_tabla(encab, filas))

    rec = eq.tubo.equivalente_recto
    if rec is not None:
        f_rec = rec.productos[rec.gobierna].flujo.factor
        s.append(f"\n**Lo que ahorra enrollarlo.** El mismo DN{t.DN} en tramos rectos usaría el "
                 f"factor de tubo recto ({f_rec:.2f} para {rec.gobierna}) y necesitaría "
                 f"{rec.L:.1f} m ({rec.volumen_L:.0f} L, Δp {rec.dp_max/1e5:.2f} bar). "
                 f"En espiral: {t.L:.1f} m ({t.volumen_L:.0f} L, Δp {t.dp_max/1e5:.2f} bar), "
                 f"{1 - t.L/rec.L:.0%} menos de largo y volumen, a cambio de más fricción por vuelta.\n")

    # Aislación del tubo
    if eq.aislacion:
        s.append(aislacion_md(r))

    # Recorrido del producto
    s.append("\n## Temperaturas del producto en el equipo (T de entrada de diseño)\n")
    filas = []
    for k in claves:
        p, v = cfg.productos[k], r.en_diseno(k)
        ph = f"{p.pH:.1f}" if p.pH is not None else "—"
        legal = f"{p.T_holding:.1f} °C × {p.t_holding/60:.2f} min"
        filas.append([f"{k} ({p.nombre})", ph, legal, p.m * 3600, v.T_entrada, v.T_regen_salida,
                      v.T_calent_salida, cfg.T_operacion(p), v.T_salida])
    s.append(_tabla(["Producto", "pH", "Tratamiento legal", "Caudal (kg/h)", "Entra (°C)",
                     "Sale de regeneración (°C)", "Entra al tubo (°C)", "Sale del tubo (°C)",
                     "Al secador (°C)"], filas))
    s.append(f"\nEl control trabaja {cfg.margen_control:.1f} °C por encima de la T legal "
             "(margen de control); el tubo pierde lo que va de la entrada a la salida.\n")

    # Operación
    s.append("\n## Operación de cada producto en el equipo\n")
    filas = []
    for k in claves:
        p = cfg.productos[k]
        for v in r.verificaciones[k]:
            filas.append([k, v.T_entrada, f"{v.T_agua_necesaria:.1f} (tope {p.T_max_agua:.0f})",
                          v.T_salida, v.Q_calent_kW, v.Q_regen_kW, v.E_calor_MWh_anio,
                          f"{v.K_salida:.3g} / {v.n_salida:.2f}"])
    s.append(_tabla(["Producto", "T entrada (°C)", "T agua necesaria (°C)", "T al secador (°C)",
                     "Q calentamiento (kW)", "Q regenerado (kW)", "Calor anual (MWh)",
                     "K (Pa·sⁿ) / n a la salida"], filas))
    total = sum(r.en_diseno(k).E_calor_MWh_anio for k in claves)
    s.append(f"\nCalor anual a T de entrada de diseño: **{total:.1f} MWh**.\n")

    # Diseño individual
    s.append("## Lo que pediría cada producto por separado\n")
    filas = []
    for k in claves:
        d = eq.disenos[k]
        for sec in (d.regen, d.calent):
            if sec:
                filas.append([k, sec.nombre, sec.Tf_in, sec.Tf_out, sec.Tc_in, sec.Tc_out,
                              sec.eps, sec.NTU, sec.Cr, sec.A, sec.Q / 1000])
    s.append(_tabla(["Producto", "Sección", "Fría entra (°C)", "Fría sale (°C)",
                     "Caliente entra (°C)", "Caliente sale (°C)", "ε", "NTU", "Cr",
                     "Área (m²)", "Q (kW)"], filas))

    # Alertas y avisos
    s.append("\n## Alertas (algo no cumple)\n")
    s += [f"- {a}" for a in r.alertas] or ["Ninguna: todos los productos cumplen en todo su rango de entrada."]
    s.append("\n## Avisos (supuestos, extrapolaciones, datos para otros equipos)\n")
    s += [f"- {a}" for a in r.avisos]
    return "\n".join(s) + "\n"


def escribir(r: Resultado, carpeta: str | Path, origen: str = "") -> Path:
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)  # crea la carpeta si no existe
    (carpeta / "reporte.md").write_text(markdown(r, origen), encoding="utf-8")

    # newline="" evita líneas en blanco de más en Windows
    with open(carpeta / "operacion.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["producto", "T_entrada_C", "T_regen_salida_C", "T_calent_salida_C",
                    "T_agua_C", "T_salida_C", "Q_calent_kW", "Q_regen_kW",
                    "E_calor_MWh_anio", "K_salida_Pa_sn", "n_salida"])
        for k, vs in r.verificaciones.items():
            for v in vs:
                w.writerow([k, v.T_entrada, round(v.T_regen_salida, 3),
                            round(v.T_calent_salida, 3), round(v.T_agua_necesaria, 3),
                            round(v.T_salida, 3), round(v.Q_calent_kW, 3), round(v.Q_regen_kW, 3),
                            round(v.E_calor_MWh_anio, 3), v.K_salida, v.n_salida])

    with open(carpeta / "tubo.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["DN", "D_int_mm", "L_m", "volumen_L", "vueltas", "alto_m", "dp_max_bar",
                    "Q_CIP_m3_h", "cumple", "producto", "regimen", "Re", "Dean", "factor",
                    "t_rapida_min", "t_media_min", "logs_min", "logs_medio", "dp_bar"])
        for op in r.equipo.tubo.opciones:
            for k, x in op.productos.items():
                w.writerow([op.DN, op.D * 1000, round(op.L, 2), round(op.volumen_L, 1),
                            round(op.vueltas, 2) if op.vueltas else "",
                            round(op.alto, 3) if op.alto else "",
                            round(op.dp_max / 1e5, 4), round(op.Q_CIP_m3_h, 2), op.factible,
                            k, x.flujo.regimen, round(x.flujo.Re), round(x.flujo.Dean, 1),
                            round(x.flujo.factor, 3), round(x.t_rapida_min, 3),
                            round(x.t_media_min, 3), round(x.logs_min, 2),
                            round(x.logs_medio, 2), round(x.dp / 1e5, 4)])
    return carpeta / "reporte.md"
