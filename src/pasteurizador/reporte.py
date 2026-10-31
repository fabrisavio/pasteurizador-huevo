"""Salidas: reporte en Markdown y tablas CSV."""

from __future__ import annotations

import csv
from datetime import date
from pathlib import Path

from .pasteurizador import Resultado


def _tabla(encabezado: list[str], filas: list[list]) -> str:
    def fmt(x):
        if isinstance(x, float):
            return f"{x:.3g}" if abs(x) < 0.01 and x != 0 else f"{x:,.2f}".replace(",", " ")
        return str(x)
    lineas = ["| " + " | ".join(encabezado) + " |",
              "|" + "|".join("---" for _ in encabezado) + "|"]
    lineas += ["| " + " | ".join(fmt(x) for x in f) + " |" for f in filas]
    return "\n".join(lineas)


def markdown(r: Resultado, origen: str = "") -> str:
    cfg, eq = r.cfg, r.equipo
    claves = list(cfg.productos)
    t = eq.tubo.elegido
    s = [f"# Pasteurizador: resultados v0.2\n",
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
            sec = getattr(eq.disenos[k], attr)
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
    if geo.es_espiral:
        forma = f"enrollado en espiral de {geo.D_espiral:.2f} m de diámetro"
        if t.alto is not None:
            forma += f": {t.vueltas:.1f} vueltas, {t.alto:.2f} m de alto"
        else:
            forma += f": {t.vueltas:.1f} vueltas"
    else:
        forma = "en tramos rectos (curvas despreciadas)"
    s.append(f"**Elegido: DN{t.DN}**, {t.L:.1f} m {forma}. {t.volumen_L:.0f} L, "
             f"Δp {t.dp_max/1e5:.2f} bar (margen {margen:.0%} contra {cfg.dp_max_holding/1e5:.1f} bar), "
             f"CIP a {cfg.v_CIP} m/s: {t.Q_CIP_m3_h:.1f} m³/h. El largo lo fija {t.gobierna}.\n")
    s.append(_tabla(
        ["DN", "D int (mm)", "Largo (m)", "Volumen (L)", "Δp máx (bar)", "Q CIP (m³/h)", "Cumple Δp"],
        [[op.DN, op.D * 1000, op.L, op.volumen_L, op.dp_max / 1e5, op.Q_CIP_m3_h,
          "sí" if op.factible else "no"] for op in eq.tubo.opciones]))
    s.append(f"\nCada producto en el DN{t.DN}:\n")
    encab = ["Producto", "Régimen", "Re (M–R)", "v_max/v_media", "t partícula rápida (min)",
             "t medio (min)", "log rápida", "log media", "Δp (bar)"]
    if geo.es_espiral:
        encab += ["Dean", "f espiral / f recto"]
    filas = []
    for k, x in t.productos.items():
        fila = [k, x.flujo.regimen.replace("transicion", "transición"), round(x.flujo.Re),
                x.flujo.factor, x.t_rapida_min, x.t_media_min, x.logs_min, x.logs_medio, x.dp / 1e5]
        if geo.es_espiral:
            fila += [round(x.flujo.Dean), x.flujo.correccion_f]
        filas.append(fila)
    s.append(_tabla(encab, filas))

    if eq.tubo.sensibilidad:
        s.append(f"\n**Sensibilidad al factor de la partícula rápida (DN{t.DN}).** El diseño usa el "
                 "de tubo recto, que es el lado seguro. Si el ensayo de trazador midiera un factor "
                 "menor, el tubo podría acortarse así (mismo factor para los tres productos):\n")
        encab = ["v_max/v_media", "Largo (m)", "Volumen (L)", "Δp máx (bar)"]
        filas = [[f"diseño ({t.productos[t.gobierna].flujo.factor:.2f}, {t.gobierna})",
                  t.L, t.volumen_L, t.dp_max / 1e5]]
        filas += [[pt.factor, pt.L, pt.volumen_L, pt.dp_max / 1e5] for pt in eq.tubo.sensibilidad]
        if geo.es_espiral:
            encab += ["Vueltas", "Alto (m)"]
            filas[0] += [t.vueltas, t.alto if t.alto is not None else "—"]
            for fila, pt in zip(filas[1:], eq.tubo.sensibilidad):
                fila += [pt.vueltas, pt.alto if pt.alto is not None else "—"]
        s.append(_tabla(encab, filas))

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
    carpeta.mkdir(parents=True, exist_ok=True)
    (carpeta / "reporte.md").write_text(markdown(r, origen), encoding="utf-8")

    with open(carpeta / "operacion.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["producto", "T_entrada_C", "T_agua_C", "T_salida_C", "Q_calent_kW",
                    "Q_regen_kW", "E_calor_MWh_anio", "K_salida_Pa_sn", "n_salida"])
        for k, vs in r.verificaciones.items():
            for v in vs:
                w.writerow([k, v.T_entrada, round(v.T_agua_necesaria, 3), round(v.T_salida, 3),
                            round(v.Q_calent_kW, 3), round(v.Q_regen_kW, 3),
                            round(v.E_calor_MWh_anio, 3), v.K_salida, v.n_salida])

    with open(carpeta / "tubo.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["DN", "D_int_mm", "L_m", "volumen_L", "dp_max_bar", "Q_CIP_m3_h", "cumple",
                    "producto", "regimen", "Re", "factor", "t_rapida_min", "t_media_min",
                    "logs_min", "logs_medio", "dp_bar"])
        for op in r.equipo.tubo.opciones:
            for k, x in op.productos.items():
                w.writerow([op.DN, op.D * 1000, round(op.L, 2), round(op.volumen_L, 1),
                            round(op.dp_max / 1e5, 4), round(op.Q_CIP_m3_h, 2), op.factible,
                            k, x.flujo.regimen, round(x.flujo.Re), round(x.flujo.factor, 3),
                            round(x.t_rapida_min, 3), round(x.t_media_min, 3),
                            round(x.logs_min, 2), round(x.logs_medio, 2), round(x.dp / 1e5, 4)])
    return carpeta / "reporte.md"
