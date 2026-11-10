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
from . import aislacion as ai
from . import dimensiones as dm
from . import placas as pl
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
    if geo.es_horquillas and t.horquillas:
        h = t.horquillas
        return (f"DN{t.DN} (Ø int {t.D*1000:.0f} mm), {t.L:.1f} m en {h.tramos} tramos de "
                f"{h.largo_tramo:.2f} m y {h.curvas} curvas de 180°, en dos planos; "
                f"sube {h.desnivel:.2f} m ({h.pendiente:.1%})")
    if not geo.es_espiral:
        return f"DN{t.DN} (Ø int {t.D*1000:.0f} mm), {t.L:.1f} m en tramos rectos"
    forma = (f"DN{t.DN} (Ø int {t.D*1000:.0f} mm), {t.L:.1f} m enrollados en espiral de "
             f"Ø {geo.D_espiral:.2f} m: {t.vueltas:.1f} vueltas")
    if t.alto is not None:
        forma += f", {t.alto:.2f} m de alto"
    return forma


NOMBRE_AISLACION = {"tubo": "Coquilla a lo largo del tubo",
                    "tubo_curvas_desnudas": "Coquilla en los tramos, curvas al aire",
                    "bloque": "Camisa en bloque"}


def aislacion_md(r: Resultado) -> str:
    """Comparación de las dos formas de aislar el tubo y lo que ahorra la elegida."""
    cfg, eq = r.cfg, r.equipo
    a, ops = cfg.aislacion, eq.aislacion
    elegida = eq.aislacion_elegida or a.forma
    s = [f"\n### Aislación del tubo\n",
         f"Lana mineral de {a.espesor*1000:.0f} mm (k = {a.k} W/m·K), sala a {a.T_amb:.0f} °C. "
         f"Diseño con: **{NOMBRE_AISLACION[elegida].lower()}**.\n"]
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
    if "tubo" in ops and "tubo_curvas_desnudas" in ops:
        t, d = ops["tubo"], ops["tubo_curvas_desnudas"]
        tubo = eq.tubo.elegido
        h, D_ext = tubo.horquillas, tubo.D_ext or tubo.D
        L_d = h.curvas * (h.largo_curva + 2 * a.union)
        UA_pelado = ai.ua_desnudo(L_d, D_ext, a.h_desnudo)[0]
        por_metro = ai.ua_desnudo(1.0, D_ext, a.h_desnudo)[0] / ai.ua_tubo(1.0, D_ext, a)[0]
        peor = max(claves, key=lambda k: d.dT[k] - t.dT[k])
        s.append(f"\n**Lo que cuesta dejar las curvas al aire.** Quedan {L_d:.1f} m de caño pelado "
                 f"({h.curvas} curvas de {h.largo_curva*1000:.0f} mm más {a.union*1000:.0f} mm de unión "
                 f"a cada lado): {L_d / tubo.L:.0%} del largo, pero {UA_pelado / d.UA:.0%} de la "
                 f"pérdida, porque un caño pelado pierde {por_metro:.0f} veces más por metro que uno "
                 f"con coquilla. El producto tiene que entrar {d.dT[peor] - t.dT[peor]:.2f} °C más "
                 f"caliente al tubo ({peor}); {(d.E_MWh_anio - t.E_MWh_anio)*1000:.0f} kWh/año más; "
                 f"placas de calentamiento {placas.get('tubo', '—')} con todo aislado y "
                 f"{placas.get('tubo_curvas_desnudas', '—')} con las curvas al aire.\n")
    return "\n".join(s)


def dimensiones_md(r: Resultado) -> str:
    """Las cotas para el plano: paquete de placas, bastidor de referencia y tubo."""
    cfg, eq = r.cfg, r.equipo
    bas = cfg.bastidor
    s = ["\n## Dimensiones para el plano\n"]
    if cfg.placa is not None:
        pq = dm.paquete(eq.placas_regen, eq.placas_calent, cfg.placa.b, cfg.placa.espesor, bas)
        lo, hi = bas.largo_rango
        s.append(f"**Intercambiador de placas.** {pq.placas_termicas} placas térmicas más "
                 f"{pq.placas_total - pq.placas_termicas} de extremo, a {pq.paso*1000:.1f} mm cada una "
                 f"(canal + chapa), más {bas.conector*1000:.0f} mm de bastidor intermedio entre "
                 f"secciones: paquete de **{pq.largo:.2f} m** apretado. Bastidor de referencia con la "
                 f"misma placa (canal de 3 mm, conexiones DN50): {bas.alto:.2f} m de alto × "
                 f"{bas.ancho:.2f} m de ancho, {lo:.2f} a {hi:.2f} m de largo según el modelo. "
                 f"Con lugar para correr la placa de presión y sacar placas (otro paquete de largo), "
                 f"hace falta un bastidor de **~{pq.largo_bastidor:.1f} m**. La cota final la da "
                 f"el proveedor.\n")
    t, geo = eq.tubo.elegido, eq.tubo.geometria
    if geo.es_horquillas and t.horquillas:
        h = t.horquillas
        s.append(f"**Tubo de retención.** {h.tramos} tramos de {h.largo_tramo:.2f} m que alternan "
                 f"entre dos planos verticales separados {h.entre_planos*1000:.0f} mm (entre ejes), "
                 f"unidos por {h.curvas} curvas de 180° de {h.largo_curva*1000:.0f} mm de desarrollo "
                 + (f"(dos codos de 90° y un niple de {h.niple*1000:.0f} mm). " if h.niple > 1e-6 else ". ")
                 + f"Sube {h.pendiente:.1%} continuo: **{h.desnivel:.2f} m** de la entrada a la salida, "
                 f"{2 * h.pendiente * (h.largo_tramo + h.largo_curva)*1000:.0f} mm entre dos tramos "
                 f"del mismo plano. Envolvente con la aislación (Ø {h.D_envolvente*1000:.0f} mm): "
                 f"**{h.largo:.2f} m × {h.ancho:.2f} m × {h.alto:.2f} m** de alto, más la altura a la "
                 "que entre el primer tramo y la estructura de soporte (fija, CFIA).\n")
    elif geo.es_espiral and t.alto:
        s.append(f"**Tubo de retención.** Espiral de Ø {geo.D_espiral:.2f} m al eje, "
                 f"{t.alto:.2f} m de alto.\n")
    return "\n".join(s)


def placas_md(r: Resultado) -> str:
    """v1: arreglos evaluados, h y Δp de cada lado, contraste y CIP."""
    cfg, eq = r.cfg, r.equipo
    claves = list(cfg.productos)
    pl_ = cfg.placa
    s = ["## Placas: arreglo, coeficientes y Δp\n",
         f"Placa chevron de {pl_.angulo:.0f}°, b = {pl_.b*1000:.1f} mm, ancho {pl_.ancho:.2f} m, "
         f"φ = {pl_.phi}, Dh = {pl_.Dh*1000:.2f} mm, largo efectivo {pl_.L:.2f} m. "
         f"h por {cfg.correlacion.capitalize()}; Δp máximo {cfg.dp_max_placas/1e5:.1f} bar por "
         "sección y corriente. Ensuciamiento: producto "
         f"{cfg.Rf_regen:.1e} m²K/W por lado en regeneración y {cfg.Rf_calent:.1e} en "
         f"calentamiento; agua {cfg.Rf_agua:.1e}.\n"]
    for nombre, ops, elegido in (("Regeneración", eq.opciones_regen, eq.arreglo_regen),
                                 ("Calentamiento", eq.opciones_calent, eq.arreglo_calent)):
        s.append(f"\n**{nombre}: arreglos evaluados.**\n")
        corrientes = list(ops[0].dp)
        filas = []
        for o in ops:
            marca = " ←" if o.arreglo == elegido else ""
            filas.append([f"{o.arreglo}{marca}", o.placas, o.gobierna]
                         + [f"{o.dp[c]/1e5:.2f}" for c in corrientes]
                         + ["sí" if o.factible else "no"])
        s.append(_tabla(["Pasos × canales", "Placas", "Gobierna"]
                        + [f"Δp {c} (bar)" for c in corrientes] + ["Cumple Δp"], filas))

    s.append("\n**Cada lado, con el equipo elegido** (T de entrada de diseño):\n")
    filas = []
    for k in claves:
        v = r.en_diseno(k)
        for sec, det in (("regeneración", v.det_regen), ("calentamiento", v.det_calent)):
            if det is None:
                continue
            for lado, c in (("caliente", det.caliente), ("fría", det.fria)):
                filas.append([k, sec, f"{lado} ({c.fluido})", c.v, c.Re, c.Pr, c.h, c.dp / 1e5,
                              round(c.gamma_pared)])
    s.append(_tabla(["Producto", "Sección", "Lado", "v (m/s)", "Re", "Pr", "h (W/m²K)",
                     "Δp (bar)", "γ̇ pared (1/s)"], filas))

    s.append("\n**De dónde sale cada U**: qué fracción de la resistencia pone cada parte.\n")
    filas = []
    for k in claves:
        v = r.en_diseno(k)
        for sec, det in (("regeneración", v.det_regen), ("calentamiento", v.det_calent)):
            if det is None:
                continue
            f = det.fraccion_resistencias
            contraste = v.U_contraste.get(sec)
            filas.append([k, sec, round(det.U), round(det.U_limpio),
                          f"{f['película caliente']:.0%}", f"{f['película fría']:.0%}",
                          f"{f['ensuciamiento']:.0%}", f"{f['pared']:.0%}",
                          round(contraste) if contraste else "—"])
    corr = (cfg.correlacion_contraste or "").capitalize() or "—"
    s.append(_tabla(["Producto", "Sección", "U (W/m²K)", "U limpio (W/m²K)", "Película caliente",
                     "Película fría", "Ensuciamiento", "Pared", f"U con {corr}"], filas))
    if cfg.correlacion_contraste:
        lo, hi = pl.RANGO_RE[cfg.correlacion_contraste]
        s.append(f"\n{corr} está ajustada para Re entre {lo:g} y {hi:g}: con los productos "
                 "(Re bajos) el contraste es orientativo.\n")

    if eq.Q_CIP_placas:
        s.append(f"\n**CIP.** Para Re ≥ {cfg.Re_CIP_placas:.0f} en los canales de producto "
                 f"(Galiński et al., 2025), con solución a {cfg.T_CIP:.0f} °C: "
                 + ", ".join(f"{k} {q:.1f} m³/h" for k, q in eq.Q_CIP_placas.items())
                 + f". El tubo pide {eq.tubo.elegido.Q_CIP_m3_h:.1f} m³/h (a {cfg.v_CIP} m/s); "
                 "el circuito de producto está en serie, así que la bomba tiene que dar el mayor.\n")
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
    for nombre, n, A, U, gob, attr, arr in [
        ("Regeneración", eq.placas_regen, eq.A_regen, cfg.U_regen, eq.gobierna_regen, "regen",
         eq.arreglo_regen),
        ("Calentamiento", eq.placas_calent, eq.A_calent, cfg.U_calent, eq.gobierna_calent, "calent",
         eq.arreglo_calent),
    ]:
        pedidas = []
        for k in claves:
            sec = getattr(eq.disenos[k], attr)  # getattr(obj, "regen") = obj.regen
            pedidas.append(f"{k} {sec.A:.2f}" if sec else f"{k} —")
        if cfg.con_correlacion:
            # U de cada producto operando en el equipo, a su T de entrada de diseño
            U = " · ".join(f"{k} {getattr(r.en_diseno(k), attr).U:.0f}" for k in claves
                           if getattr(r.en_diseno(k), attr) is not None)
            filas.append([nombre, str(arr), n, A, U, gob, " · ".join(pedidas)])
        else:
            filas.append([nombre, n, A, U, gob, " · ".join(pedidas)])
    if cfg.con_correlacion:
        s.append(_tabla(["Sección", "Pasos × canales", "Placas térmicas", "Área instalada (m²)",
                         "U en operación (W/m²K)", "Gobierna", "Área que pide cada producto (m²)"],
                        filas))
    else:
        s.append(_tabla(["Sección", "Placas térmicas", "Área instalada (m²)", "U (W/m²K)",
                         "Gobierna", "Área que pide cada producto (m²)"], filas))
    s.append(f"\nPlaca de {cfg.area_placa} m². Bomba de agua caliente: "
             f"{eq.m_agua:.2f} kg/s ({eq.m_agua*3.6:.1f} m³/h aprox.).\n")
    if cfg.con_correlacion:
        s.append(placas_md(r))

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
    if geo.es_horquillas:
        s.append(f"Horquillas: partícula más rápida como en tubo recto, v_max/v_media = (3n+1)/(n+1) "
                 f"en laminar. Subida continua mínima {geo.pendiente:.0%}, curvas incluidas (CFIA). "
                 f"Curvas de 180° de R = {geo.radio_curva*1000:.0f} mm; Δp por el método 2-K "
                 f"(Hooper, 1981) con K1 = {geo.K1_curva:.0f} y K∞ = {geo.Kinf_curva}.\n")
    encab = ["DN", "D int (mm)", "Largo (m)", "Volumen (L)", "Δp máx (bar)", "Q CIP (m³/h)"]
    if geo.es_espiral:
        encab += ["Vueltas", "Alto (m)"]
    if geo.es_horquillas:
        encab += ["Tramos", "Desnivel (m)"]
    encab.append("Cumple Δp")
    filas = []
    for op in eq.tubo.opciones:
        # Δp con 3 decimales: si no, un 1,004 bar se ve "1.00" y parece que cumple
        fila = [op.DN, op.D * 1000, op.L, op.volumen_L, f"{op.dp_max / 1e5:.3f}", op.Q_CIP_m3_h]
        if geo.es_espiral:
            fila += [op.vueltas, op.alto if op.alto is not None else "—"]
        if geo.es_horquillas and op.horquillas:
            fila += [op.horquillas.tramos, op.horquillas.desnivel]
        fila.append("sí" if op.factible else "no")
        filas.append(fila)
    s.append(_tabla(encab, filas))

    s.append(f"\nCada producto en el DN{t.DN}:\n")
    encab = ["Producto", "Régimen", "Re (M–R)"]
    if geo.es_espiral:
        encab += ["Dean", "f espiral / f recto"]
    encab += ["v_max/v_media", "t partícula rápida (min)", "t medio (min)",
              "log rápida", "log media", "Δp (bar)"]
    if geo.es_horquillas:
        encab.append("Δp de las curvas (bar)")
    filas = []
    for k, x in t.productos.items():
        fila = [k, x.flujo.regimen.replace("transicion", "transición"), round(x.flujo.Re)]
        if geo.es_espiral:
            fila += [round(x.flujo.Dean), x.flujo.correccion_f]
        fila += [x.flujo.factor, x.t_rapida_min, x.t_media_min, x.logs_min, x.logs_medio,
                 x.dp / 1e5]
        if geo.es_horquillas:
            fila.append(f"{x.dp_curvas / 1e5:.4f}")
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

    s.append(dimensiones_md(r))

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
        w.writerow(["DN", "D_int_mm", "L_m", "volumen_L", "vueltas", "alto_m", "tramos",
                    "desnivel_m", "dp_max_bar",
                    "Q_CIP_m3_h", "cumple", "producto", "regimen", "Re", "Dean", "factor",
                    "t_rapida_min", "t_media_min", "logs_min", "logs_medio", "dp_bar"])
        for op in r.equipo.tubo.opciones:
            h = op.horquillas
            for k, x in op.productos.items():
                w.writerow([op.DN, op.D * 1000, round(op.L, 2), round(op.volumen_L, 1),
                            round(op.vueltas, 2) if op.vueltas else "",
                            round(op.alto, 3) if op.alto else "",
                            h.tramos if h else "", round(h.desnivel, 3) if h else "",
                            round(op.dp_max / 1e5, 4), round(op.Q_CIP_m3_h, 2), op.factible,
                            k, x.flujo.regimen, round(x.flujo.Re), round(x.flujo.Dean, 1),
                            round(x.flujo.factor, 3), round(x.t_rapida_min, 3),
                            round(x.t_media_min, 3), round(x.logs_min, 2),
                            round(x.logs_medio, 2), round(x.dp / 1e5, 4)])
    return carpeta / "reporte.md"
