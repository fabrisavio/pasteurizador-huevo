"""Punto de entrada.

    python correr.py                         # diseña con config/base.yaml
    python correr.py config/mi_caso.yaml     # diseña con otra config (copia de base.yaml)
    python correr.py --validar               # reproduce el caso de Gut et al. (2003)
"""

import argparse
import sys
from pathlib import Path

from pasteurizador import __version__, config, reporte
from pasteurizador.pasteurizador import correr
from pasteurizador.validacion import chequear_U_gut, tabla_chequeo_U, validar

RAIZ = Path(__file__).parent  # la carpeta donde está este archivo


def _sin_config(ruta: Path) -> int:
    """Mensaje claro cuando la config no existe, en vez de un traceback."""
    disponibles = sorted(p.name for p in (RAIZ / "config").glob("*.yaml"))
    print(f"No encuentro la configuración '{ruta}'.\n")
    print("Configuraciones disponibles en config/: " + ", ".join(disponibles))
    print("\nPara probar otro caso, copiá la base y editá la copia:")
    print("    cp config/base.yaml config/mi_caso.yaml")
    print("    python correr.py config/mi_caso.yaml")
    return 2  # código de salida distinto de 0 = algo salió mal


def main() -> int:
    ap = argparse.ArgumentParser(description="Modelo paramétrico del pasteurizador")
    ap.add_argument("config", nargs="?", default=RAIZ / "config" / "base.yaml",
                    help="archivo YAML de diseño (por defecto config/base.yaml)")
    ap.add_argument("--validar", action="store_true", help="corre la validación contra Gut 2003")
    ap.add_argument("--salida", default=RAIZ / "resultados", help="carpeta de resultados")
    args = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")  # Git Bash en Windows no siempre usa UTF-8

    if args.validar:
        res = validar(RAIZ / "config" / "validacion_gut2003.yaml")
        print("Validación contra Gut et al. (2003), Tabla 7, contracorriente pura\n")
        print("  T1 yema cruda (5 °C)")
        print("   ▼")
        print("  [regeneración, lado frío] ─T2─► [calentamiento] ─T3─► tubo ─T4─┐")
        print("  T6 ◄─ [enfriamiento] ◄─T5─ [regeneración, lado caliente] ◄─────┘\n")
        print(res.tabla())
        print("\nChequeo del U (v1): con la correlación y los datos que publica el paper")
        print("(placas planas, Sharifi et al. 1995; ensuciamientos de Gut; T del modelo de U constante)\n")
        print(tabla_chequeo_U(chequear_U_gut(RAIZ / "config" / "validacion_gut2003.yaml")))
        print("\nEl U de la Tabla 7 no sale de los datos publicados: por eso v1 no lo usa.")
        return 0 if res.ok else 1

    ruta_cfg = Path(args.config)
    if not ruta_cfg.exists():
        return _sin_config(ruta_cfg)
    try:
        cfg = config.cargar(ruta_cfg)
    except config.ConfigInvalida as e:
        print(e)
        return 2

    r = correr(cfg)
    ruta = reporte.escribir(r, args.salida, origen=ruta_cfg.name)
    eq, t = r.equipo, r.equipo.tubo.elegido

    print(f"Pasteurizador v{__version__} · config {ruta_cfg.name}\n")
    print("EQUIPO")
    for nombre, n, A, gob, arr, attr in (
            ("Regeneración: ", eq.placas_regen, eq.A_regen, eq.gobierna_regen, eq.arreglo_regen, "regen"),
            ("Calentamiento:", eq.placas_calent, eq.A_calent, eq.gobierna_calent, eq.arreglo_calent, "calent")):
        linea = f"  {nombre} {n} placas ({A:.2f} m²), gobierna {gob}"
        if arr is not None:
            # U de cada producto en el equipo, a su T de entrada de diseño
            Us = [getattr(r.en_diseno(k), attr).U for k in r.verificaciones
                  if getattr(r.en_diseno(k), attr) is not None]
            linea += f"; {arr.pasos} pasos × {arr.canales} canales, U {min(Us):.0f}–{max(Us):.0f} W/m²K"
        print(linea)
    print(f"  Agua caliente: {eq.m_agua*3.6:.1f} m³/h")
    print(f"  Tubo:          {reporte.describir_tubo(r)}")
    print(f"                 {t.volumen_L:.0f} L, Δp {t.dp_max/1e5:.2f} bar "
          f"(máx {cfg.dp_max_holding/1e5:.1f}), CIP {t.Q_CIP_m3_h:.1f} m³/h, largo fijado por {t.gobierna}")
    if t.horquillas:
        h = t.horquillas
        print(f"                 envolvente {h.largo:.2f} × {h.ancho:.2f} m en planta, "
              f"sube {h.desnivel:.2f} m de la entrada a la salida")
    if eq.aislacion:
        a = cfg.aislacion
        elegida = eq.aislacion[eq.aislacion_elegida]
        peor = max(eq.perdida_tubo, key=eq.perdida_tubo.get)
        print(f"  Aislación:     {reporte.NOMBRE_AISLACION[elegida.forma].lower()}, "
              f"{a.espesor*1000:.0f} mm: pierde hasta {eq.perdida_tubo[peor]:.2f} °C ({peor}), "
              f"{elegida.E_MWh_anio:.2f} MWh/año")
        for otra in (o for f, o in eq.aislacion.items() if f != elegida.forma):
            print(f"                 {reporte.NOMBRE_AISLACION[otra.forma].lower()}: "
                  f"{otra.dT[peor]:.2f} °C, {otra.E_MWh_anio:.2f} MWh/año, "
                  f"{eq.placas_calent_por_aislacion[otra.forma]} placas de calentamiento "
                  f"contra {eq.placas_calent_por_aislacion[elegida.forma]}")
    else:
        print(f"  Pérdida en el tubo: {cfg.perdida_fija} °C (fija)")
    print(f"  Control:       {cfg.margen_control:.1f} °C sobre la T legal")
    print(f"  U:             {'por correlación de ' + cfg.correlacion.capitalize() if cfg.con_correlacion else 'fijo (v0)'}")

    # Tabla de ancho fijo: {x:>7} alinea a la derecha en 7 caracteres
    print("\nPRODUCTO (a su T de entrada de diseño; temperaturas en °C)")
    print(f"  {'':4}{'kg/h':>6} {'entra':>6} {'regen':>6} {'tubo':>11} {'secador':>8} "
          f"{'agua':>6} {'Q cal':>7} {'t rápida':>9}")
    for k in r.verificaciones:
        p, v, x = cfg.productos[k], r.en_diseno(k), t.productos[k]
        tubo = f"{v.T_calent_salida:.1f}→{cfg.T_operacion(p):.1f}"
        print(f"  {k:4}{p.m*3600:>6.0f} {v.T_entrada:>6.1f} {v.T_regen_salida:>6.1f} {tubo:>11} "
              f"{v.T_salida:>8.1f} {v.T_agua_necesaria:>6.1f} {v.Q_calent_kW:>5.1f} kW "
              f"{x.t_rapida_min:>5.2f} min")
    print("  regen: sale de la regeneración · tubo: entra→sale (T legal + margen) · agua: T de entrada del agua caliente")

    estado = "cumple" if r.cumple else f"{len(r.alertas)} ALERTAS"
    print(f"\n{estado} · {len(r.avisos)} avisos · reporte completo en {ruta}")
    return 0 if r.cumple else 1


if __name__ == "__main__":
    # Solo corre si se ejecuta este archivo directamente (no si se importa)
    sys.exit(main())
