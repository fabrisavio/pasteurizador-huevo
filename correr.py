"""Punto de entrada.

    python correr.py                      # diseña con config/base.yaml
    python correr.py config/otra.yaml     # diseña con otra configuración
    python correr.py --validar            # reproduce el caso de Gut et al. (2003)
"""

import argparse
import sys
from pathlib import Path

from pasteurizador import config, reporte
from pasteurizador.pasteurizador import correr
from pasteurizador.validacion import validar

RAIZ = Path(__file__).parent


def main() -> int:
    ap = argparse.ArgumentParser(description="Modelo paramétrico del pasteurizador")
    ap.add_argument("config", nargs="?", default=RAIZ / "config" / "base.yaml")
    ap.add_argument("--validar", action="store_true", help="corre la validación contra Gut 2003")
    ap.add_argument("--salida", default=RAIZ / "resultados", help="carpeta de resultados")
    args = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")  # Git Bash en Windows no siempre usa UTF-8

    if args.validar:
        res = validar(RAIZ / "config" / "validacion_gut2003.yaml")
        print("Validación contra Gut et al. (2003), Tabla 7, contracorriente pura\n")
        print(res.tabla())
        return 0 if res.ok else 1

    r = correr(config.cargar(args.config))
    ruta = reporte.escribir(r, args.salida, origen=str(Path(args.config).name))
    eq, t = r.equipo, r.equipo.tubo.elegido
    print(f"Regeneración:  {eq.placas_regen} placas ({eq.A_regen:.2f} m²), gobierna {eq.gobierna_regen}")
    print(f"Calentamiento: {eq.placas_calent} placas ({eq.A_calent:.2f} m²), gobierna {eq.gobierna_calent}")
    print(f"Tubo:          DN{t.DN}, {t.L:.1f} m, Δp {t.dp_max/1e5:.2f} bar, gobierna {t.gobierna}")
    for k in r.verificaciones:
        v = r.en_diseno(k)
        print(f"  {k}: agua a {v.T_agua_necesaria:.1f} °C, sale a {v.T_salida:.1f} °C, "
              f"{v.Q_calent_kW:.1f} kW")
    estado = "cumple" if r.cumple else f"{len(r.alertas)} ALERTAS"
    print(f"\n{estado} · {len(r.avisos)} avisos · reporte completo en {ruta}")
    return 0 if r.cumple else 1


if __name__ == "__main__":
    sys.exit(main())
