import csv
import math
import sys
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, List

# ZetaPhase v1.0 frozen constants
GAMMA_1 = 14.134725141734693790457251983562  # first non-trivial zeta zero imag part
LAMBDA = 10.0 * GAMMA_1
BANDS = [0.0, 1.0 / 3.0, 1.0 / 2.0, 2.0 / 3.0]

# SI Planck units (CODATA exact/standard values used for normalization)
PLANCK_UNITS = {
    "mP": 2.176434e-8,               # kg
    "lP": 1.616255e-35,             # m
    "tP": 5.391247e-44,             # s
    "TP": 1.416784e32,              # K
    "EP": 1.9561e9,                 # J
    "rhoP": 5.155e96,               # kg/m^3
}

# Exact conversion factors to SI where practical
EV_TO_J = 1.602176634e-19
GEV_TO_J = 1.602176634e-10
MEV_TO_J = 1.602176634e-13
GEV_C2_TO_KG = 1.78266192e-27
MEV_C2_TO_KG = 1.78266192e-30
MPC_TO_M = 3.0856775814913673e22
KM_TO_M = 1000.0

# Symbols/constant names that should be treated as exactly 1 in their own natural Planck normalization.
SELF_NORMALIZED_IDENTIFIERS = {
    "mp": "mP",  # proton mass is not self-normalized; keep this list for true Planck units only via symbol match below
}

PLANCK_SYMBOLS = {"mP", "lP", "tP", "TP"}


@dataclass
class RowResult:
    normalized_value: Optional[float]
    log_value: Optional[float]
    phase: Optional[float]
    nearest_band: Optional[str]
    residual: Optional[float]
    note: str = ""


def circular_distance(a: float, b: float) -> float:
    d = abs(a - b) % 1.0
    return min(d, 1.0 - d)


def nearest_band(phase: float) -> Tuple[float, str, float]:
    labels = {
        0.0: "0",
        1.0 / 3.0: "1/3",
        1.0 / 2.0: "1/2",
        2.0 / 3.0: "2/3",
    }
    best = min(BANDS, key=lambda q: circular_distance(phase, q))
    return best, labels[best], circular_distance(phase, best)


def phase_from_normalized(x: float) -> Tuple[float, float, str, float]:
    logv = math.log(x)
    phase = (logv / LAMBDA) % 1.0
    _, label, residual = nearest_band(phase)
    return logv, phase, label, residual


def parse_float(text: str) -> Optional[float]:
    if text is None:
        return None
    s = str(text).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def mass_to_kg(value: float, units: str) -> Optional[float]:
    u = units.strip().lower()
    if u == "kg":
        return value
    if u == "gev":
        return value * GEV_C2_TO_KG
    if u == "mev":
        return value * MEV_C2_TO_KG
    if u == "ev":
        return value * 1.78266192e-36
    return None


def energy_to_j(value: float, units: str) -> Optional[float]:
    u = units.strip().lower()
    if u == "j":
        return value
    if u == "gev":
        return value * GEV_TO_J
    if u == "mev":
        return value * MEV_TO_J
    if u == "ev":
        return value * EV_TO_J
    return None


def length_to_m(value: float, units: str) -> Optional[float]:
    u = units.strip().lower()
    if u == "m":
        return value
    if u == "mpc":
        return value * MPC_TO_M
    return None


def time_to_s(value: float, units: str) -> Optional[float]:
    u = units.strip().lower()
    if u == "s":
        return value
    return None


def temp_to_k(value: float, units: str) -> Optional[float]:
    u = units.strip().lower()
    if u == "k":
        return value
    return None


def density_to_si(value: float, units: str) -> Optional[float]:
    u = units.strip().lower()
    if u in {"kg/m^3", "kg/m3"}:
        return value
    return None


def area_to_m2(value: float, units: str) -> Optional[float]:
    u = units.strip().lower()
    if u in {"m^2", "m2"}:
        return value
    return None


def convert_hubble_to_dimensionless(value: float, units: str) -> Optional[float]:
    # Expected template units: km/s/Mpc
    u = units.strip().lower().replace(" ", "")
    if u in {"km/s/mpc", "km/s/ mpc", "km/s/mpc"}:
        h_si = value * KM_TO_M / MPC_TO_M  # s^-1
        return h_si * PLANCK_UNITS["tP"]
    if u in {"s^-1", "1/s", "hz"}:
        return value * PLANCK_UNITS["tP"]
    return None


def compute_normalized_value(row: Dict[str, str]) -> Tuple[Optional[float], str]:
    raw_value = parse_float(row.get("raw_value", ""))
    units = (row.get("raw_units") or "").strip()
    rule = (row.get("normalization_rule") or "").strip()
    symbol = (row.get("symbol") or "").strip()
    const_name = (row.get("constant_name") or "").strip()
    category = (row.get("category") or "").strip()

    if raw_value is None:
        return None, "missing raw_value"

    # Special-case natural self-normalized Planck units.
    if symbol in PLANCK_SYMBOLS and rule == "direct":
        return 1.0, "self-normalized Planck unit"

    # Special-case Newton's G if left as direct in the template: treat as 1 in Planck-unit normalization.
    if symbol == "G" and rule == "direct":
        return 1.0, "G treated as unity in Planck units"

    if rule == "direct":
        if units.lower() == "dimensionless":
            return raw_value, "dimensionless direct"
        # Safe fallback for already-dimensionless values entered with blank/other label.
        if category in {"coupling", "cosmology", "neutrino", "ratio", "scale", "gravity", "composite", "atomic_constant"}:
            return raw_value, "direct as provided"
        return None, f"unsupported direct dimensional quantity with units '{units}'"

    if rule == "m/mP" or rule == "fpi/mP":
        kg = mass_to_kg(raw_value, units)
        if kg is None:
            return None, f"unsupported mass units '{units}'"
        return kg / PLANCK_UNITS["mP"], "mass normalized by Planck mass"

    if rule == "L/lP":
        meters = length_to_m(raw_value, units)
        if meters is None:
            return None, f"unsupported length units '{units}'"
        return meters / PLANCK_UNITS["lP"], "length normalized by Planck length"

    if rule == "A/lP^2":
        m2 = area_to_m2(raw_value, units)
        if m2 is None:
            return None, f"unsupported area units '{units}'"
        return m2 / (PLANCK_UNITS["lP"] ** 2), "area normalized by Planck length squared"

    if rule == "E/EP":
        joules = energy_to_j(raw_value, units)
        if joules is None:
            return None, f"unsupported energy units '{units}'"
        return joules / PLANCK_UNITS["EP"], "energy normalized by Planck energy"

    if rule == "T/TP":
        kelvin = temp_to_k(raw_value, units)
        if kelvin is None:
            return None, f"unsupported temperature units '{units}'"
        return kelvin / PLANCK_UNITS["TP"], "temperature normalized by Planck temperature"

    if rule == "t/tP":
        seconds = time_to_s(raw_value, units)
        if seconds is None:
            return None, f"unsupported time units '{units}'"
        return seconds / PLANCK_UNITS["tP"], "time normalized by Planck time"

    if rule == "rho/rhoP":
        rho = density_to_si(raw_value, units)
        if rho is None:
            return None, f"unsupported density units '{units}'"
        return rho / PLANCK_UNITS["rhoP"], "density normalized by Planck density"

    if rule == "GF*mP^2":
        # In natural units, G_F has units GeV^-2. Multiply by m_P(GeV)^2.
        u = units.strip().lower().replace(" ", "")
        if u not in {"gev^-2", "gev-2", "1/gev^2"}:
            return None, f"unsupported Fermi constant units '{units}'"
        m_planck_gev = PLANCK_UNITS["mP"] / GEV_C2_TO_KG
        return raw_value * (m_planck_gev ** 2), "Fermi constant normalized by Planck mass squared"

    if rule == "H0*tP":
        dimless = convert_hubble_to_dimensionless(raw_value, units)
        if dimless is None:
            return None, f"unsupported H0 units '{units}'"
        return dimless, "Hubble constant normalized by Planck time"

    if rule == "Lambda*lP^2":
        u = units.strip().lower().replace(" ", "")
        if u not in {"1/m^2", "m^-2", "1/m2", "m-2"}:
            return None, f"unsupported cosmological constant units '{units}'"
        return raw_value * (PLANCK_UNITS["lP"] ** 2), "curvature normalized by Planck length squared"

    return None, f"unknown normalization rule '{rule}'"


def compute_row(row: Dict[str, str]) -> RowResult:
    xnorm, note = compute_normalized_value(row)
    if xnorm is None:
        return RowResult(None, None, None, None, None, note)
    if xnorm <= 0:
        return RowResult(xnorm, None, None, None, None, f"non-positive normalized value ({xnorm})")
    logv, phase, band_label, residual = phase_from_normalized(xnorm)
    return RowResult(xnorm, logv, phase, band_label, residual, note)


def write_output(input_path: str, output_path: str) -> None:
    with open(input_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    if "notes" not in fieldnames:
        fieldnames.append("notes")

    output_rows: List[Dict[str, str]] = []
    included_count = 0
    computed_count = 0
    skipped_count = 0

    for row in rows:
        result = compute_row(row)
        out = dict(row)
        out["normalized_value"] = "" if result.normalized_value is None else f"{result.normalized_value:.16e}"
        out["log_value"] = "" if result.log_value is None else f"{result.log_value:.16e}"
        out["phase"] = "" if result.phase is None else f"{result.phase:.16e}"
        out["nearest_band"] = "" if result.nearest_band is None else result.nearest_band
        out["residual"] = "" if result.residual is None else f"{result.residual:.16e}"
        out["notes"] = result.note
        output_rows.append(out)

        if result.phase is None:
            skipped_count += 1
        else:
            computed_count += 1
            include_flag = (row.get("include_in_primary_stats") or "").strip().lower()
            if include_flag == "yes":
                included_count += 1

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Input: {input_path}")
    print(f"Output: {output_path}")
    print(f"Computed rows: {computed_count}")
    print(f"Skipped rows: {skipped_count}")
    print(f"Primary rows with computed values: {included_count}")
    print(f"Lambda = {LAMBDA:.12f}")


def main() -> int:
    input_path = sys.argv[1] if len(sys.argv) > 1 else "constants_master.csv"
    output_path = sys.argv[2] if len(sys.argv) > 2 else "constants_phase_results.csv"
    try:
        write_output(input_path, output_path)
        return 0
    except FileNotFoundError as e:
        print(f"File not found: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
