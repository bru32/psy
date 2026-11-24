# Dealing with SI and IP
# In addition, allowing user interaction in preferred unit (e.g g/kg and %)

from dataclasses import dataclass, fields
from typing import ClassVar, Callable

# Helper factories for linear scale vs. temperature offset conversions
def scale(factor: float) -> tuple[Callable[[float], float], Callable[[float], float]]:
  # (SI -> Display, Display -> SI)
  return (lambda si: si * factor, lambda disp: disp / factor)

C_TO_F = (lambda c: c * 9.0 / 5.0 + 32.0, lambda f: (f - 32.0) * 5.0 / 9.0)
IDENTITY = (lambda x: x, lambda x: x)


@dataclass(slots=True)
class TAir:
  # Always stored internally in strict SI base units
  pa: float = 101.325  # kPa
  db: float = 23.0     # °C
  wb: float = 16.0     # °C
  ha: float = 45.0     # kJ/kg
  wa: float = 8.6e-3   # kg/kg
  rh: float = 0.5      # fraction
  td: float = 13.0     # °C
  pv: float = 1.2      # kPa
  ps: float = 2.4      # kPa
  va: float = 0.86     # m³/kg
  ro: float = 1.19     # kg/m³
  cp: float = 1.02     # kJ/kg-K
  za: float = 44.0     # kJ/kg
  units: str = "SI"    # "SI" or "IP" (controls input/output boundary)

  # Map: system -> field -> (unit_label, (si_to_display, display_to_si))
  _SYSTEMS: ClassVar[dict[str, dict[str, tuple[str, tuple[Callable, Callable]]]]] = {
      "SI": {
          "pa": ("kPa",     IDENTITY),
          "db": ("°C",      IDENTITY),
          "wb": ("°C",      IDENTITY),
          "ha": ("kJ/kg",   IDENTITY),
          "wa": ("g/kg",    scale(1000.0)),   # User enters/reads g/kg -> stored as kg/kg
          "rh": ("%",       scale(100.0)),    # User enters/reads %    -> stored as fraction
          "td": ("°C",      IDENTITY),
          "pv": ("kPa",     IDENTITY),
          "ps": ("kPa",     IDENTITY),
          "va": ("m³/kg",   IDENTITY),
          "ro": ("kg/m³",   IDENTITY),
          "cp": ("kJ/kg-K", IDENTITY),
          "za": ("kJ/kg",   IDENTITY),
      },
      "IP": {
          "pa": ("psi",       scale(0.1450377)),
          "db": ("°F",        C_TO_F),
          "wb": ("°F",        C_TO_F),
          "ha": ("Btu/lb",    scale(0.4299226)),
          "wa": ("gr/lb",     scale(7000.0)),   # Grains of moisture per lb dry air -> kg/kg
          "rh": ("%",         scale(100.0)),
          "td": ("°F",        C_TO_F),
          "pv": ("psi",       scale(0.1450377)),
          "ps": ("psi",       scale(0.1450377)),
          "va": ("ft³/lb",    scale(16.01846)),
          "ro": ("lb/ft³",    scale(0.06242796)),
          "cp": ("Btu/lb-°F", scale(0.2388459)),
          "za": ("Btu/lb",    scale(0.4299226)),
      },
  }

  def __post_init__(self):
    self.units = self.units.upper()
    if self.units not in self._SYSTEMS:
      raise ValueError(f"Unsupported unit system: '{self.units}'. Use 'SI' or 'IP'.")

    # Convert user-supplied values from chosen unit system into internal SI
    # Note: If instantiated with no args, defaults above are SI—so only convert passed args,
    # OR store defaults in SI and skip conversion unless passed (handled cleanly via factory below).
    conv_map = self._SYSTEMS[self.units]
    for k, (_, (_, to_si)) in conv_map.items():
      # Convert defaults if IP was selected, or better: define SI defaults and convert them to IP first
      pass

  @classmethod
  def from_inputs(cls, units: str = "SI", **kwargs) -> "TAir":
    """Recommended constructor: starts from exact SI defaults, then applies user inputs in 'units'."""
    units = units.upper()
    if units not in cls._SYSTEMS:
      raise ValueError(f"Unsupported unit system: '{units}'. Use 'SI' or 'IP'.")
    obj = cls(units=units)
    obj.set(**kwargs)
    return obj

  def set(self, **kwargs) -> None:
    conv_map = self._SYSTEMS[self.units]
    for k, v in kwargs.items():
      if k not in conv_map:
        raise AttributeError(f"Unknown air property: '{k}'")
      _, (_, to_si) = conv_map[k]
      setattr(self, k, to_si(v))

  def get_display(self, name: str) -> float:
    """Returns internal SI value converted to the active display unit system."""
    _, (to_display, _) = self._SYSTEMS[self.units][name]
    return to_display(getattr(self, name))

  def report(self, title: str = "Moist Air State") -> str:
    conv_map = self._SYSTEMS[self.units]
    lines = [f"{title} ({self.units})"]
    for key, (unit_label, (to_display, _)) in conv_map.items():
      val = to_display(getattr(self, key))
      lines.append(f"  {key:<3} = {val:0.6g} {unit_label}")
    return "\n".join(lines)


# --- Usage ---

# 1. SI User (enters db in °C, wa in g/kg, rh in %)
air_si = TAir.from_inputs(units="SI", db=32.0, wa=9.123, rh=45.5)
print(air_si.report())
print("Internal SI:", air_si.db, air_si.wa, air_si.rh)  # 32.0, 0.009123, 0.455

# 2. US / IP User (enters db in °F, pa in psi, wa in grains/lb, rh in %)
air_ip = TAir.from_inputs(units="IP", pa=14.696, db=89.6, wa=63.86, rh=45.5)
print(air_ip.report())
print("Internal SI:", round(air_ip.db, 2), round(air_ip.wa, 6), air_ip.rh)  # 32.0, 0.009123, 0.455

# 3. Seamless switching: compute in SI, then switch display units anytime!
air_si.units = "IP"
print(air_si.report("Same Air State Viewed in IP"))
