# psy_wernick.py
#
# Consistent and accurate set of Psychrometric Equations.
#
# Uses: 
#   ASHRAE sat pressure correlation coefficients.
#   Constants and enhancement factor from Gatley.
#
# Interface units: 
#   Temperature     [°C]
#   Pressure        [kPa]
#   Abs Humidity    [kg_w / kg_da]
#   Rel Humidity    [fraction]
#   Spec Enthalpy   [kJ/kg]
#   Spec Heat       [kJ/kg-K]
#   Specific Vol    [m³/kg_a]
#   Density         [kg_av / m³]
#   Entropy         [kJ/kg]
#   Sigma Heat      [kJ/kg]
#   Degree of Sat   [fraction]
#
# Critical Equation:
# Efficient and safe solution of the adiabatic 
# saturation relationship.
#
# Special Features:
# Includes general air state solver 
#   air = calcstate(pa, ha=45, wa=8.5e-3)
#   air = calcstate(pa, db=23, td=12)
# calcstate can take ANY combination of properties.
# It uses the following as an error function:
#   air = astate(pa, db, wb)
# and this means that ALL combinations of calculation 
# give consistent results.
#
# Psychrometric Chart Plotter:
#   Includes interactive mouse
#   Linear in energy and mass:
#     Therefore mix falls exactly on line between end points.
#
# Author: Bruce Wernick
# 24 November 2026
# Version: 0.0.0.1
#
# Copyright © (2026) Bruce Wernick,  All rights reserved...
# ---------------------------------------------------------

import sys
import math
from nu_meth import brent, broyden_nd as fsolve
from dataclasses import dataclass, astuple, asdict
from typing import ClassVar, NamedTuple

# =========================================================
# ASHRAE Fundamentals Constants 
# (Internal SI: [K], [Pa]; Interface: [°C], [kPa], [g/kg])
# =========================================================

# Thermodynamic constants (ASHRAE Fundamentals Chapter 1)
P0 = 101.325          # Barom pressure sea level [kPa]
K0 = 273.15           # °C to Kelvin offset
PS0 = 0.611213        # Saturation pressure at 0°C [kPa]
fs = 1.0044           # molecular interaction coeff
TMIN = -100.0         # Minimum valid temperature [°C]
TMAX = 200.0          # Maximum valid temperature [°C]
TRIPLE = 0.01         # Triple point temperature [°C]
MIN_WA = 1e-7         # Minimum humidity ratio [kg/kg]
LN_1000 = math.log(1000.0)

MA = 28.9645077       # molecular weight of air
MW = 18.01528         # molecular weight of water
MWA = MW / MA         # Ratio of molar masses (water / dry air)
RU = 8.31446261815324 # universal gas constant
RA = RU / MA          # air gas constant
RW = RU / MW          # water gas constant
RAW = RA / RW         # gas constant ratio

# Specific heats and latent heats in [kJ/kg] and [kJ/kg-K]
CPda = 1.006          # Dry air specific heat
CPv = 1.82            # Water vapor specific heat
CPw = 4.186           # Liquid water specific heat
CPi = 2.093           # Ice specific heat
Hfgo = 2500.77        # Latent heat of vaporization at 0°C [kJ/kg]
Hsgo = 2834.2         # Latent heat of sublimation at 0°C [kJ/kg]

# =========================================================
# Saturation Pressure and Temperature
# =========================================================

class PhaseCoeffs(NamedTuple):
  inv_t: float # C1 or C8 (1/T term)
  poly: tuple[float, ...] # C2..C6 or C9..C12 (powers of T)
  log_t: float # C7 or C13 (ln(T) term)

# 173.15 K to 273.15 K (-100°C to 0°C)
ICE_COEFFS = PhaseCoeffs(
  inv_t=-5.6745359e3,
  poly=(6.3925247, -9.6778430e-3, 6.2215701e-7, 2.0747825e-9, -9.4840240e-13),
  log_t=4.1635019)

# 273.15 K to 473.15 K (0°C to 200°C)
LIQUID_COEFFS = PhaseCoeffs(
  inv_t=-5.8002206e3,
  poly=(1.3914993, -4.8640239e-2, 4.1764768e-5, -1.4452093e-8),
  log_t=6.5459673)

def horner(coeffs: tuple[float, ...], x: float) -> float:
  """ Evaluate polynomial with ascending-order 
      coefficients using Horner's method.
  """
  result = 0.0
  for c in reversed(coeffs):
    result = result * x + c
  return result

def ln_psat(t: float) -> float:
  """ Natural log of saturation pressure [kPa] 
      for t in [°C].
  """
  tk = t + K0
  coeffs = ICE_COEFFS if t < 0.0 else LIQUID_COEFFS
  ln_pa = (coeffs.inv_t / tk
    + horner(coeffs.poly, tk)
    + coeffs.log_t * math.log(tk))
  return ln_pa - LN_1000

def psat(t):
  """ Water vapor saturation pressure [kPa] 
      for t in [°C].
  """
  if not (TMIN <= t <= TMAX):
    raise ValueError(f"t out of range [{TMIN}, {TMAX}]!")
  return math.exp(ln_psat(t))

def tsat(p):
  """ Invert psat(t) = p safely within [TMIN, TMAX] 
      returns t in °C.
  """
  if p <= psat(TMIN):
    return TMIN

  if p >= psat(TMAX):
    return TMAX

  lnp = math.log(p)
  t, its, ok = brent(lambda t: ln_psat(t) - lnp, TMIN, TMAX)
  return t

def wsat(pa, t):
  """ Saturation humidity ratio [kg_w/kg_da] 
      at temperature t [°C] and pressure p [kPa].
  """
  ps = fs * psat(t)
  if ps >= pa:
    raise ValueError(f"ps >= pa!")
  return MWA * ps / (pa - ps)

# =========================================================
# Adiabatic Saturation: Wet-Bulb Solver
# =========================================================

def wa_from_dbwb(pa, db, wb, over_ice=None):
  """ Calculate humidity ratio W [kg_w/kg_da] 
      from pressure in [kPa], dry-bulb and 
      wet-bulb temperatures in [°C] 
      (ASHRAE Ch. 1).
  """
  if over_ice is None:
    over_ice = wb < 0.0

  wsb = wsat(pa, wb)
  if over_ice:
    # Sublimation energy balance (wb < 0°C)
    num = (Hsgo - (CPi - CPv) * wb) * wsb - CPda * (db - wb)
    den = Hsgo + CPv * db - CPi * wb
  else:
    # Vaporization energy balance (wb >= 0°C)
    num = (Hfgo - (CPw - CPv) * wb) * wsb - CPda * (db - wb)
    den = Hfgo + CPv * db - CPw * wb
  return num / den

def wb_from_dbwa(pa, db, wa):
  """ Calculate thermodynamic wet-bulb temperature [°C]
      given barometric pressure [kPa], 
      dry-bulb temperature [°C] and 
      humidity ratio [kg_w/kg_da].
  """
  if wa < 0.0:
    raise ValueError(f"wa negative!")
  
  if not (TMIN <= db <= TMAX):
    raise ValueError(f"db out of range [{TMIN}, {TMAX}]!")

  # 1. Establish safe upper temperature 
  # bound strictly below boiling point at pressure P
  t_boil = tsat(pa)
  t_upper_limit = min(db, t_boil - 1e-4)

  if t_upper_limit <= TMIN:
    return TMIN

  # 2. Check upper bound using the exact root function
  use_ice_upper = t_upper_limit <= TRIPLE
  wmax = wa_from_dbwb(pa, db, t_upper_limit, over_ice=use_ice_upper)
  if wa >= wmax:
    return t_upper_limit

  # 3. Check lower physical bound at TMIN
  wmin = wa_from_dbwb(pa, db, TMIN, over_ice=True)
  if wa <= wmin:
    return TMIN

  # 4. Split search domain at triple point
  if t_upper_limit <= TRIPLE:
    f = lambda t: wa_from_dbwb(pa, db, t, over_ice=True) - wa
    wb, its, ok = brent(f, TMIN, t_upper_limit)
    return wb

  w_triple_liq = wa_from_dbwb(pa, db, TRIPLE, over_ice=False)
  if wa >= w_triple_liq:
    if abs(t_upper_limit - TRIPLE) < 1e-9:
      return TRIPLE
    f = lambda t: wa_from_dbwb(pa, db, t, over_ice=False) - wa
    wb, its, ok = brent(f, TRIPLE, t_upper_limit)
    return wb

  w_triple_ice = wa_from_dbwb(pa, db, TRIPLE, over_ice=True)
  if wa >= w_triple_ice:
    return TRIPLE

  f = lambda t: wa_from_dbwb(pa, db, t, over_ice=True) - wa
  wb, its, ok = brent(f, TMIN, TRIPLE)
  return wb

# ---------------------------------------------------------

def pv_from_wa(pa, wa):
  """ Vapor pressure from humidity.
  """
  pv = (pa * wa) / (MWA + wa)
  return pv / fs

def wa_from_pv(pa, pv):
  """ Humidity from vapor pressure.
  """
  p = fs * pv
  return MWA * p / (pa - p)

# =========================================================
# Dew-Point / Frost-Point Temperature
# =========================================================

def td_from_wa(pa, wa):
  """ Calculate dew-point (or frost-point below 0°C) 
      temperature [°C] given barometric pressure [kPa] 
      and humidity ratio [kg_w/kg_da].
  """
  if wa < 0.0:
    raise ValueError(f"wa negative!")

  if wa == 0.0:
    return TMIN

  # Calculate water vapour partial pressure (kPa)
  pv = pv_from_wa(pa, wa)

  # Invert psat(t) = pw safely using Brent solver
  return tsat(pv)

# =========================================================
# Relative Humidity [fraction]
# =========================================================

def rh_from_dbwa(pa, db, wa):
  """ Calculate relative humidity [fraction] from 
      barometric pressure [kPa], 
      dry-bulb temperature [°C] and 
      humidity ratio [kg_w/kg_da].
  """
  if wa < 0.0:
    raise ValueError(f"wa negative!")

  if wa == 0.0:
    return 0.0

  # Water vapour partial pressure [kPa]
  pv = pv_from_wa(pa, wa) 

  # Saturation pressure at dry-bulb temperature [kPa]
  ps = psat(db)

  return pv / ps

def wa_from_dbrh(pa, db, rh):
  """ Calculate humidity ratio [kg_w/kg_da] from 
      barometric pressure [kPa], 
      dry-bulb temperature [°C] and 
      relative humidity [-].
  """
  if not (0.0 <= rh <= 1.0):
    raise ValueError(f"rh out of range [0..1]!")

  if rh <= 0.0:
    return 0.0

  # Water vapour partial pressure [kPa]
  pv = rh * psat(db)

  if pv >= pa:
    raise ValueError(f"pv >= pa!")

  # Convert partial pressure to humidity ratio in [kg/kg]
  return wa_from_pv(pa, pv)

# =========================================================
# Moist Air Enthalpy (kJ/kg_da)
# =========================================================

def ha_from_dbwa(db, wa):
  """ Calculate moist air specific enthalpy [kJ/kg_da] 
      referenced to saturated liquid water at 0°C 
      (ASHRAE Ch. 1 standard datum).
  """
  if wa < 0.0:
    raise ValueError("wa negative!")
  
  if not (TMIN <= db <= TMAX):
    raise ValueError(f"db out of range [{TMIN}, {TMAX}]!")

  return CPda * db + wa * (Hfgo + CPv * db)

def ha_ice_from_dbwa(db, wa):
  """ Calculate moist air specific enthalpy [kJ/kg_da] 
      referenced to solid ice at 0°C 
      (for sub-zero frost/sublimation energy balances).
  """
  if wa < 0.0:
    raise ValueError(f"wa negative!")
  
  if not (TMIN <= db <= TMAX):
    raise ValueError(f"db out of range [{TMIN}, {TMAX}]!")

  return CPda * db + wa * (Hsgo + CPv * db)

def db_from_hawa(ha, wa):
  """ Calculate dry-bulb temperature [°C]
      from liquid-referenced enthalpy [kJ/kg_da] and 
      humidity ratio [kg_w/kg_da].
  """
  if wa < 0.0:
    raise ValueError(f"wa nagative!")

  db = (ha - wa * Hfgo) / (CPda + wa * CPv)

  if not (TMIN <= db <= TMAX):
    raise ValueError(f"db out of range [{TMIN}, {TMAX}]!")

  return db

def db_from_ha_ice_wa(ha_ice, wa):
  """ Calculate dry-bulb temperature [°C] 
      from ice-referenced enthalpy [kJ/kg_da] and 
      humidity ratio [kg_w/kg_da].
  """
  if wa < 0.0:
    raise ValueError(f"wa negative!")

  db = (ha_ice - wa * Hsgo) / (CPda + wa * CPv)

  if not (TMIN <= db <= TMAX):
    raise ValueError(f"db out of range [{TMIN}, {TMAX}]!")

  return db

def sigma_heat(db, wb, wa):
  """ Sigma heat [kJ/kg]. 
      sg = h_moist_air - h_liquid_water_at_wb
      Sigma heat is primarily used in mining engineering 
      to regulate and calculate mine air temperatures and 
      heat transfer in underground ventilation systems.  
      za = (cpa * wb) + (ws * hfg)   
  """

  if wb >= 0:
    return CPda * db + wa * (Hfgo + CPv * db - CPw * wb)
  
  # below freezing
  return CPda * db + wa * (Hsgo + CPv * db - CPi * wb)

def za_from_wb(pa, wb):
  """ Calculate Sigma Heat [kJ/kg dry air] 
      from pa [kPa] and wb [°C].
  """
  # Saturation humidity ratio [kg/kg] at wb
  wsb = wsat(pa, wb)

  # Latent heat of vaporization/sublimation [kJ/kg] at wb
  if wb >= 0:
    hfg = Hfgo - (CPw - CPv) * wb
  else:
    hfg = Hsgo - (CPi - CPv) * wb
  
  # Sigma heat [kJ/kg dry air]
  za = (CPda * wb) + (wsb * hfg)
  
  return za

# ---------------------------------------------------------

def cpmoist(wa):
  """ Specific heat of moist air [kJ/kg-K].
  """
  return CPda + CPv * wa

def db_from_hawa(ha, wa):
  
  if wa < 0.0:
    raise ValueError(f"wa negative!")

  # Enthalpy threshold at db = 0°C for the liquid branch
  h_zero_liq = wa * Hfgo
  
  if ha > h_zero_liq:
    db = (ha - wa * Hfgo) / (CPda + wa * CPv)
  else:
    # Note: If ha falls in the jump gap 
    # [wa * Hfgo, wa * Hsgo], clamping db to 0.0°C 
    # prevents a positive db from the ice formula.
    db = min(0.0, (ha - wa * Hsgo) / (CPda + wa * CPv))

  if not (TMIN <= db <= TMAX):
    raise ValueError(f"db out of range [{TMIN}, {TMAX}]!")

  return db

# =========================================================
# Specific Volume [m³/kg_da]
# =========================================================

def va_from_dbwa(pa, db, wa):
  """ Calculate moist air specific volume [m³/kg_da] 
      from barometric pressure [kPa], 
      dry-bulb temperature [°C] and 
      humidity ratio [kg_w/kg_da].
  """
  if wa < 0.0:
    raise ValueError(f"wa negative!")

  if not (TMIN <= db <= TMAX):
    raise ValueError(f"db out of range [{TMIN}, {TMAX}]!")

  # Convert temperature to [Kelvin]
  tk = db + K0

  # ASHRAE Ch. 1: v = R_da * T * (1 + 1.607858 * W) / P
  # Note: RA in [kJ/(kg-K] and pa in [kPa], 
  # yielding [m³/kg_da] directly
  return (RA * tk * (1.0 + (wa / MWA))) / pa

def wa_from_dbva(pa, db, va):
  """ Calculate humidity ratio [kg_w/kg_da] 
      from barometric pressure [kPa], 
      dry-bulb temperature [°C] and 
      specific volume [m³/kg_da].
  """

  if va <= 0.0:
    raise ValueError(f"va <= 0!")

  if not (TMIN <= db <= TMAX):
    raise ValueError(f"db out of range [{TMIN}, {TMAX}]!")

  tk = db + K0
  wa = MWA * ((va * pa) / (RA * tk) - 1.0)

  if wa < 0.0:
    raise ValueError(f"wa negative!")

  return wa

def ro_from_wava(wa, va):
  """ True density of moist air 
      [kg (DA + WV) / m³]
  """
  return (1.0 + wa) / va

# ---------------------------------------------------------

def db_from_ws(pa, ws, db0=20):
  """ tsat [°C] from wsat [kg/kg].
      Main purpose is to locate the 
      top point of the sat curve.
  """
  pv = pv_from_wa(pa, ws)
  db, its, ok = brent(lambda x: psat(x)-pv, -100.0, 200.0)
  return db

# =========================================================
# TAir Dataclass
# =========================================================

keys = "pa,db,wb,ha,wa,rh,td,pv,ps,va,ro,cp,za".split(',')
units = "kPa,°C,°C,kJ/kg,g/kg,%,°C,kPa,kPa,m³/kg,kg/m³,kJ/kg-K,kJ/kg".split(',')

@dataclass(slots=True)
class TAir:
  pa: float = 101.325     # kPa
  db: float = 23.0        # °C
  wb: float = 16.0        # °C
  ha: float = 44.7861     # kJ/kg
  wa: float = 0.00851405  # kg/kg
  rh: float = 0.484721    # -
  td: float = 11.5586     # °C
  pv: float = 1.36228     # kPa
  ps: float = 2.81044     # kPa
  va: float = 0.850487    # m³/kg
  ro: float = 1.18581     # kg/m³
  cp: float = 1.0215      # kJ/kg-K
  za: float = 44.2159     # kJ/kg

  UNITS: ClassVar[tuple[str, ...]] = (
    "kPa", "°C", "°C", "kJ/kg", "kg/kg", "-", "°C",
    "kPa", "kPa", "m³/kg", "kg/m³", "kJ/kg-K", "kJ/kg")

  REPORTUNITS: ClassVar[tuple[str, ...]] = (
    "kPa", "°C", "°C", "kJ/kg", "g/kg", "%", "°C",
    "kPa", "kPa", "m³/kg", "kg/m³", "kJ/kg-K", "kJ/kg")

  REPORTFACTOR: ClassVar[tuple[float, ...]] = (
    1.0, 1.0, 1.0, 1.0, 1e3, 1e2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

  def __getitem__(self, key: str) -> float:
    return getattr(self, key)

  def update(self, values: dict[str, float]) -> None:
    for key, value in values.items():
      setattr(self, key, value)

  def is_close(self, other: "TAir", rel_tol: float = 1e-6, abs_tol: float = 1e-9) -> bool:
    if not isinstance(other, TAir):
      return False
    return all(
      math.isclose(getattr(self, slot), getattr(other, slot), rel_tol=rel_tol, abs_tol=abs_tol)
      for slot in self.__slots__)

  def report(self, title: str = "Moist Air State") -> str:
    lines = [f"{title:^20}"]
    for key, unit, factor in zip(self.__slots__, self.REPORTUNITS, self.REPORTFACTOR, strict=True):
      lines.append(f"{key} = {factor * getattr(self, key):0.6g} {unit}")
    return "\n".join(lines)

  def short_report(self) -> str:
    return f"{self.db:8.1f} {self.wb:8.1f} {self.ha:8.1f} {1e3*self.wa:8.1f} {1e2*self.rh:8.2f} {self.td:8.1f}"

  def canvas_report(self) -> str:
    return "\n".join(
      f"{key}: {factor * getattr(self, key):8.4f} {unit}"
      for key, unit, factor in zip(self.__slots__, self.REPORTUNITS, self.REPORTFACTOR, strict=True))

# =========================================================
# Psychrometric Solvers
# =========================================================

def astate(pa, db, wb):
  """All moist air properties directly from (pa, db, wb).
  """
  ps = psat(db)
  wa = wa_from_dbwb(pa, db, wb)
  pv = pv_from_wa(pa, wa)
  ha = ha_from_dbwa(db, wa)
  rh = rh_from_dbwa(pa, db, wa)
  td = td_from_wa(pa, wa)
  va = va_from_dbwa(pa, db, wa)
  ro = ro_from_wava(wa, va)
  cp = cpmoist(wa)
  za = sigma_heat(db, wb, wa)
  return TAir(pa, db, wb, ha, wa, rh, td, pv, ps, va, ro, cp, za)

def astate2(pa, db, wa):
  """ All moist air properties 
      directly from (pa, db, wa).
  """
  ps = psat(db)
  wb = wb_from_dbwa(pa, db, wa)
  pv = pv_from_wa(pa, wa)
  ha = ha_from_dbwa(db, wa)
  rh = rh_from_dbwa(pa, db, wa)
  td = td_from_wa(pa, wa)
  va = va_from_dbwa(pa, db, wa)
  ro = ro_from_wava(wa, va)
  cp = cpmoist(wa)
  za = sigma_heat(db, wb, wa)
  return TAir(pa, db, wb, ha, wa, rh, td, pv, ps, va, ro, cp, za)

def astate3(pa, ha, wa):
  """ All moist air properties 
      directly from (pa, ha, wa).
  """
  db = db_from_hawa(ha, wa)
  ps = psat(db)
  wb = wb_from_dbwa(pa, db, wa)
  pv = pv_from_wa(pa, wa)
  rh = rh_from_dbwa(pa, db, wa)
  td = td_from_wa(pa, wa)
  va = va_from_dbwa(pa, db, wa)
  ro = ro_from_wava(wa, va)
  cp = cpmoist(wa)
  za = sigma_heat(db, wb, wa)
  return TAir(pa, db, wb, ha, wa, rh, td, pv, ps, va, ro, cp, za)

def calcstate(pa=101.325, db0=23.0, wb0=16.0, **targets):
  """ Solve any two independent moist air properties 
      using fsolve. db0 and wb0 are starting guesses.
  """

  valid_keys = {"db", "wb", "ha", "wa", "rh", "td", "pv", "ps", "va", "ro", "cp", "za"}
  if len(targets) != 2 or not set(targets).issubset(valid_keys):
    raise ValueError(f"Provide exactly two valid properties from {sorted(valid_keys)}!")

  invalid_combs = {frozenset({"wb", "ha"}), 
    frozenset({"wb", "za"}), frozenset({"wa", "td"}),
    frozenset({"td", "pv"}), frozenset({"td", "cp"})}
  target_keys = frozenset(targets.keys())
  if target_keys in invalid_combs:
    raise ValueError(f"Invalid combination {tuple(targets.keys())}: properties are not independent!")

  k0, k1 = list(targets.keys())
  v0, v1 = float(targets[k0]), float(targets[k1])

  # 1. Estimate Dry Bulb (db) from available targets
  if "db" in targets:
    db_guess = targets["db"]
  elif "wb" in targets:
    # Dry bulb is always >= wet bulb; add a small buffer (e.g., +2°C)
    db_guess = targets["wb"] + 2.0
  elif "td" in targets:
    # Dry bulb is always >= dew point; add a buffer (e.g., +5°C)
    db_guess = targets["td"] + 5.0
  elif "ha" in targets:
    # Quick rough estimate: h ≈ 1.006*db + w*(2501 + 1.86*db)
    # Assuming ~50% RH, db is roughly ha / 2.5 (in °C and kJ/kg)
    db_guess = max(0.0, targets["ha"] / 2.5)
  else:
    db_guess = db0  # Fallback (e.g., 23.0°C)

  # 2. Estimate Wet Bulb (wb) from available targets
  if "wb" in targets:
    wb_guess = targets["wb"]
  elif "td" in targets:
    # Wet bulb always lies between dew point and dry bulb
    wb_guess = 0.5 * (db_guess + targets["td"])
  elif "rh" in targets:
    # Rough linear approximation: at 100% RH, wb = db; at 0% RH, wb ≈ db - 10
    wb_guess = db_guess - (1.0 - targets["rh"]) * 10.0
  else:
    wb_guess = min(wb0, db_guess - 1.0)

  # 3. Enforce physical constraints: db >= wb (and db >= td if known)
  if "td" in targets:
    db_guess = max(db_guess, targets["td"] + 0.5)
    wb_guess = max(wb_guess, targets["td"] + 0.2)

  wb_guess = min(wb_guess, db_guess)

  # Characteristic positive spans for each property
  scale_map = {
    "db": 20.0, "wb": 15.0, "ha": 40.0, "wa": 8.0e-3,
    "rh": 0.50, "td": 15.0, "pv": 1.5,  "ps": 2.5,
    "va": 0.85, "ro": 1.2,  "cp": 1.0,  "za": 40.0}

  # Safeguard against zero targets while adapting to large target magnitudes
  s0 = max(scale_map[k0], abs(v0))
  s1 = max(scale_map[k1], abs(v1))

  def residual(x):
    st = astate(pa, x[0], x[1])
    return (st[k0] - v0) / scale_map[k0], (st[k1] - v1) / scale_map[k1]

  # db, wb = fsolve(residual, guess=(db_guess, wb_guess),
  #   bounds=((TMIN + 1e-3, TMAX - 1e-3), (0.0, None)),
  #   scales=(scale_map[k0], scale_map[k1]))
  (db, wb), it_c, ok_c = fsolve(residual, (db_guess, wb_guess))

  air = astate(pa, db, wb)
  return air

def frange(start: float, stop: float, step: float):
  """ Yield floating-point numbers 
      from start up to stop (inclusive).
  """
  steps = int(round((stop - start) / step))
  for i in range(steps + 1):
    yield round(start + i * step, 10)

# =========================================================
# Psychrometric Chart Plotter
# =========================================================

# Define constant property line colors
colors = {
  "sat": "indian red",
  "db": "cadet blue",
  "wb": "dark sea green",
  "ha": "dark sea green",
  "wa": "cadet blue",
  "rh": "light gray",
  "va": "light gray",
  "za": "dark sea green",
}

def plotchart(**kwargs):
  """Plot toy psychrometric chart using tkinter.
  """
  import tkinter as tk
  import os
  if os.name == "nt":
    import ctypes
    try:
      ctypes.windll.shcore.SetProcessDpiAwareness(2)
    except Exception:
      try:
        ctypes.windll.user32.SetProcessDPIAware()
      except Exception:
        pass

  width, height = 1600, 900
  padx, pady = 80, 80
  x0, x1 = padx, width - padx
  y0, y1 = height - pady, pady

  pa = kwargs.get('pa', 101.325)
  db0, db1 = kwargs.get('db_range', (-10, 45))
  wa0, wa1 = kwargs.get('wa_range', (0, 0.030))

  dtdx = (db1 - db0) / (x1 - x0)
  dwdy = (wa1 - wa0) / (y1 - y0)

  # 1. Define sheared horizontal coordinate u(db, wa)
  h_bot_max = ha_from_dbwa(db1, wa0)
  h_top_max = ha_from_dbwa(db1, wa1)
  k_shear = (h_top_max - h_bot_max) / (wa1 - wa0)


  def toU(db, wa):
    return ha_from_dbwa(db, wa) - k_shear * wa

  u0 = toU(db0, wa0)  # Bottom-left corner (db0, wa0)
  u1 = toU(db1, wa0)  # Bottom-right corner (db1, wa0)
  dudx = (u1 - u0) / (x1 - x0)
  dwdy = (wa1 - wa0) / (y1 - y0)

  def toX(db, wa):
    return x0 + (toU(db, wa) - u0) / dudx

  def toY(wa):
    return y0 + (wa - wa0) / dwdy

  def toW(y):
    return wa0 + (y - y0) * dwdy

  # 2. Screen / chart conversions (db, wa) <-> (x, y)

  def toScreen(db, wa):
    return toX(db, wa), toY(wa)

  def toChart(x, y):
    wa = toW(y)
    u = u0 + (x - x0) * dudx
    ha = u + k_shear * wa
    db = db_from_hawa(ha, wa)
    return db, wa

  def draw_small_point(db, wa):
    x, y = toScreen(db, wa)
    r = 5
    coords = x-r, y-r, x+r, y+r
    canvas.create_oval(coords, fill="pink")

  def draw_point(db, wa):
    x, y = toScreen(db, wa)
    r = 7
    canvas.create_line(x, y-2*r, x, y+2*r)
    canvas.create_line(x-2*r, y, x+2*r, y)
    coords = x-r, y-r, x+r, y+r
    canvas.create_oval(coords, fill="light yellow")

  def draw_process(db1, wa1, db2, wa2):
    x1, y1 = toScreen(db1, wa1)
    x2, y2 = toScreen(db2, wa2)
    canvas.create_line(x1, y1, x2, y2, 
      width=3, fill="green")

  def move_mouse(event):
    canvas.delete("label")
    x, y = event.x, event.y
    if x < x0 or x > x1 or y < y1 or y > y0:
      return
    db, wa = toChart(x, y)
    try:
      air = astate2(pa, db, wa)
      color = "grey" if air["wb"] >= air["db"] else "black"
      text_id = canvas.create_text(padx+5, pady+5, 
        text=air.canvas_report(), 
        fill=color, anchor="nw", 
        font=("Consolas", 12), tags="label")
      xa, ya, xb, yb = canvas.bbox(text_id)
      pad = 5
      box_coords = xa-pad, ya-pad, xb+pad, yb+pad
      rect_id = canvas.create_rectangle(box_coords, 
        fill="white", outline="#006064", 
        width=1, tags="label")
      canvas.tag_lower(rect_id, text_id)
    except:
      ...

  def left_click(event):
    x, y = event.x, event.y
    if x < x0 or x > x1 or y < y1 or y > y0:
      return
    db, wa = toChart(x, y)
    draw_point(db, wa)


  def draw_db_lines():
    """ Draw constant dry-bulb lines 
        in sheared (h, W) space.
    """
    db_step = 5
    db_start = math.ceil(db0 / db_step) * db_step
    db_end = math.floor(db1 / db_step) * db_step
    for db in frange(db_start, db_end, db_step):
      w_top = min(wa1, wsat(pa, db))
      if w_top <= wa0:
        continue
      # Two endpoints are exact for db > 0 (or db <= 0) since h is linear in W
      p_bot = toScreen(db, wa0)
      p_top = toScreen(db, w_top)
      canvas.create_line(*p_bot, *p_top, 
        fill=colors["db"], width=1)
      canvas.create_text(p_bot[0], y0 + 12, 
        text=f"{db:g}", fill="black", 
        font=("Segoe UI", 8))

  def draw_wa_lines():
    """ Draw constant abs humidity lines [kg_w/kg_da], 
        labeled in [g/kg].
    """
    # Scale step automatically in kg/kg (e.g. 0.005 for 30 g/kg, 0.010 for 100 g/kg)
    span_g = (wa1 - wa0) * 1000.0
    if span_g <= 5.0:
      wa_step = 0.0005   # 0.5 g/kg for low-temp charts
    elif span_g <= 40.0:
      wa_step = 0.005    # 5 g/kg for normal/med-temp charts
    else:
      wa_step = 0.010    # 10 g/kg for high-temp charts

    wa_start = math.ceil(wa0 / wa_step) * wa_step
    wa_end = math.floor(wa1 / wa_step) * wa_step

    for wa in frange(wa_start, wa_end, wa_step):
      # Left end stops at the saturation curve unless clipped by db0
      db_left = db0 if wa <= wsat(pa, db0) else db_from_ws(pa, wa)
      if db_left >= db1:
        continue
      p_left = toScreen(db_left, wa)
      p_right = toScreen(db1, wa)
      canvas.create_line(*p_left, *p_right,
        fill=colors["wa"], width=1)
      if wa > wa0:
        canvas.create_text(
          x1 + 6, p_right[1], text=f"{1e3*wa:g}",
          fill="black", anchor="w", 
          font=("Segoe UI", 8))

  def draw_wb_lines():
    """ Draw constant wet bulb lines across range 
        (db0, wa0) to (db1, wa1).
    """
    wb0 = wb_from_dbwa(pa, db0, wa0)
    wb1 = wb_from_dbwa(pa, db1, wa1)
    wb_step = 5
    wb_start = math.ceil(wb0 / wb_step) * wb_step
    wb_end = math.floor(wb1 / wb_step) * wb_step    
    db_step = 0.5
    for wb in frange(wb_start, wb_end, wb_step):
      points = []
      db = max(wb, db0)
      while db <= db1:
        wa = wa_from_dbwb(pa, db, wb)
        if wa0 <= wa <= wa1:
          points.extend(toScreen(db, wa))
        elif wa < wa0:
          w_prev = wa_from_dbwb(pa, db - db_step, wb)
          if w_prev > wa0:
            t_zero = (db - db_step) + db_step * ((w_prev - wa0) / (w_prev - wa))
            points.extend(toScreen(t_zero, wa0))
          break
        db += db_step
      if len(points) >= 4:
        canvas.create_line(points, fill=colors["wb"], 
          width=1, dash=(4, 3), smooth=True)
        start_x, start_y = points[0], points[1]

        # Draw label only if the line actually starts on the saturation curve inside the chart
        if db0 <= wb <= db1 and wa0 <= wa_from_dbwb(pa, wb, wb) <= wa1:
          start_x, start_y = points[0], points[1]
          canvas.create_text(start_x - 4, start_y - 2,
            text=f"{wb:g}°C", fill=colors["wb"], 
            anchor="se",
            font=("Segoe UI", 8))

  def draw_rh_lines():
    """ Draw constant relative humidity curves.
    """
    db_step = 0.5
    for rh in [0.05, 0.25, 0.50, 0.75]:
      points = []
      db = db0
      while db <= db1 + 1e-9:
        wa = wa_from_dbrh(pa, db, rh)
        if wa <= wa1:
          if wa >= wa0:
            points.extend(toScreen(db, wa))
        else:
          # Interpolate exact db where this RH curve exits 
          # the top border (wa1)
          if db > db0:
            db_top, its, ok = brent(lambda t: wa_from_dbrh(pa, t, rh) - wa1, db - db_step, db)
            points.extend(toScreen(db_top, wa1))
          break
        db += db_step

      if len(points) >= 4:
        canvas.create_line(points, fill=colors["rh"], width=1, smooth=True)
        # Place RH label near the end of the curve
        end_x, end_y = points[-2], points[-1]
        if end_x < x1 - 15:
          canvas.create_text(
            end_x, end_y - 4, text=f"{100*rh}%", fill=colors["rh"],
            anchor="s", font=("Segoe UI", 8))
        else:
          canvas.create_text(
            end_x - 4, end_y - 2, text=f"{100*rh}%", fill=colors["rh"],
            anchor="se", font=("Segoe UI", 8))

  def draw_ha_lines():
    """ Draw constant specific enthalpy lines.
    """
    def raw_db(ha, wa):
      return (ha - wa * Hfgo) / (CPda + wa * CPv)

    ha_min = ha_from_dbwa(db0, wa0)
    ha_max = ha_from_dbwa(db1, wa1)
    ha_step = 10
    ha_start = math.ceil(ha_min / ha_step) * ha_step
    ha_end = math.floor(ha_max / ha_step) * ha_step

    for ha in frange(ha_start, ha_end, ha_step):
      # 1. Lower-right endpoint: bottom axis (wa0) or right axis (db1)
      db_bot = raw_db(ha, wa0)
      if db_bot < db0:
        continue
      if db_bot <= db1:
        db_start, wa_start = db_bot, wa0
      else:
        wa_right = (ha - CPda * db1) / (Hfgo + CPv * db1)
        if not (wa0 <= wa_right <= wa1):
          continue
        db_start, wa_start = db1, wa_right

      # 2. Upper-left endpoint: top axis (wa1) or saturation curve
      db_top = raw_db(ha, wa1)
      if db_top >= db0 and wa1 <= wsat(pa, min(db_top, db1)):
        if db_top > db1:
          continue
        db_end, wa_end = db_top, wa1
      else:
        # Intersect constant-h line with saturation curve
        t_low = db0
        t_high = min(db1, db_from_ws(pa, wa1))
        if ha_from_dbwa(t_low, wsat(pa, t_low)) <= ha <= ha_from_dbwa(t_high, wsat(pa, t_high)):
          db_end, its, ok = brent(lambda t: ha_from_dbwa(t, wsat(pa, t)) - ha, t_low, t_high)
          wa_end = wsat(pa, db_end)
        else:
          continue

      p1 = toScreen(db_start, wa_start)
      p2 = toScreen(db_end, wa_end)
      canvas.create_line(*p1, *p2, fill=colors["ha"], width=1, dash=(2, 2))
      canvas.create_text(p2[0]-4, p2[1]-2, text=f"{ha:g}", fill=colors["ha"],
        anchor="se", font=("Segoe UI", 8))

  def draw_va_lines():
    """ Draw constant specific volume lines [m³/kg_da].
    """
    va_min = va_from_dbwa(pa, db0, wa0)
    va_max = va_from_dbwa(pa, db1, wa1)
    va_step = 0.02
    va_start = math.ceil(va_min / va_step) * va_step
    va_end = math.floor(va_max / va_step) * va_step

    for va in frange(va_start, va_end, va_step):
      # 1. Lower-right endpoint: where constant-va meets wa0 or db1
      # At wa = wa0: va = RA * (db + K0) * (1 + wa0/MWA) / pa
      db_bot = (va * pa) / (RA * (1.0 + wa0 / MWA)) - K0
      if db_bot <= db1:
        if db_bot < db0:
          continue
        db_start, wa_start = db_bot, wa0
      else:
        wa_right = wa_from_dbva(pa, db1, va)
        if not (wa0 <= wa_right <= wa1):
          continue
        db_start, wa_start = db1, wa_right

      # 2. Upper-left endpoint: 
      #    where constant-va meets top border (wa1) 
      #    or saturation curve
      db_top = (va * pa) / (RA * (1.0 + wa1 / MWA)) - K0
      if db_top >= db0 and wa1 <= wsat(pa, db_top):
        db_end, wa_end = db_top, wa1
      else:
        # Intersect constant-va curve with saturation curve: 
        # wa_from_dbva(pa, t, va) == wsat(pa, t)
        t_low = max(db0, (va * pa) / (RA * (1.0 + wa1 / MWA)) - K0)
        t_high = min(db1, db_bot)
        if t_low < t_high:
          f_sat = lambda t: wa_from_dbva(pa, t, va) - wsat(pa, t)
          if f_sat(t_low) * f_sat(t_high) <= 0:
            db_end, its, ok = brent(f_sat, t_low, t_high)
            wa_end = wsat(pa, db_end)
          else:
            continue
        else:
          continue

      # Sample points along the curve 
      # (slightly curved in sheared h-W coordinates)
      points = []
      n_pts = 15
      for i in range(n_pts + 1):
        try:
          db = db_start + (db_end - db_start) * i / n_pts
          wa = wa_from_dbva(pa, db, va)
          points.extend(toScreen(db, wa))
        except:
          ...

      if len(points) >= 4:
        canvas.create_line(points, fill=colors["va"], width=1, smooth=True)
        # Place label slightly above the midpoint of the line
        mid_idx = (len(points) // 4) * 2
        mx, my = points[mid_idx], points[mid_idx + 1]
        canvas.create_text(mx+3, my-3, text=f"{va:.2f}", fill=colors["va"],
          anchor="sw", font=("Segoe UI", 8))
  
  def draw_za_lines():
    """ Draw constant sigma heat lines [kJ/kg_da].
    """
    
    # Sigma heat depends 
    # strictly on wet-bulb temperature 
    # (za = za_from_wb(pa, wb))
    wb_min = wb_from_dbwa(pa, db0, wa0)
    wb_max = wb_from_dbwa(pa, db1, wa1)
    za_min = za_from_wb(pa, wb_min)
    za_max = za_from_wb(pa, wb_max)

    za_step = 10
    za_start = math.ceil(za_min / za_step) * za_step
    za_end = math.floor(za_max / za_step) * za_step

    for za in frange(za_start, za_end, za_step):
      # 1. Solve for the corresponding wet-bulb temperature on the saturation curve
      wb, its, ok = brent(lambda t: za_from_wb(pa, t) - za, wb_min - 1.0, wb_max + 1.0)

      # 2. Upper-left endpoint: saturation curve (db = wb) or left border (db0)
      if wb >= db0:
        db_start, wa_start = wb, wsat(pa, wb)
      else:
        db_start = db0
        wa_start = wa_from_dbwb(pa, db0, wb)

      # Skip if starting point is outside the chart's humidity bounds
      if not (wa0 <= wa_start <= wa1) or db_start >= db1:
        continue

      # 3. Lower-right endpoint: bottom axis (wa0) or right axis (db1)
      # From za = CPda * db + wa0 * (Hfgo + CPv * db - CPw * wb)  [or ice branch for wb < 0]
      if wb >= 0.0:
        db_bot = (za - wa0 * (Hfgo - CPw * wb)) / (CPda + wa0 * CPv)
      else:
        db_bot = (za - wa0 * (Hsgo - CPi * wb)) / (CPda + wa0 * CPv)

      if db_bot <= db1:
        db_end, wa_end = db_bot, wa0
      else:
        db_end = db1
        wa_end = wa_from_dbwb(pa, db1, wb)
        if not (wa0 <= wa_end <= wa1):
          continue

      p1 = toScreen(db_start, wa_start)
      p2 = toScreen(db_end, wa_end)
      canvas.create_line(*p1, *p2, fill=colors["za"], width=1, dash=(6, 3))

      # Label along the saturation curve if inside the chart bounds
      if db0 <= wb <= db1 and wa0 <= wa_start <= wa1:
        canvas.create_text(
          p1[0] - 4, p1[1] - 2,
          text=f"{za:g}",
          fill=colors["za"],
          anchor="se",
          font=("Segoe UI", 8))

  root = tk.Tk()
  root.title("PsyChart - true energy vs mass plot")
  root.geometry(f"{width}x{height}+20+20") 
  canvas = tk.Canvas(root, background="ivory")
  canvas.config(bd=0, highlightthickness=0)
  canvas.pack(expand=True, fill="both")
  canvas.bind("<Motion>", move_mouse)
  canvas.bind("<Button-1>", left_click)

  # Constant lines ordered in layers.
  draw_va_lines()
  draw_rh_lines()
  draw_ha_lines() # typically ha only, no wb or za lines.
  #draw_wb_lines()
  #draw_za_lines()
  draw_wa_lines()
  draw_db_lines() # draw db last, so on top

  # Box chart area.
  #canvas.create_rectangle(x0, y0, x1, y1, width=1, outline="grey")
  
  # Draw both axes.
  coords = x1, y1, x1, y0, x0, y0
  canvas.create_line(coords, width=2) 

  # Draw saturation curve.
  db3 = db_from_ws(pa, wa1)
  db = db0
  npoints = 19
  dt = (db3 - db0) / (npoints + 1)
  points = []
  while db <= db3 + 1e-3:
    wa = wsat(pa, db)
    x, y = toScreen(db, wa)
    points.extend((x,y))
    db += dt
  canvas.create_line(points, width=2, fill=colors["sat"], smooth=True)

  root.mainloop()
