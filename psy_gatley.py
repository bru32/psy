# psy_gatley.py
#
# Psychrometric routines
# The constants are based on those published by 
# Don Gatley "Understanding Psychrometrics"
# The root finders, astate, calcstate and dataclass
# definitions and coding structure my work.
#
# Author: Bruce Wernick
# last updated: 24 November 2026
#
# Copyright © (2026) Bruce Wernick,  All rights reserved...
# ---------------------------------------------------------


import sys
from math import hypot, exp, log as ln, isclose
from collections import namedtuple
from dataclasses import dataclass, fields
from typing import Optional
# ---------------------------------------------------------
from root_finder import secant
from root_finder_2d import dnewton2
# ---------------------------------------------------------

# Psychrometric Constants

p0 = 101.325    # barometric pressure at sea level [kPa]
k0 = 273.15     # degC to Kelvin offset
ps0 = 0.611213  # Saturation pressure at 0°C [kPa]
fs = 1.0044     # molecular interaction coeff used by Gatley
cpa = 1.006     # specific heat of dry air
cpv = 1.82      # specific heat of vapor (see Gatley page 228)
cpw = 4.186     # specific heat of water
cpi = 2.093     # specific heat of ice
hgo = 2500.77   # saturated enthalpy of gas at 0degC
hgi = 2834.2    # saturation enthalpy over ice at 0degC
sgo = hgo / k0  # saturated entropy of gas at 0degC
ma = 28.9645    # molecular weight of air
mw = 18.01528   # molecular weight of water
ru = 8.31441    # universal gas constant
ra = ru / ma    # air gas constant
rw = ru / mw    # water gas constant
raw = ra / rw   # gas constant ratio

# ---------------------------------------------------------

# Define constants once at the module level 
# so they aren't recreated on every call.

C_ICE = (-9.4840240e-13, 2.0747825e-9, 6.2215701e-7, 
         -9.6778430e-3, 6.3925247, -5.674359e3, 4.1635019)

C_WATER = (-1.4452093e-8, 4.1764768e-5, -4.8640239e-2, 
            1.3914993, -5.8002206e3, 6.5459673)

def psat(t):
  k = t + k0
  # Select the right tuple and unpack the log coefficient from the polynomial coefficients
  *poly_coeffs, log_coeff = C_WATER if t >= 0.0 else C_ICE

  # Evaluate the polynomial cleanly using Horner's method
  poly = 0.0
  for c in poly_coeffs:
    poly = poly * k + c

  return 1e-3 * exp(log_coeff * ln(k) + poly / k)

def tsat(p):
  """ASHRAE empirical inverse of psat p [kPa], t [degC]
  """
  a = ln(p)
  if p>=ps0:
    # Over Liquid Water (Valid for 0°C to 93°C)
    c = (6.54, 14.526, 0.7389, 0.09486, 0.4569)
    return c[0] + a*(c[1] + a*(c[2] + a*c[3])) + c[4]*p**0.1984
  if p <= 0:
    #raise ValueError("Sat pressure must be positive!")
    return 0
  else:
    # Over Ice / Frost (Valid for -60°C to 0°C)
    return 6.09 + a*(12.608 + a*0.4959)

def tdew(p):
  """Dew Point Temperature from pressure (degC).
     Curve fit by Bruce Wernick.
  """
  if p>=ps0:
    # liquid water (0 to 93 degC)
    a,b,c,d = 208.082023334, -6.53101910963, -0.0675213217071, 0.00101000996947
  elif p<=0:
    return 0
  else:
    # over ice (-60 to 0 degC)
    a,b,c,d = 212.569301264, 7.02161885813, -0.00152929927821, -0.00114243639171
  p = ln(1.0e3*p)
  return (a+b*p)/(1.0+p*(c+p*d))-k0

def tsat_inv(p):
  """Reverse ASHRAE psat equation using 
     Secant method with error function.
  """
  guess = tdew(p)
  func = lambda x: psat(x) - p
  ts = secant(func, guess, maxits=124, tol=1e-8)
  return ts

def wa_papv(pa,pv):
  """Absolute Humidity (kg/kg) 
     wa = R*p/(pa-p)
  """
  p = fs*pv
  return raw * p / (pa - p)

def humidity(pa, db, wb):
  """Calc wa from (db,wb) 
  Energy balance of adiabatic saturator.

  ... Extract from PSYCHROLIB ...
  ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 33 and 35
  wa = (wsw*(hgo+(cpw-cpv)*wb)-cpa*(db-wb))/(hgo+cpv*db-cpw*wb)
  """
  wsw = wa_papv(pa, psat(wb))
  if wb >= 0.0:
    num = wsw * (hgo + (cpv - cpw) * wb) - cpa * (db - wb)
    den = hgo + cpv * db - cpw * wb
  else:
    num = wsw * (hgi + (cpv - cpi) * wb) - cpa * (db - wb)
    den = hgi + cpv * db - cpi * wb
  return num / den

def pvap(pa,wa):
  """Vapor pressure (kPa), 
     p = wa*pa / (R + wa)
  """
  return wa*pa/(raw+wa)/fs

def enthalpy(db,wa):
  """Enthalpy of moist air (kJ/kg)
  """
  return cpa * db + wa * (hgo + cpv * db)
  #return cpa * db + wa * (hgi + cpv * db)

def sigheat(db,wb,wa):
  """Sigma heat (kJ/kg). 
     sg = h_moist_air - h_liquid_water_at_wb
  """
  if wb >= 0:
    return cpa*db + wa*(hgo + cpv*db - cpw*wb)
  return cpa*db + wa*(hgi + cpv*db - cpi*wb)

def specvol(pa,db,wa):
  """Specific volume based on perfect gas equation. 
     pv = RT
  """
  return ra*(db + k0) * (1.0 + wa*ma / mw) / pa

def trudens(va,wa):
  """True density of moist air (kg(DA+WV)/m3)
  """
  return (1.0+wa)/va

def cpmoist(wa):
  """Specific heat of moist air (kJ/kg-K)
  """
  return cpa+cpv*wa

def entropy(pa, db, pv=0):
  """Moist air entropy [kJ/kg-K]
  """
  T = db + k0

  # 1. Dry air entropy 
  p_da = pa - pv
  sd = cpa * ln(T / k0) - ra * ln((pa - pv) / p0)
  if pv == 0.0:
    # Completely dry air
    return sd

  # 2. Water vapor entropy
  sg = sgo + cpv * ln(T / k0) - raw * ln(pv / ps0)

  return sd + wa_papv(pa, pv) * sg  

# ---------------------------------------------------------

def db_sat(pa, ws, **kwargs):
  """Return db for given ws [degC]
     Uses secant root finder.
  """
  f = lambda x: 1e3*(ws - humidity(pa, x, x))
  db0 = kwargs.get("db0", 23)
  db = secant(f, db0)
  return db

def wetbulb(pa, db, wa, **kwargs):
  """wb from pa, db and wa.
     Give an initial guess for a better start
     wetbulb(101.325, 23, 0.095, wb=16)
     Uses secant root finder.
  """
  def f(x): 
    return wa - humidity(pa, db, x)
  def wb_guess():
    pv = pvap(pa, wa)
    td = tsat_inv(pv)
    return (db+td)/2
  wb = kwargs.pop("wb") if "wb" in kwargs else wb_guess()
  wb = secant(f, wb)
  return wb

# ---------------------------------------------------------

# Some structure definitions.

propstr = "pa,db,wb,ha,wa,rh,td,pv,ps,va,ro,cp,za"
airprops = propstr.split(",")
sfact = (1,1,1,1,1e3,1e2,1,1,1,1,1,1,1)
Air = namedtuple("Air", airprops)

# Moist Air State Meta Data
air_spec = {
  'pa': {'label': 'Barometric Pressure', 'name': 'Barom', 'value': 101.325, 'minvalue': 78, 'maxvalue': 119, 'ndigs': 3, 'sf': 1.0, 'units': 'kPa'}, 
  'db': {'label': 'Dry Bulb Temperature', 'name': 'DryBulb', 'value': 23.0, 'minvalue': -60, 'maxvalue': 90, 'ndigs': 3, 'sf': 1.0, 'units': 'degC'}, 
  'wb': {'label': 'Wet Bulb Temperature', 'name': 'WetBulb', 'value': 16.0, 'minvalue': -60, 'maxvalue': 90, 'ndigs': 3, 'sf': 1.0, 'units': 'degC'}, 
  'ha': {'label': 'Specific Enthalpy', 'name': 'Enthalpy', 'value': 44.786, 'minvalue': -60, 'maxvalue': 355, 'ndigs': 3, 'sf': 1.0, 'units': 'kJ/kg'}, 
  'wa': {'label': 'Absolute Humidity', 'name': 'Humidity', 'value': 8.514, 'minvalue': 0.0, 'maxvalue': 100, 'ndigs': 3, 'sf': 1000.0, 'units': 'g/kg'}, 
  'rh': {'label': 'Relative Humidity', 'name': 'RelHumid', 'value': 48.47, 'minvalue': 0.0, 'maxvalue': 100.0, 'ndigs': 3, 'sf': 100.0, 'units': '%'}, 
  'td': {'label': 'Dew Point Temperature', 'name': 'DewPoint', 'value': 11.558, 'minvalue': -60, 'maxvalue': 90, 'ndigs': 3, 'sf': 1.0, 'units': 'degC'}, 
  'pv': {'label': 'Vapor Pressure', 'name': 'VapPress', 'value': 1.36228, 'minvalue': 0.001, 'maxvalue': 12.4, 'ndigs': 3, 'sf': 1.0, 'units': 'kPa'}, 
  'ps': {'label': 'Saturation Pressure', 'name': 'SatPress', 'value': 2.81, 'minvalue': 0.001, 'maxvalue': 12.4, 'ndigs': 3, 'sf': 1.0, 'units': 'kPa'}, 
  'va': {'label': 'Specific Volume', 'name': 'SpecVol', 'value': 0.85048, 'minvalue': 0.514, 'maxvalue': 1.41, 'ndigs': 4, 'sf': 1.0, 'units': 'm3/kg'}, 
  'ro': {'label': 'True Moist Air Density', 'name': 'Density', 'value': 1.1858147, 'minvalue': 0.792, 'maxvalue': 1.95, 'ndigs': 4, 'sf': 1.0, 'units': 'kg/m3'}, 
  'cp': {'label': 'Specific Heat', 'name': 'SpecHeat', 'value': 1.02149, 'minvalue': 0.5, 'maxvalue': 5.0, 'ndigs': 4, 'sf': 1.0, 'units': 'kJ/kg-K'}, 
  'za': {'label': 'Sigma Heat', 'name': 'SigmaHeat', 'value': 44.21585, 'minvalue': -60, 'maxvalue': 330, 'ndigs': 3, 'sf': 1.0, 'units': 'kJ/kg'}, 
  'ws': {'label': 'Saturated Abs Humidity', 'name': 'SatHumidity', 'value': 17.8242, 'minvalue': 0, 'maxvalue': 117, 'ndigs': 3, 'sf': 1000.0, 'units': 'g/kg'}, 
  'ds': {'label': 'Degree of Saturation', 'name': 'DegSat', 'value': 0.4776678, 'minvalue': 0.0, 'maxvalue': 1.0, 'ndigs': 3, 'sf': 1.0, 'units': '-'}, 
  'sa': {'label': 'Specific Entropy', 'name': 'Entropy', 'value': 0.1601728, 'minvalue': -0.30, 'maxvalue': 1.25, 'ndigs': 3, 'sf': 1.0, 'units': 'kJ/kg-K'}
 }

# Complete air class with optional extra fields
@dataclass
class TAir:
  pa: float  # Barometric Pressure [kPa]
  db: float  # Dry Bulb Temperature [degC]
  wb: float  # Wet Bulb Temperature [degC]
  ha: float  # Specific Enthalpy [kJ/kg]
  wa: float  # Absolute Humidity [kg water vapor / kg dry air]
  rh: float  # Relative Humidity []
  td: float  # Dew Point Temperature [degC]
  pv: float  # Vapor Pressure [kPa]
  ps: float  # Saturation Pressure [kPa]
  va: float  # Specific Volume [m3/kg]
  ro: float  # True Moist Air Density [kg/m3]
  cp: float  # Specific Heat [kJ/kg-K]
  # --- Unique / Optional Variables (Must have defaults) ---
  za: Optional[float] = None  # Sigma Heat [kJ/kg]
  ws: Optional[float] = None  # Saturation Humidity [kg/kg]
  ds: Optional[float] = None  # Degree of Saturation []
  sa: Optional[float] = None  # Moist Air Entropy [kJ/kg-K]

  def to_spec(self) -> dict:
    """Returns a populated copy of air_spec for all non-None fields."""
    return {
        key: {**meta, "value": getattr(self, key) * meta["sf"]}
        for key, meta in air_spec.items()
        if getattr(self, key) is not None
    }

  def report(self, title: str = "Moist Air Report (Gatley)") -> None:
    """Prints a cleanly aligned report using air_spec metadata."""
    print(f"\n{title}")
    print("-" * 48)
    for key, meta in self.to_spec().items():
      val_str = f"{meta['value']:.{meta['ndigs']}f}"
      unit_str = "" if meta["units"] == "-" else f"[{meta['units']}]"
      print(f"{meta['label']:<24} {val_str:>10} {unit_str}")
    print("-" * 48)

  def is_close(self, other: "TAir", rel_tol: float = 1e-6, abs_tol: float = 1e-9) -> bool:
    if not isinstance(other, TAir):
      return False

    for f in fields(self):
      v1 = getattr(self, f.name)
      v2 = getattr(other, f.name)

      # If one is None and the other isn't, they don't match
      if (v1 is None) != (v2 is None):
        return False

      # If both are numbers (not None), compare them with tolerance
      if v1 is not None and not isclose(v1, v2, rel_tol=rel_tol, abs_tol=abs_tol):
        return False

    return True

# ---------------------------------------------------------

def astate(pa, db, wb):
  """Returns the primary properties from (pa,db,wb)
  Result is a tuple in a critical order.  The reason is 
  because it is used in a 2D root finder to calculate
  the error function for any given property pair.
  """
  ps = psat(db)
  wa = humidity(pa, db, wb)
  ha = enthalpy(db, wa)
  pv = pvap(pa, wa)
  rh = pv / ps
  td = tsat_inv(pv)
  va = specvol(pa, db, wa)
  ro = trudens(va, wa)
  cp = cpmoist(wa)
  za = sigheat(db, wb, wa)
  return (pa, db, wb, ha, wa, rh, td, pv, ps, va, ro, cp, za)

def astate2(pa, db, wb):
  """Calculate extended properties from (pa,db,wb)
     NOT for use by calcstate
  """
  ps = psat(db)
  wa = humidity(pa, db, wb)
  ha = enthalpy(db, wa)
  pv = pvap(pa, wa)
  rh = pv / ps
  td = tsat_inv(pv)
  va = specvol(pa, db, wa)
  ro = trudens(va, wa)
  cp = cpmoist(wa)
  za = sigheat(db, wb, wa)
  ws = humidity(pa, db, db)
  ds = wa / ws # degree of saturation
  sa = entropy(pa, db, pv)
  return (pa, db, wb, ha, wa, rh, td, pv, ps, va, ro, cp, za, ws, ds, sa)

def psy_dbwa(pa, db, wa):
  """Calculate all props from (pa,db,wa).
  """
  ps = psat(db)
  pv = pvap(pa, wa)
  rh = pv / ps
  td = tsat_inv(pv)
  wb = wetbulb(pa, db, wa, wb=(db+td)/2) # <<< non linear function
  ha = enthalpy(db, wa)
  va = specvol(pa, db, wa)
  ro = trudens(va, wa)
  cp = cpmoist(wa)
  za = sigheat(db, wb, wa)
  return (pa, db, wb, ha, wa, rh, td, pv, ps, va, ro, cp, za)

def calcstate(pa, db0=23, wb0=16, **kwargs):
  """Return primary props for pa and any other pair.
     air = calcstate(101.325, ha=45, wa=8.5e-3)

     *** Original Newton 2D method. ***
  """
  def efunc(dbx, wbx, arg1, arg2):
    """2D root finder error function.
       Returns a tuple of two errors.
    """
    ax = astate(pa, dbx, wbx)
    i, pi = arg1
    j, pj = arg2
    return sfact[i]*(pi-ax[i]), sfact[j]*(pj-ax[j])

  args = [(airprops.index(key),value) for key,value in kwargs.items()]
  db, wb = dnewton2(efunc, (db0, wb0), args)
  return astate(pa, db, wb)

# ---------------------------------------------------------------------
# EXPERIMENTAL
#   Test with different root finders methods     
#     - Damped Newton   <<<-- method of choice
#     - Trust Region Newton
#     - Broyden
#     - Secant  xxx not suitable for psy problems 
# ---------------------------------------------------------------------

def calcstate2(pa, db0=23, wb0=16, **kwargs):
  """Return primary props for pa and any other pair.
     air = calcstate2(101.325, ha=45, wa=8.5e-3)
  """
  def efunc(dbx, wbx, arg1, arg2):
    """2D root finder error function.
       Returns a tuple of two errors.
    """
    ax = astate(pa, dbx, wbx)
    i, pi = arg1
    j, pj = arg2
    return sfact[i]*(pi-ax[i]), sfact[j]*(pj-ax[j])

  args = [(airprops.index(key),value) for key,value in kwargs.items()]

  # With simple combinations and near guesses, 
  # all of these methods converge in 3 its...
  # I didn't do a really thorough test but dnewton2
  # seems to converge quickly for extreme variations.
  db, wb = dnewton2(efunc, (db0, wb0), args)
  #db, wb = broyden2(efunc, (db0, wb0), args)
  #db, wb = lm_newton2(efunc, (db0, wb0), args)

  return astate(pa, db, wb)

# ---------------------------------------------------------

def get_spec(pa, db, wb):
  air = TAir(*astate2(pa, db, wb))
  return {
      key: {**meta, "value": getattr(air, key) * meta["sf"]}
      for key, meta in air_spec.items()
  }

def print_spec(spec: dict, title: str = 'Moist Air Properties') -> None:
  """Prints a clean, aligned summary table from a spec dictionary."""
  print(f'\n{title}')
  print('-' * 52)
  for key, meta in spec.items():
    val_str = f"{meta['value']:.{meta['ndigs']}f}"
    unit_str = '' if meta['units'] == '-' else meta['units']
    print(f"{key}: {meta['label']:<24} {val_str:>12} {unit_str}")
  print('-' * 52)

# ---------------------------------------------------------------------

# Some spot tests during development.

def test1(): 
  """Test the astate tuple, namedtuple and TAir
  """
  pa, db, wb = 101.325, 23.0, 16.0
  air = astate(pa, db, wb) # air is tuple
  air = Air(*astate(pa, db, wb)) # air is namedtuple
  air = TAir(*astate(pa, db, wb)) # air is dataclass
  print(f"ha = {air.ha:0.3f} kJ/kg")

def test2():
  """psat and tsat comparison.
  """
  t = 20.0
  p = psat(t)
  print(f"psat = {p:0.6g} kPa")
  print(f"t={t:0.5g}, tsat={tsat(p):0.5g}, tdew={tdew(p):0.5g}, tsec={tsat_inv(p):0.9g}")

def test3():
  """Test the astate with TAir.
  """
  pa, db, wb = 101.325, 23.0, 16.0
  air = TAir(*astate(pa, db, wb)) # air is dataclass
  air.report()

def test4():
  """Test the db_sat function.
  """
  pa = 101.325
  db = 30.0
  ws = humidity(pa, db, db)
  db = db_sat(pa, ws)
  print(f"db={db:0.6g} degC, ws={1e3*ws:0.6g} g/kg")

def test5():
  """Test the get_spec and print_spec functions.
  """
  a_spec = get_spec(101.325, 23, 16)
  print_spec(a_spec, title="Moist Air State")

def test6():
  """Test the calcstate2 function.
     In calcstate2, change the 2D root finder method.
  """
  air = TAir(*calcstate2(101.325, ha=45, wa=8.5e-3))
  air.report()

def test7():
  air = astate(101.325, 23.0, 16.0)
  print(air)


# ---------------------------------------------------------------------

if __name__ == "__main__":

  # choose some air db/wb
  pa, db, wb = 101.325, 23, 16
  
  # use astate(pa, db, wb)
  air1 = TAir(*astate(pa, db, wb))
  air1.report()
  print()

  # Now, take the (ha, wa) from air1
  # and use it to calculate the state
  # with calcstate(pa, ha= , wa= )
  air2 = TAir(*calcstate(pa, ha=air1.ha, wa=air1.wa))
  air2.report()
  print()

  # Compare the reports
  if air1.is_close(air2):
    print('Success: Exactly the same ✔')
  else:
    print('Failed: Different! ✘')
