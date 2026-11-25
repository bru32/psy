# psy_broyden_demo.py
#
# Write ANY psychrometic solver using a 2D root finder.
#
# 1. Create a function named ab_from_cd(c, d)
# 2. Write the error function that calculates 
#    (cx, dx) for the current (ax,bx).
#    return the error (cx-c) and (dx-d).
# 3. One small tweak is to apply the scale factor.
#
# Author: Bruce Wernick
# Date: 25 November 2026
#
# Copyright © (2026) Bruce Wernick,  All rights reserved...
# ---------------------------------------------------------

from nu_meth import broyden_nd
from psy_wernick import *

# ---------------------------------------------------------

# plug in newton 2D solver for testing

def jac(f, x, y, f1, f2, h=1e-6):
  f1_x, f2_x = f((x + h, y))
  f1_y, f2_y = f((x, y + h))
  j11 = (f1_x - f1) / h
  j21 = (f2_x - f2) / h
  j12 = (f1_y - f1) / h
  j22 = (f2_y - f2) / h
  det = j11 * j22 - j12 * j21
  return j11, j21, j12, j22, det

def newton2d(f, guess, tol=1e-9, maxits=48):
  x, y = guess
  for its in range(maxits):
    f1, f2 = f((x, y))
    if abs(f1) < tol and abs(f2) < tol:
      return x, y
    j11, j21, j12, j22, det = jac(f, x, y, f1, f2)
    if abs(det) < 1e-15:
      raise ZeroDivisionError("Singular Jacobian encountered")
    x -= (j22 * f1 - j12 * f2) / det
    y -= (j11 * f2 - j21 * f1) / det
  raise RuntimeError(f"Did not converge within {maxits} iterations")

# ---------------------------------------------------------

pa = 101.325
scale = {"ha":45, "wa":8.5e-3, "rh": 0.5, "pv": 1.2}

def dbwb_from_hawa(ha, wa):
  """Calculate dry bulb and wet bulb temperature
     from specific enthalpy and abs. humidity.
  """

  def f(x):
    """Error function. 
       System is solved when f(x)=0.
    """
    dbx, wbx = x
    wax = wa_from_dbwb(pa, dbx, wbx)
    hax = ha_from_dbwa(dbx, wax)
    return (hax - ha) / scale["ha"], (wax - wa) / scale["wa"]

  guess = (23, 16)
  #x, its, ok = broyden_nd(f, guess)
  x = newton2d(f, guess)

  return x

def dbwb_from_rhpv(rh, pv):
  """Calculate dry bulb and wet bulb temperature
     from relative humidity and vapor pressure.
  """
  def f(x):
    dbx, wbx = x
    psx = psat(dbx)
    wax = wa_from_dbwb(pa, dbx, wbx)
    pvx = pv_from_wa(pa, wax)
    rhx = pvx / psx
    return (rhx - rh) / scale["rh"], (pvx - pv) / scale["pv"]
  
  guess = (23, 16)
  #x, its, ok = broyden_nd(f, guess)
  x = newton2d(f, guess)

  return x

# Example 1. Calculate db/wb from (ha, wa)
db, wb = dbwb_from_hawa(45.678, 8.7654e-3)
print(f"Result: db = {db:.3f} °C, wb = {wb:.3f} °C")

# Example 2. Calculate db/wb from (rh, pv)
db, wb = dbwb_from_rhpv(0.45, 1.23)
print(f"Result: db = {db:.3f} °C, wb = {wb:.3f} °C")
