# broyden_demo.py
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

pa = 101.325
scale = {"ha":45, "wa":8.5e-3}

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
  x, its, ok = broyden_nd(f, guess)
  return x

db, wb = dbwb_from_hawa(45.678, 8.7654e-3)
print(f"Result: db = {db:.3f} °C, wb = {wb:.3f} °C")
