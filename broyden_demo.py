# broyden_demo.py

# How to write ANY psychrometic solver
# Create a function named ab_from_cd(c, d)
# Write the error function that calculates 
# (cx, dx) for the current (ax,bx).
# return the error (cx-c) and (dx-d).
# One small tweak is to apply the scale factor.


from nu_meth import broyden_nd
from psy_wernick import *

pa = 101.325
scale = {"ha":45, "wa":8.5e-3}

def dbwb_from_hawa(ha, wa):
  """Calculate (db, wb) from (ha, wa).
  """

  def func(x):
    """error function"""
    dbx, wbx = x
    wax = wa_from_dbwb(pa, dbx, wbx)
    hax = ha_from_dbwa(dbx, wax)
    return (hax - ha) / scale["ha"], (wax - wa) / scale["wa"]

  guess = (23, 16)
  x, its, ok = broyden_nd(func, guess)
  return x

db, wb = dbwb_from_hawa(45.678, 8.7654e-3)
print(f"Result: db = {db:.3f} °C, wb = {wb:.3f} °C")
