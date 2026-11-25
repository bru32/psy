from mpmath import mp
from psy_wernick import *

mp.dps = 9

pa = 101.325
scale = {"ha": 45, "wa": 8.5e-3}

def dbwb_from_hawa(ha, wa):

  def f(dbx, wbx):
    wax = wa_from_dbwb(pa, dbx, wbx)
    hax = ha_from_dbwa(dbx, wax)
    return (hax - ha) / scale["ha"], (wax - wa) / scale["wa"]

  guess = (23, 16)
  db, wb = mp.findroot(f, guess, solver='mdnewton')

  return db, wb

db, wb = dbwb_from_hawa(45.678, 8.7654e-3)
print(f"Result: db = {db:.3f} °C, wb = {wb:.3f} °C")
