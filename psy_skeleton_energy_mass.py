# VERY simplified script to demonstrate how to 
# plot psy chart that is linear in energy and mass.

from math import exp
import tkinter as tk

cpa = 1.006
cpv = 1.86
hgo = 2501

def psat(db):
  """Magnus"""
  return 0.61078 * exp((17.27 * db) / (db + 237.3))

def humidity(pa, db, rh):
  ps = psat(db)
  pv = rh * ps
  if pv >= pa:
    return 0.050 # max humidity
  return 0.621945 * (pv / (pa - pv))

def enthalpy(db, wa):
  return cpa * db + wa * (hgo + cpv * db)

def toScreen(ha, wa):
  W = width - 2*pad
  H = height - 2*pad
  xs = ha - wa * sf # x-sheared
  x = pad + W * (xs - xmin) / (xmax - xmin)
  y = (height - pad) - H * (wa - wmin) / (wmax - wmin)
  return x, y

def toChart(x, y):
  W = width - 2*pad
  H = height - 2*pad
  wa = wmin + ((height - pad - y) / H) * (wmax - wmin)
  xs = xmin + ((x - pad) / W) * (xmax - xmin)
  ha = xs + wa * sf
  db = (ha - wa * hgo) / (cpa + wa * cpv)
  return db, wa

def draw_wa_lines():
  for i in range(0, 31, 5):
    wa = 1e-3 * i
    x0, y0 = pad, toScreen(0, wa)[1]
    x1, y1 = width - pad, y0
    canvas.create_line(x0, y0, x1, y1)

def draw_db_lines():
  for db in range(0, 51, 5):
    wa0 = 0.0
    ha0 = enthalpy(db, wa0)
    wa1 = min(humidity(pa, db, 1.0), wmax)
    ha1 = enthalpy(db, wa1)
    x0, y0 = toScreen(ha0, wa0)
    x1, y1 = toScreen(ha1, wa1)
    width = 2.0 if db == int(tref) else 1.2
    canvas.create_line(x0, y0, x1, y1, width=width, fill="maroon")

def draw_rh_lines():
  for pc in range(20, 101, 20):
    rh = pc / 100.0
    points = []
    for t in range(0, 505, 5):
      db = t / 10.0
      wa = humidity(pa, db, rh)
      ha = enthalpy(db, wa)
      xs = ha - wa * sf # x-sheared
      if (xmin <= xs <= xmax) and (wmin <= wa <= wmax):
        points.extend(toScreen(ha, wa))
    if len(points) >= 4:
      width = 2.5 if pc == 100 else 1.0
      canvas.create_line(*points, width=width, smooth=True)

def move_mouse(event):
  canvas.delete("label")
  x, y = event.x, event.y
  db, wa = toChart(x, y)
  canvas.create_text(pad+10, pad+10, 
    text=f"db {db:4.1f} degC\nwa {1e3*wa:4.1f} g/kg",
    font=("Consolas", 14),
    anchor="nw", tags="label")

# ---------------------------------------------------------------------

root = tk.Tk()
root.title("Energy Mass PsyChart")
width, height = 800, 500
root.geometry(f"{width}x{height}+20+20")
canvas = tk.Canvas(root)
canvas.pack(fill="both", expand=True)
canvas.bind("<Motion>", move_mouse)
pad = 40
pa = 101.325
tref = 50.0 # max dry bulb temp
sf = hgo + cpv * tref # shear factor
xmin, xmax = 0.0, cpa * tref
wmin, wmax = 0.0, 0.030
draw_wa_lines()
draw_db_lines()
draw_rh_lines()
root.mainloop()