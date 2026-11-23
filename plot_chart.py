
from psy_wernick import *

# Sea level 101.325 kPa
# 700 m above sea level 93.19 kPa
# 1400 m above sea level 85.6 kPa
# 1700 m above sea level 82.5 kPa
# low temp (-40 to 10 °C) (0 to 2 g/kg)
# med temp (-10 to 45 °C) (0 to 30 g/kg)
# high temp (10 to 90 °C) (0 to 100 g/kg)
# extreme temperature chart (-60 to 90 °C) (0 to 10 g/kg)

# plot chart
plotchart(pa=101.325, db_range=(5, 40), wa_range=(0, 0.030))
