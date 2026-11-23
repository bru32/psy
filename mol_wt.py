# Molecular Weight of Air and Water

# Standard Air
gases = [
    ("Nitrogen (N2)", 0.78084, 28.0134),
    ("Oxygen (O2)", 0.209476, 31.9988),
    ("Argon (Ar)", 0.00934, 39.948),
    ("Carbon Dioxide (CO2)", 0.000314, 44.0095),
    ("Neon (Ne)", 0.00001818, 20.1797),
    ("Helium (He)", 0.00000524, 4.0026),
    ("Methane (CH4)", 0.000002, 16.0425),
    ("Krypton (Kr)", 0.00000114, 83.798),
    ("Hydrogen (H2)", 0.0000005, 2.0159),
    ("Xenon (Xe)", 0.000000087, 131.293)
]

total_frac = sum(g[1] for g in gases)
total_mw = sum(g[1] * g[2] for g in gases)
MA = total_mw / total_frac
print("Standard Air")
print("Normalized MA:", MA)


# With modern CO2 (~350-400 ppm, adjusted O2)
gases_mod = [
    ("Nitrogen (N2)", 0.78084, 28.0134),
    ("Oxygen (O2)", 0.20946, 31.9988),
    ("Argon (Ar)", 0.00934, 39.948),
    ("Carbon Dioxide (CO2)", 0.00033, 44.0095),
    ("Neon (Ne)", 0.00001818, 20.1797),
    ("Helium (He)", 0.00000524, 4.0026),
    ("Methane (CH4)", 0.00000179, 16.0425),
    ("Krypton (Kr)", 0.00000114, 83.798),
    ("Hydrogen (H2)", 0.0000005, 2.0159),
    ("Xenon (Xe)", 0.000000087, 131.293)
]

print("Modified Air (with modern CO2 levels)")
print("Total frac mod:", sum(g[1] for g in gases_mod))
print("MA mod:", sum(g[1] * g[2] for g in gases_mod))
print("Normalized MA mod:", sum(g[1] * g[2] for g in gases_mod) / sum(g[1] for g in gases_mod))


# Water
h = 1.00794
o = 15.9994
MW = 2 * h + o
print(f"Mol Weight of water: MW = {MW:.5f}")

R = MW / MA 
print(f"{R}")