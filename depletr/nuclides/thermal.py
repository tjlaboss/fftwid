# Thermal Spectrum (PWR)

from .nuclide import Nuclide


# Uranium
u235 = Nuclide('U', 235)
u235.sigma_f = 388
u235.sigma_y = 8.7
u238 = Nuclide('U', 238)
u238.sigma_f = 0.103
u238.sigma_y = 0.86
# Neptunium
np237 = Nuclide('Np', 237)
np237.sigma_f = 0.52
np237.sigma_y = 33
np238 = Nuclide('Np', 238)
np238.sigma_f = 134
np238.sigma_y = 13.6
# Plutonium
pu238 = Nuclide('Pu', 238)
pu238.sigma_f = 2.4
pu238.sigma_y = 27.7
pu239 = Nuclide('Pu', 239)
pu239.sigma_f = 102
pu239.sigma_y = 58.7
pu240 = Nuclide('Pu', 240)
pu240.sigma_f = 0.53
pu240.sigma_y = 210.2
pu241 = Nuclide('Pu', 241)
pu241.sigma_f = 102.2
pu241.sigma_y = 40.9
pu242 = Nuclide('Pu', 242)
pu242.sigma_f = 0.44
pu242.sigma_y = 28.8
# Americium
am241 = Nuclide('Am', 241)
am241.sigma_f = 1.1
am241.sigma_y = 110
am242 = Nuclide('Am', 242)
am242.sigma_f = 159
am242.sigma_y = 301
am242m = Nuclide('Am', 242)
am242m.name += 'm'  # metastable
am242m.sigma_f = 595
am242m.sigma_y = 137
am243 = Nuclide('Am', 243)
am243.sigma_f = 0.44
am243.sigma_y = 49
# Curium
cm242 = Nuclide('Cm', 242)
cm242.sigma_f = 1.14
cm242.sigma_y = 4.5
cm243 = Nuclide('Cm', 243)
cm243.sigma_f = 88
cm243.sigma_y = 14
cm244 = Nuclide('Cm', 244)
cm244.sigma_f = 1.0
cm244.sigma_y = 16
cm245 = Nuclide('Cm', 245)
cm245.sigma_f = 116
cm245.sigma_y = 17
