# Fast Spectrum

from .nuclide import Nuclide


# Uranium
u235 = Nuclide('U', 235)
u235.sigma_f = 1.98
u235.sigma_y = 0.57
u238 = Nuclide('U', 238)
u238.sigma_f = 0.04
u238.sigma_y = 0.30
# Neptunium
np237 = Nuclide('Np', 237)
np237.sigma_f = 0.32
np237.sigma_y = 1.7
np238 = Nuclide('Np', 238)
np238.sigma_f = 3.6
np238.sigma_y = 0.2
# Plutonium
pu238 = Nuclide('Pu', 238)
pu238.sigma_f = 1.1
pu238.sigma_y = 0.58
pu239 = Nuclide('Pu', 239)
pu239.sigma_f = 1.86
pu239.sigma_y = 0.56
pu240 = Nuclide('Pu', 240)
pu240.sigma_f = 0.36
pu240.sigma_y = 0.57
pu241 = Nuclide('Pu', 241)
pu241.sigma_f = 2.49
pu241.sigma_y = 0.47
pu242 = Nuclide('Pu', 242)
pu242.sigma_f = 0.24
pu242.sigma_y = 0.44
# Americium
am241 = Nuclide('Am', 241)
am241.sigma_f = 0.27
am241.sigma_y = 2.0
am242 = Nuclide('Am', 242)
am242.sigma_f = 3.2
am242.sigma_y = 0.6
am242m = Nuclide('Am', 242)
am242m.name += 'm'  # metastable
am242m.sigma_f = 3.3
am242m.sigma_y = 0.6
am243 = Nuclide('Am', 243)
am243.sigma_f = 0.21
am243.sigma_y = 1.8
# Curium
cm242 = Nuclide('Cm', 242)
cm242.sigma_f = 0.58
cm242.sigma_y = 1.0
cm243 = Nuclide('Cm', 243)
cm243.sigma_f = 7.2
cm243.sigma_y = 1.0
cm244 = Nuclide('Cm', 244)
cm244.sigma_f = 0.42
cm244.sigma_y = 0.6
cm245 = Nuclide('Cm', 245)
cm245.sigma_f = 5.1
cm245.sigma_y = 0.9
