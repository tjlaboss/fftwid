# Fast Spectrum

from .nuclide import Nuclide


ALL_NUCLIDES = []

# Uranium
u235 = Nuclide('U', 235)
u235.sigma_f = 1.98
u235.sigma_y = 0.57
ALL_NUCLIDES.append(u235)
u238 = Nuclide('U', 238)
u238.sigma_f = 0.04
u238.sigma_y = 0.30
ALL_NUCLIDES.append(u238)
# Neptunium
np237 = Nuclide('Np', 237)
np237.sigma_f = 0.32
np237.sigma_y = 1.7
ALL_NUCLIDES.append(np237)
np238 = Nuclide('Np', 238)
np238.sigma_f = 3.6
np238.sigma_y = 0.2
ALL_NUCLIDES.append(np238)
# Plutonium
pu238 = Nuclide('Pu', 238)
pu238.sigma_f = 1.1
pu238.sigma_y = 0.58
ALL_NUCLIDES.append(pu238)
pu239 = Nuclide('Pu', 239)
pu239.sigma_f = 1.86
pu239.sigma_y = 0.56
ALL_NUCLIDES.append(pu239)
pu240 = Nuclide('Pu', 240)
pu240.sigma_f = 0.36
pu240.sigma_y = 0.57
ALL_NUCLIDES.append(pu240)
pu241 = Nuclide('Pu', 241)
pu241.sigma_f = 2.49
pu241.sigma_y = 0.47
ALL_NUCLIDES.append(pu241)
pu242 = Nuclide('Pu', 242)
pu242.sigma_f = 0.24
pu242.sigma_y = 0.44
ALL_NUCLIDES.append(pu242)
# Americium
am241 = Nuclide('Am', 241)
am241.sigma_f = 0.27
am241.sigma_y = 2.0
ALL_NUCLIDES.append(am241)
am242 = Nuclide('Am', 242)
am242.sigma_f = 3.2
am242.sigma_y = 0.6
ALL_NUCLIDES.append(am242)
am242m = Nuclide('Am', 242)
am242m.name += 'm'  # metastable
am242m.sigma_f = 3.3
am242m.sigma_y = 0.6
ALL_NUCLIDES.append(am242m)
am243 = Nuclide('Am', 243)
am243.sigma_f = 0.21
am243.sigma_y = 1.8
ALL_NUCLIDES.append(am243)
# Curium
cm242 = Nuclide('Cm', 242)
cm242.sigma_f = 0.58
cm242.sigma_y = 1.0
ALL_NUCLIDES.append(cm242)
cm243 = Nuclide('Cm', 243)
cm243.sigma_f = 7.2
cm243.sigma_y = 1.0
ALL_NUCLIDES.append(cm243)
cm244 = Nuclide('Cm', 244)
cm244.sigma_f = 0.42
cm244.sigma_y = 0.6
ALL_NUCLIDES.append(cm244)
cm245 = Nuclide('Cm', 245)
cm245.sigma_f = 5.1
cm245.sigma_y = 0.9
ALL_NUCLIDES.append(cm245)
