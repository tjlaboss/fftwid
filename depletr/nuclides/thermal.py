# Thermal Spectrum (PWR)

from .nuclide import Nuclide
from .half_lives import *


ALL_NUCLIDES = []

# Uranium
u235 = Nuclide('U', 235)
u235.sigma_f = 388
u235.sigma_y = 8.7
ALL_NUCLIDES.append(u235)
u238 = Nuclide('U', 238)
u238.sigma_f = 0.103
u238.sigma_y = 0.86
ALL_NUCLIDES.append(u238)
u239 = Nuclide('u', 239)
ALL_NUCLIDES.append(u239)
# Neptunium
np237 = Nuclide('Np', 237)
np237.sigma_f = 0.52
np237.sigma_y = 33
ALL_NUCLIDES.append(np237)
np238 = Nuclide('Np', 238)
np238.sigma_f = 134
np238.sigma_y = 13.6
ALL_NUCLIDES.append(np238)
np239 = Nuclide('Np', 239)
ALL_NUCLIDES.append(np239)
# Plutonium
pu238 = Nuclide('Pu', 238)
pu238.sigma_f = 2.4
pu238.sigma_y = 27.7
ALL_NUCLIDES.append(pu238)
pu239 = Nuclide('Pu', 239)
pu239.sigma_f = 102
pu239.sigma_y = 58.7
ALL_NUCLIDES.append(pu239)
pu240 = Nuclide('Pu', 240)
pu240.sigma_f = 0.53
pu240.sigma_y = 210.2
ALL_NUCLIDES.append(pu240)
pu241 = Nuclide('Pu', 241)
pu241.sigma_f = 102.2
pu241.sigma_y = 40.9
ALL_NUCLIDES.append(pu241)
pu242 = Nuclide('Pu', 242)
pu242.sigma_f = 0.44
pu242.sigma_y = 28.8
ALL_NUCLIDES.append(pu242)
# Americium
am241 = Nuclide('Am', 241)
am241.sigma_f = 1.1
am241.sigma_y = 110
ALL_NUCLIDES.append(am241)
am242 = Nuclide('Am', 242)
am242.sigma_f = 159
am242.sigma_y = 301
ALL_NUCLIDES.append(am242)
am242m = Nuclide('Am', 242)
am242m.name += 'm'  # metastable
am242m.latex += '$_m$'
am242m.sigma_f = 595
am242m.sigma_y = 137
ALL_NUCLIDES.append(am242m)
am243 = Nuclide('Am', 243)
am243.sigma_f = 0.44
am243.sigma_y = 49
ALL_NUCLIDES.append(am243)
# Curium
cm242 = Nuclide('Cm', 242)
cm242.sigma_f = 1.14
cm242.sigma_y = 4.5
ALL_NUCLIDES.append(cm242)
cm243 = Nuclide('Cm', 243)
cm243.sigma_f = 88
cm243.sigma_y = 14
ALL_NUCLIDES.append(cm243)
cm244 = Nuclide('Cm', 244)
cm244.sigma_f = 1.0
cm244.sigma_y = 16
ALL_NUCLIDES.append(cm244)
cm245 = Nuclide('Cm', 245)
cm245.sigma_f = 116
cm245.sigma_y = 17
ALL_NUCLIDES.append(cm245)

for _nuclide in ALL_NUCLIDES:
	_n = _nuclide.name
	if _n in ALPHA:
		_nuclide.lambda_alpha = LN2/ALPHA[_n]
	if _n in BETAM:
		_nuclide.lambda_betam = LN2/BETAM[_n]
	if _n in BETAP:
		_nuclide.lambda_betap = LN2/BETAP[_n]
	if _n in GAMMA:
		_nuclide.lambda_gamma = LN2/GAMMA[_n]
