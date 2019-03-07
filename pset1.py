# 22.215, PSet01
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Write a depletion code.
# Produce plots of nuclide inventory and k-infinity during burnup,
# and nuclide inventories vs. time after burnup.

import numpy as np
import depletr

POWER = 35      # W/g
ENRICHMENT = 5  # wt% U235
BURNUP = 60     # MW-d/kg-HM
DECAY_TIMES = np.logspace(-1, 4, 6) # Years
SPECTRA = ("thermal", "fast")
# For each spectrum:  Set `PLOT` to 0 for no plots,
# 1 for subplots on the same figure, or >1 for separate figures.
PLOT = 1

for spectrum in SPECTRA:
	print()
	header = "{} Reactor".format(spectrum.title())
	header += "\n" + "-"*len(header)
	print(header)
	
	study = depletr.Depleter(POWER, ENRICHMENT, max_burnup=BURNUP, spectrum=spectrum)
	print(">Depleting to {} MW-d/kg-HM:".format(BURNUP))
	spent_fuel = study.deplete(nsteps=1000, plots=PLOT)
	print("\n>Decaying in spent fuel pool for {} years".format(int(DECAY_TIMES.max())))
	times_in_seconds = DECAY_TIMES*depletr.nuclides.half_lives.YEAR
	respository = study.decay(spent_fuel, nsteps=1000, times=times_in_seconds)
	nuclide_names = study.get_all_nuclide_names()
	depletr.printer.print_decay_results(nuclide_names, DECAY_TIMES, respository)
	if PLOT:
		study.show()
