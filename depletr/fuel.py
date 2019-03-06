# Fuel


import scipy
from scipy.constants import eV
from nuclides.half_lives import DAY


ELEMENTS = ("U", "Pu", "Uranium", "Plutonium")


def give_me_fuel(element, enrichment, num_nuclides):
	assert element in ELEMENTS, "Implemented fuel types: {}".format(ELEMENTS)
	assert num_nuclides >= 2, "You cannot have fewer than 2 nuclides."
	fuel_nuclides = scipy.zeros(num_nuclides)
	if element[0] == "U":
		at235 = wt_to_at_uranium(enrichment)
		at238 = 1 - at235
		fuel_nuclides[0] = at235
		fuel_nuclides[1] = at238
	else:
		# Plutonium, mox?
		raise NotImplementedError(element)
	return fuel_nuclides


def give_me_fire(watts_per_gram, burnup):
	"""Get the time to achieve a required burnup"""
	burnup *= DAY*1E3  # convert from MWd/kg to W-s/g
	return burnup/watts_per_gram


def get_fission_rate(megawatts, mev_per_fission=200):
	"""Get the number of required fissions"""
	mevs = megawatts/eV
	return mevs/mev_per_fission


def wt_to_at_uranium(wt235):
	assert 0 <= wt235 <= 1, "Weight fraction must be on [0, 1]."
	wt238 = 1 - wt235
	at238 = wt238/238
	at235 = wt235/235
	return at235/(at235 + at238)


def at_to_wt_uranium(at235):
	assert 0 <= at235 <= 1, "Atom fraction must be on [0, 1]."
	at238 = 1 - at235
	wt238 = at238*238
	wt235 = at235*235
	return wt235/(wt235 + wt238)
	

def wt_to_at_plutonium(wt239):
	return NotImplementedError("Plutonium fuel")
