# Depleter
#
# Class to deplete nuclides

import scipy as sp
from scipy.linalg import expm as matrexp
from scipy.constants import N_A
from . import fuel, nuclides
from .dataset import DataSet

SPECTRA = {"fast"   : nuclides.fast,
           "thermal": nuclides.thermal}


class Depleter:
	def __init__(self, power, enrichment, max_burnup, spectrum,
				 mass=1E8, fuel_type="U"):
		self.power = power
		self.enrichment = enrichment / 100
		self.max_burnup = max_burnup
		self.mass = mass
		self.fuel_type = fuel_type
		assert spectrum in SPECTRA, \
			"Spectrum must be one of: {}".format(tuple(SPECTRA.keys()))
		self.data = SPECTRA[spectrum]
		self._scale = mass*N_A/238*1E-24


	def deplete(self, nsteps):
		power_megawatt = self.power*self.mass*1E-6
		fission_rate = fuel.get_fission_rate(power_megawatt)
		quantities = fuel.give_me_fuel(self.fuel_type, self.enrichment,
		                               len(self.data.ALL_NUCLIDES))
		time = fuel.give_me_fire(self.power, self.max_burnup)
		dt = time/nsteps
		
		ds = DataSet()
		ds.add_nuclides(self.data.ALL_NUCLIDES, quantities)
		m, l = ds.build_matrix()
		fv = ds.build_fission_vector()
		
		fluxvals = sp.zeros(nsteps)
		enrichvals = sp.zeros(nsteps)
		concentrations = sp.zeros((ds.size, nsteps))
		enrichvals[0] = self.enrichment
		c0 = ds.get_initial_quantities()*self._scale
		concentrations[:, 0] = c0
		fission_xs = (c0*fv).sum()
		flux = fission_rate/fission_xs
		
		for k in range(nsteps):
			a = m*flux + l
			dn = matrexp(-a*dt)
			if k > 0:
				concentrations[:, k] = concentrations[:, k-1].dot(dn)
				enrichvals[k] = concentrations[0, k]/(concentrations[0, k] + concentrations[1, k])
			fission_xs = (concentrations[:, k]*fv).sum()
			flux = fission_rate/fission_xs*1E-24
			fluxvals[k] = flux
		
		mass_reduct_238 = (concentrations[1, -1] - c0[1])/c0[1]
		print("U235 mass change: {:5.2%}".format(mass_reduct_238))
		mass_reduct_235 = (concentrations[0, -1] - c0[0])/c0[0]
		print("U235 mass change: {:5.2%}".format(mass_reduct_235))
		est235 = fuel.at_to_wt_uranium(enrichvals[-1])
		print("Final enrichment estimate: {:5.2%}".format(est235))
