# Depleter
#
# Class to deplete nuclides

import matplotlib.pyplot as plt
import scipy as sp
from scipy.linalg import expm as matrexp
from scipy.constants import N_A
from . import fuel, nuclides
from . import plotter
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
	
	
	def get_all_nuclides(self):
		return self.data.ALL_NUCLIDES


	def deplete(self, nsteps, plots=1):
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
		absv = ds.build_xs_vector(rxn="absorption")
		nufv = ds.build_xs_vector(rxn="nu-fission")
		
		kinfvals = sp.zeros(nsteps)
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
			ck = concentrations[:, k]
			fission_xs = (ck*fv).sum()
			flux = fission_rate/fission_xs*1E-24
			fluxvals[k] = flux
			kinfvals[k] = (ck*nufv).sum()/(ck*absv).sum()
		
		mass_reduct_238 = (concentrations[1, -1] - c0[1])/c0[1]
		print("U238 mass change: {:5.2%}".format(mass_reduct_238))
		mass_reduct_235 = (concentrations[0, -1] - c0[0])/c0[0]
		print("U235 mass change: {:5.2%}".format(mass_reduct_235))
		est235 = fuel.at_to_wt_uranium(enrichvals[-1])
		print("Final enrichment estimate: {:5.2%}".format(est235))
		
		if plots:
			tvals = sp.arange(0, nsteps*dt, dt)
			if plots == 1:
				fig = plt.figure()
				axa = fig.add_subplot(221)
				axf = fig.add_subplot(222)
				axk = fig.add_subplot(223)
				axu = fig.add_subplot(224)
			else:
				axa = plt.figure().add_subplot(111)
				axf = plt.figure().add_subplot(111)
				axk = plt.figure().add_subplot(111)
				axu = plt.figure().add_subplot(111)
			# Top left: Actinide depletion
			plotter.make_actinides_plot(tvals, concentrations, self.data.ALL_NUCLIDES, axa,
			                            fission_products=True, deadend_actinides=True)
			# Top Right: Enrichment and flux
			plotter.make_enrichment_flux_plot(tvals, enrichvals, fluxvals, axf)
			# Bottom left: Approximate k-infinity
			plotter.make_kinf_plot(tvals, kinfvals, axk)
			# Bottom Right: Relative uranium depletion
			cvals = concentrations[0:2]
			nucs = (self.data.u235, self.data.u238)
			plotter.make_element_depletion_plot(tvals, cvals, nucs, axu,
			                                    relative=True, element="Uranium")
			if plots == 1:
				# Hide redundant xlabels to avoi crowding
				axa.set_xlabel("")
				axf.set_xlabel("")
		
		return concentrations[:, -1]
	
	def show(self):
		return plt.show()
