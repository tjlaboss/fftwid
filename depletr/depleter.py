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
		self._last_dataset = None
	
	
	def get_all_nuclides(self):
		return self.data.ALL_NUCLIDES
	
	
	def get_all_nuclide_names(self):
		return [n.name for n in self.data.ALL_NUCLIDES]
	
	
	def _deplete(self, nsteps, dt, ds, c0, fission_rate):
		m = ds.m
		l = ds.l
		kinfvals = sp.zeros(nsteps)
		fluxvals = sp.zeros(nsteps)
		enrichvals = sp.zeros(nsteps)
		concentrations = sp.zeros((ds.size, nsteps))
		enrichvals[0] = self.enrichment
		
		concentrations[:, 0] = c0
		fv = ds.get_fission_vector()
		nufv = ds.get_xs_vector("nu-fission")
		absv = ds.get_xs_vector("absorption")
		fission_xs = (c0*fv).sum()
		flux = fission_rate/fission_xs
		
		for k in range(nsteps):
			if k > 0:
				a = m*flux + l
				dn = matrexp(-a*dt)
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
		return concentrations, kinfvals, fluxvals, enrichvals

	
	def deplete_fresh(self, nsteps, plots=1):
		power_megawatt = self.power*self.mass*1E-6
		fission_rate = fuel.get_fission_rate(power_megawatt)
		quantities = fuel.give_me_fuel(self.fuel_type, self.enrichment,
		                               len(self.data.ALL_NUCLIDES))
		time = fuel.give_me_fire(self.power, self.max_burnup)
		dt = time/nsteps
		
		ds = DataSet()
		ds.add_nuclides(self.data.ALL_NUCLIDES, quantities)
		ds.build()
		c0 = ds.get_initial_quantities()*self._scale
		
		concs, kinf, flux, enrich = self._deplete(nsteps, dt, ds, c0, fission_rate)
		
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
			plotter.make_actinides_plot(tvals, concs, self.data.ALL_NUCLIDES, axa,
			                            fission_products=True, deadend_actinides=True)
			# Top Right: Enrichment and flux
			plotter.make_enrichment_flux_plot(tvals, enrich, flux, axf)
			# Bottom left: Approximate k-infinity
			plotter.make_kinf_plot(tvals, kinf, axk)
			# Bottom Right: Relative uranium depletion
			cvals = concs[0:2]
			nucs = (self.data.u235, self.data.u238)
			plotter.make_element_depletion_plot(tvals, cvals, nucs, axu,
			                                    relative=True, element="Uranium")
			if plots == 1:
				# Hide redundant xlabels to avoid crowding
				axa.set_xlabel("")
				axf.set_xlabel("")
		
		self._last_dataset = ds
		return concs[:, -1]
	
	
	def decay(self, quantities, nsteps, times):
		ds = DataSet()
		ds.add_nuclides(self.data.ALL_NUCLIDES, quantities[:-2])
		ds.build()
		# Ignore fission products
		l = ds.l[:-1, :-1]
		c0 = quantities[:-1]
		# Deplete over the intervals of interest with nsteps each
		nt = len(times)
		nc = len(c0)
		results = sp.zeros((nc, nt + 1))
		results[:, 0] = c0
		concentrations = sp.zeros((nc, nt*nsteps + 1))
		concentrations[:, 0] = c0
		elapsed = 0
		i = 0
		for interval, time in enumerate(times):
			dt = (time - elapsed)/nsteps
			for k in range(nsteps):
				elapsed += dt
				dn = matrexp(-l*dt)
				concentrations[:, i+1] = concentrations[:, i].dot(dn)
			results[:, interval + 1] = concentrations[:, i+1]
			i += 1
		
		return results
	
	
	def show(self):
		return plt.show()
