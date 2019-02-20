# Spectral data

from . import Nuclide


class DataSet:
	def __init__(self, name):
		self.name = name
		self._nuclides = {}
		self._q = {}  # quantities
		self._built = False
		
	
	def add_nuclide(self, nuc, q):
		"""Add a nuclide to the set
		
		Parameters:
		-----------
		nuc:        Nuclide
		
		"""
		if self._built:
			raise ValueError("Cannot add nuclides; matrix already built.")
		key = nuc.name
		if key in self._nuclides:
			errstr = "Nuclide {} already exists.".format(key)
			raise ValueError(errstr)
		self._nuclides[key] = nuc
		self._q[key] = q
	
	
	def add_nuclides(self, list_of_nuclides, list_of_quantities):
		assert len(list_of_nuclides) == len(list_of_quantities)
		for i, (nuc, q) in zip(list_of_nuclides, list_of_quantities):
			self.add_nuclide(nuc, q)
		
	def build_matrix(self):
		# TODO: Add method to build the matrix.
		# Complain if nuclides reference others that don't exist
		# in their decay chains.
		#self._built = True
		return
		



# Thermal Spectrum (PWR)
u235_t = Nuclide('U', 235)
u235_t.sigma_f = 388
u235_t.sigma_y = 8.7
u238_t = Nuclide('U', 238)
u238_t.sigma_f = 0.103
u238_t.sigma_y = 0.86
pu239_t = Nuclide('Pu', 239)
pu239_t.sigma_f = 102
pu239_t.sigma_y = 58.7
