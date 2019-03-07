# Data Set

from collections import OrderedDict
import scipy as sp
from . import nuclides


class DataSet:
	def __init__(self):
		self._nuclides = OrderedDict()
		self._q = OrderedDict()  # quantities
		self._size = None
		self._built = False
	
	@property
	def size(self):
		return self._size
	
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
		for nuc, q in zip(list_of_nuclides, list_of_quantities):
			self.add_nuclide(nuc, q)
		
	def build_matrix(self):
		indices = dict(zip(self._nuclides.keys(), range(len(self._nuclides))))
		indices[nuclides.FISSION_PRODUCT] = -1
		indices[nuclides.DEADEND_ACTINIDE] = -2
		
		n = len(indices)  # including the lumped boys
		A = sp.zeros((n, n))  # sigma
		L = sp.zeros((n, n))  # lambda
		
		warnstr = "{} daughter {} of nuclide {} not in data set."
		for i, (key, nuclide) in enumerate(self._nuclides.items()):
			# Populate the diagonal with the absorption xs
			A[i, i] = nuclide.sigma_a
			# Include the capture daughters...
			for daughter, branch_ratio in nuclide.capture():
				if not branch_ratio:
					continue
				if daughter in self._nuclides:
					j = indices[daughter]
					A[i, j] = -nuclide.sigma_y*branch_ratio
				elif nuclide.sigma_y:
					print(warnstr.format("capture", daughter, nuclide.name))
					j = indices[nuclides.DEADEND_ACTINIDE]
					A[i, j] = -nuclide.sigma_y*branch_ratio
			# ...and the fission products.
			if nuclide.sigma_f:
				j = indices[nuclides.FISSION_PRODUCT]
				A[i, j] = -nuclide.sigma_f
			
			# Populate the decay matrix with the lambdas
			L[i, i] = nuclide.lambda_total
			if nuclide.lambda_betam:
				daughter = nuclide.decay_betam()
				if daughter in self._nuclides:
					j = indices[daughter]
				else:
					print(warnstr.format("Beta-", daughter, nuclide.name))
					j = indices[nuclides.DEADEND_ACTINIDE]
				L[i, j] += -nuclide.lambda_betam
			if nuclide.lambda_betap:
				daughter = nuclide.decay_betap()
				if daughter in self._nuclides:
					j = indices[daughter]
				else:
					print(warnstr.format("Beta+", daughter, nuclide.name))
					j = indices[nuclides.DEADEND_ACTINIDE]
				L[i, j] += -nuclide.lambda_betap
			if nuclide.lambda_alpha:
				daughter = nuclide.decay_alpha()
				if daughter in self._nuclides:
					j = indices[daughter]
				else:
					print(warnstr.format("alpha", daughter, nuclide.name))
					j = indices[nuclides.DEADEND_ACTINIDE]
				L[i, j] += -nuclide.lambda_alpha
			if nuclide.lambda_gamma:
				# e.g. if it's Am-242m
				daughter = nuclide.decay_gamma()
				if daughter in self._nuclides:
					j = indices[daughter]
				else:
					j = indices[nuclides.DEADEND_ACTINIDE]
				L[i, j] += -nuclide.lambda_gamma
		
		
		self._built = True
		self._size = n
		return A, L
	
	
	def get_initial_quantities(self):
		if self._built:
			if not nuclides.DEADEND_ACTINIDE in self._nuclides:
				self._q[nuclides.DEADEND_ACTINIDE] = 0
			if not nuclides.FISSION_PRODUCT in self._nuclides:
				self._q[nuclides.FISSION_PRODUCT] = 0
			vals = tuple(self._q.values())
			return sp.array(vals)
	
	
	def build_xs_vector(self, rxn):
		assert self._built, "You must finalize the depletion matrix first."
		vector = sp.zeros(self._size)
		for i, nuclide in enumerate(self._nuclides.values()):
			if rxn == "fission":
				xs = nuclide.sigma_f
			elif rxn == "nu-fission":
				xs = nuclide.nu_sigma_f
			elif rxn == "absorption":
				xs = nuclide.sigma_a
			elif rxn == "capture":
				xs = nuclide.sigma_y
			else:
				raise NotImplementedError(xs)
			vector[i] = xs
		# Leave deadend actinides and fission products at 0.
		return vector
	
	
	def build_fission_vector(self):
		return self.build_xs_vector("fission")
