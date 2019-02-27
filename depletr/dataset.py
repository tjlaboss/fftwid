# Data Set

from collections import OrderedDict
import scipy as sp
from nuclides import Nuclide


class DataSet:
	def __init__(self, name):
		self.name = name
		self._nuclides = OrderedDict()
		self._q = OrderedDict()  # quantities
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
		for nuc, q in zip(list_of_nuclides, list_of_quantities):
			self.add_nuclide(nuc, q)
		
	def build_matrix(self):
		
		n = len(self._nuclides)
		A = sp.zeros((n, n))  # sigma
		L = sp.zeros((n, n))  # lambda
		
		indices = dict(zip(self._nuclides.keys(), range(n)))
		
		warnstr = "{} daughter {} of nuclide {} not in data set."
		for i, (key, nuclide) in enumerate(self._nuclides.items()):
			# Populate the diagonal with the absorption xs
			A[i, i] = nuclide.sigma_a
			daughter = nuclide.capture()
			if daughter in self._nuclides:
				j = indices[daughter]
				A[j, i] = -nuclide.sigma_y
			elif nuclide.sigma_y:
				print(warnstr.format("capture", daughter, nuclide.name))
			
			# Populate the decay matrix with the lambdas
			L[i, i] = nuclide.lambda_total
			if nuclide.lambda_betam:
				daughter = nuclide.decay_betam()
				if daughter in self._nuclides:
					j = indices[daughter]
					L[j, i] = -nuclide.lambda_betam
				elif nuclide.lambda_betam:
					print(warnstr.format("Beta-", daughter, nuclide.name))
			if nuclide.lambda_betap:
				daughter = nuclide.decay_betap()
				if daughter in self._nuclides:
					j = indices[daughter]
					L[j, i] = -nuclide.lambda_betap
				elif nuclide.lambda_betap:
					print(warnstr.format("Beta+", daughter, nuclide.name))
			if nuclide.lambda_alpha:
				daughter = nuclide.decay_alpha()
				if daughter in self._nuclides:
					j = indices[daughter]
					L[j, i] = -nuclide.lambda_alpha
				# OK to ignore alpha decay products
			# TODO: enable Am-242m later
			'''
			if nuclide.lambda_gamma:
				# e.g. if it's Am-242m
				daughter = nuclide.decay_gamma()
				if daughter in self._nuclides:
					j = indices[daughter]
					L[j, i] = -nuclide.lambda_gamma
			'''
		
		
		self._built = True
		return A, L
		

if __name__ == "__main__":
	st = DataSet("test set")
	'''
	from nuclides.thermal import u238, pu239, pu240
	st.add_nuclides([u238, pu239, pu240], [.333]*3)
	'''
	from nuclides.thermal import ALL_NUCLIDES
	quants = sp.ones(len(ALL_NUCLIDES))
	st.add_nuclides(ALL_NUCLIDES, quants)
	
	matA, matL = st.build_matrix()
	print(matA)
	print(matL)
	
	

