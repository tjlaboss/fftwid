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
			else:
				print(warnstr.format("capture", daughter, nuclide.name))
			
			# Populate the decay matrix with the lambdas
			L[i, i] = nuclide.lambda_total
			if nuclide.lambda_betam:
				daughter = nuclide.decay_betam()
				if daughter in self._nuclides:
					j = indices[daughter]
					L[j, i] = -nuclide.lambda_betam
				else:
					print(warnstr.format("Beta-", daughter, nuclide.name))
			if nuclide.lambda_betap:
				daughter = nuclide.decay_betap()
				if daughter in self._nuclides:
					j = indices[daughter]
					L[j, i] = -nuclide.lambda_betap
				else:
					print(warnstr.format("Beta-", daughter, nuclide.name))
			# Ignore alpha decay products
		
		
		self._built = True
		return A, L
		

if __name__ == "__main__":
	st = DataSet("test set")
	from nuclides.thermal import u238, pu239, pu240
	st.add_nuclides([u238, pu239, pu240], [.333]*3)
	matA, matL = st.build_matrix()
	print(matA)

