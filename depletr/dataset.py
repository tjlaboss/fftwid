# Data Set

from collections import OrderedDict
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
		for i, (nuc, q) in zip(list_of_nuclides, list_of_quantities):
			self.add_nuclide(nuc, q)
		
	def build_matrix(self):
		# TODO: Add method to build the matrix.
		# Complain if nuclides reference others that don't exist
		# in their decay chains.
		#self._built = True
		return
		



