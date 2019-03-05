import scipy
from . import elements
from . import constants


class Nuclide:
	"""A nuclide with several useful physical quantities

	Parameters:
	-----------
	element:        str; atomic symbol for the element, e.g. "U"
	a:              int; mass number of the nucleus,    e.g. 235

	Attributes:
	-----------
	name:           str;   element+A, e.g., "U235" for Uranium-235
	alpha:          float; maximum fractional energy loss per collision
	xi:             float; mean logarithmic energy decrement
	lambda:         float, s^-1; decay constant
	halflife:       float, s; decay half-life
	sigma_n:        float, barns; (n, n) micro xs
	sigma_y:        float, barns; (n, y) micro xs
	sigma_f:        float, barns; (n, f) micro xs
	"""
	def __init__(self, element, a):
		element = element.title()
		self.element = element
		self.z = elements.Z[element]
		self.a = a
		self.name = element + str(a)
		self.latex = "${{}}^{{{}}}${}".format(a, element)
		# Scattering physics
		self._alpha = ((a - 1)/(a + 1))**2
		if a == 1:
			self._xi = 1
		else:
			self._xi = 1 + (self.alpha*scipy.log(self.alpha))/(1 - self.alpha)
		# Cross sections
		self.sigma_n = 0  # scatter
		self.sigma_y = 0  # capture
		self.sigma_f = 0  # fission
		# Decay
		self.lambda_alpha = 0  # alpha decay
		self.lambda_betap = 0  # beta+ decay
		self.lambda_betam = 0  # beta- decay
		self.lambda_gamma = 0  # internal conversion
		self._lambda_total = None
	
	@property
	def alpha(self):
		return self._alpha
	
	@property
	def xi(self):
		return self._xi
	
	@property
	def sigma_a(self):
		return self.sigma_f + self.sigma_y
	
	@property
	def lambda_total(self):
		if self._lambda_total is None:
			return self.lambda_betam + self.lambda_betap + \
			       self.lambda_alpha + self.lambda_gamma
	
	@lambda_total.setter
	def lambda_total(self, l):
		self._lambda_total = l
	
	@property
	def halflife(self):
		return scipy.log(2)/self._lambda_total
	
	@halflife.setter
	def halflife(self, t12):
		self._lambda_total = scipy.log(2)/t12

	def capture(self):
		"""Get the daughter nuclide from a neutron capture"""
		return self.element + str(self.a + 1)
	
	def decay_betap(self):
		"""Get the daughter nuclide from a Beta+ decay"""
		ep = elements.SYMBOL[self.z - 1]
		return ep + str(self.a)
	
	def decay_betam(self):
		"""Get the daughter nuclide from a Beta- decay"""
		em = elements.SYMBOL[self.z + 1]
		return em + str(self.a)
	
	def decay_alpha(self):
		"""Get the daughter nuclide from an alpha decay"""
		e = elements.SYMBOL[self.z - 2]
		a = self.a - 4
		return e + str(a)
	
	def decay_gamma(self):
		"""Get the daughter nuclide from internal conversion or gamma decay"""
		if self.name[-1] == "m":
			return self.name[:-1]
	
	def fission(self):
		if self.sigma_f:
			return [constants.FISSION_PRODUCT]
	
	def get_all_daughters(self):
		"""Nuclides, lock up your daughters
		
		Returns:
		--------
		daughters:  list of str; daughter nuclides from all possible decays
		"""
		daughters = []
		if self.lambda_alpha:
			daughters.append(self.decay_alpha())
		if self.lambda_betam:
			daughters.append(self.decay_betam())
		if self.lambda_betap:
			daughters.append(self.decay_betap())
		if self.lambda_gamma:
			daughters.append(self.decay_gamma())
		return daughters
