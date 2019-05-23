# Unit tests for exponential decay

from assert_test import assert_test as _assert
from assert_test import assert_test_array as _assert_array
import numpy as np
from scipy.linalg import expm
import sys; sys.path.append('..')
import depletr
from depletr.nuclides.thermal import u235, u238, u239, np239, pu238, pu239, pu242

_ND = 4  # Number of decimal places to round matrix exponential calcs to


def test_half_life():
	dt = depletr.nuclides.half_lives.BETAM[u239.name]
	test_result = np.exp(-u239.lambda_total*dt)
	true_result = 0.5
	_assert("U239 half-life", test_result, true_result)


def test_matrix_half_life():
	dt = depletr.nuclides.half_lives.ALPHA[pu242.name]
	ds = depletr.DataSet()
	ds.add_nuclide(u238, 0)
	ds.add_nuclide(pu242, 1)
	ds.build()
	q0 = ds.get_initial_quantities()
	q1 = q0.dot(expm(-ds.l*dt))
	test_result = (q1[0], q1[1])
	true_result = (0.5,     0.5)
	_assert("Pu242 matrix half-life", test_result, true_result)


def _get_depleter_decay():
	dt = depletr.nuclides.half_lives.ALPHA[pu239.name]
	depper = depletr.Depleter(power=1, enrichment=0, spectrum="thermal", max_burnup=0)
	nuclides_names = depper.get_all_nuclide_names()
	nuclides_names += [depletr.nuclides.constants.DEADEND_ACTINIDE,
	                   depletr.nuclides.constants.FISSION_PRODUCT]
	num = len(nuclides_names)
	indices = dict(zip(nuclides_names, range(num)))
	q0 = np.zeros(num)
	i_np239 = indices[np239.name]
	q0[i_np239] = 1
	i_pu238 = indices[pu238.name]
	q0[i_pu238] = 1
	qhist = depper.decay(q0, 1, [dt, 2*dt])
	return indices, qhist


def test_depleter_conservation():
	qhist = _get_depleter_decay()[1]
	test_result = np.round(qhist.sum(axis=0), _ND)
	true_result = 2*np.ones(len(test_result))
	_assert_array("Conservation", test_result, true_result)


def test_depleter_np239():
	indices, qhist = _get_depleter_decay()
	i_np239 = indices[np239.name]
	test_result = np.round(qhist[i_np239], _ND)
	true_result = np.array([1, 0, 0], dtype=float)
	_assert_array("Np239 history", test_result, true_result)


def test_depleter_pu239():
	indices, qhist = _get_depleter_decay()
	i_pu239 = indices[pu239.name]
	test_result = np.round(qhist[i_pu239], _ND)
	true_result = np.array([0, 0.5, 0.25], dtype=float)
	_assert_array("Pu239 history", test_result, true_result)


def test_depleter_u235():
	indices, qhist = _get_depleter_decay()
	i_u235 = indices[u235.name]
	test_result = np.round(qhist[i_u235], _ND)
	true_result = np.array([0, 0.5, 0.75], dtype=float)
	_assert_array("U235 history", test_result, true_result)
	


