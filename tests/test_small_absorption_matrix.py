# Unit tests for populating the fission/capture (absorption) matrix


from assert_test import assert_test as _assert
import sys; sys.path.append('..')
import depletr
from depletr.nuclides.thermal import u238, u239


def _get_m_matrix():
	decayset = depletr.DataSet()
	decayset.add_nuclide(u238, 1)
	decayset.add_nuclide(u239, 0)
	decayset.build()
	return decayset.m


def test_a_matrix_dimensions():
	a = _get_m_matrix()
	test_result = a.shape
	true_result = (4, 4)
	_assert("A matrix shape", test_result, true_result)


def test_u239_capture():
	a = _get_m_matrix()
	test_result = a[0, 1]
	true_result = -u238.sigma_y
	_assert("Uranium 238 -> Uranium 239", test_result, true_result)


def test_u238_absorption():
	a = _get_m_matrix()
	test_result = a[0, 0]
	true_result = u238.sigma_f + u238.sigma_y
	_assert("Uranium 238 absorption", test_result, true_result)


def test_lumped_actinide():
	a = _get_m_matrix()
	test_result = a[-1, -2]
	true_result = -u239.sigma_y
	_assert("Uranium 239 capture", test_result, true_result)


def test_lumped_fission_product():
	a = _get_m_matrix()
	test_result = a[0, -1]
	true_result = -u238.sigma_f
	_assert("Uranium 238 fission", test_result, true_result)
