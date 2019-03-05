# Unit tests for populating the decay matrix


from assert_test import assert_test as _assert
import sys; sys.path.append('..')
import depletr
from depletr.nuclides.thermal import am242m, am242


def _get_l_matrix():
	decayset = depletr.DataSet()
	decayset.add_nuclide(am242m, 1)
	decayset.add_nuclide(am242, 0)
	a, l = decayset.build_matrix()
	return l


def test_l_matrix_dimensions():
	l = _get_l_matrix()
	test_result = l.shape
	true_result = (4, 4)
	_assert("L matrix shape", test_result, true_result)


def test_am242m_lambda_total():
	l = _get_l_matrix()
	test_result = l[0, 0]
	true_result = am242m.lambda_total
	_assert("Americium 242 metastable lambda total", test_result, true_result)


def test_am242m_lambda_gamma():
	l = _get_l_matrix()
	test_result = l[1, 0]
	true_result = -am242m.lambda_gamma
	_assert("Americium 242 metastable lambda gamma", test_result, true_result)


def test_lumped_actinide():
	l = _get_l_matrix()
	test_result = l[-2, 1]
	true_result = -(am242.lambda_betap + am242.lambda_betam)
	_assert("Americium 242 lumped actinide product", test_result, true_result)
