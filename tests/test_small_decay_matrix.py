# Unit tests for populating the fission and decay matrices


import sys; sys.path.append('..')
import depletr
from depletr.nuclides.thermal import am242m, am242


def _assert(name, test_result, true_result):
	errstr = "{} test failed:\n{} != {}".format(name, test_result, true_result)
	assert test_result == true_result, errstr


def _get_l_matrix():
	decayset = depletr.DataSet()
	decayset.add_nuclide(am242m, 1)
	decayset.add_nuclide(am242, 0)
	a, l = decayset.build_matrix()
	return l


def test_l_matrix_dimensions():
	l = _get_l_matrix()
	_assert("L matrix shape", l.shape, (4, 4))


def test_am242m_lambda_total():
	l = _get_l_matrix()
	_assert("Americium 242 metastable lambda total",
	        l[0, 0], am242m.lambda_total)


def test_am242m_lambda_gamma():
	l = _get_l_matrix()
	_assert("Americium 242 metastable lambda gamma",
	        l[1, 0], -am242m.lambda_gamma)


def test_lumped_actinide():
	l = _get_l_matrix()
	print(l)
	beta_lam = am242.lambda_betap + am242.lambda_betam
	_assert("Americium 242 lumped actinide product (beta decay)",
	        l[-2, 1], -beta_lam)

