# Unit tests for nuclide absorptions and ecays


from assert_test import assert_test as _assert
import sys; sys.path.append('..')
from depletr.nuclides.thermal import np238, pu238, pu242, am242, am242m


def test_neutron_capture():
	test_result = '{}(n, y){}'.format(pu238.name, pu238.capture())
	true_result = "Pu238(n, y)Pu239"
	_assert("Neutron capture", test_result, true_result)


def test_alpha_decay():
	test_result = '{}(n, a){}'.format(pu242.name, pu242.decay_alpha())
	true_result = "Pu242(n, a)U238"
	_assert("Alpha decay", test_result, true_result)
	
	
def test_betap_decay():
	test_result = '{}(n, e+){}'.format(am242.name, am242.decay_betap())
	true_result = "Am242(n, e+)Pu242"
	_assert("Beta+ decay", test_result, true_result)


def test_betam_decay():
	test_result = '{}(n, e-){}'.format(np238.name, np238.decay_betam())
	true_result = "Np238(n, e-)Pu238"
	_assert("Beta- decay", test_result, true_result)
	

def test_gamma_decay():
	test_result = '{}(*, y){}'.format(am242m.name, am242m.decay_gamma())
	true_result = "Am242m(*, y)Am242"
	_assert("Internal conversion", test_result, true_result)
