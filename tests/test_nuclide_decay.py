# Unit tests for nuclide absorptions and ecays


import sys; sys.path.append('..')
from depletr.nuclides.thermal import np238, pu238, pu242, am242


def _assert(name, test_result, true_result):
	errstr = "{} test failed:\n{} != {}".format(name, test_result, true_result)
	assert test_result == true_result, errstr


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
	
