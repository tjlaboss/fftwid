# Assert Test
#
# Just a function to make unit testing easier


def assert_test(name, test_result, true_result):
	errstr = "{} test failed:\n{} != {}".format(name, test_result, true_result)
	assert test_result == true_result, errstr


def assert_test_array(name, test_array, true_array):
	errstr = "{} test failed:\n{} != {}".format(name, test_array, true_array)
	assert (test_array == true_array).all(), errstr
