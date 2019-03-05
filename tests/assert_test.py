# Assert Test
#
# Just a function to make unit testing easier


def assert_test(name, test_result, true_result):
	errstr = "{} test failed:\n{} != {}".format(name, test_result, true_result)
	assert test_result == true_result, errstr
