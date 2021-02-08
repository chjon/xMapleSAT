import subexpressions

class TestResult:
	def __init__(self, testResults = 0, totalNumTests = 0):
		self.testResults   = testResults
		self.totalNumTests = totalNumTests
	
	def __iadd__(self, other):
		self.testResults   += other.testResults
		self.totalNumTests += other.totalNumTests
		return self
	
	def result_str(self):
		return "{}/{} tests passed".format(self.totalNumTests - self.testResults, self.totalNumTests)

TEST_OKAY = TestResult(0, 1)
TEST_FAIL = TestResult(1, 1)

def expect(testResult, actual, expected):
	if actual != expected:
		print("Test {} failed: actual = {}; expected = {}".format(testResult.totalNumTests, actual, expected))
		testResult += TEST_FAIL
	else: testResult += TEST_OKAY
	
def test_intersect():
	testResult = TestResult(0, 0)

	expect(testResult, subexpressions.intersect([ ], [ ]), [ ])
	expect(testResult, subexpressions.intersect([1], [ ]), [ ])
	expect(testResult, subexpressions.intersect([ ], [1]), [ ])
	expect(testResult, subexpressions.intersect([1], [1]), [1])
	expect(testResult, subexpressions.intersect([1], [2]), [ ])
	expect(testResult, subexpressions.intersect([1], [1, 2]), [1])
	expect(testResult, subexpressions.intersect([1, 2], [1]), [1])
	expect(testResult, subexpressions.intersect([1, 2], [1, 2]), [1, 2])
	expect(testResult, subexpressions.intersect([1, 2, 3], [1]), [1])
	expect(testResult, subexpressions.intersect([1, 2, 3], [2]), [2])
	expect(testResult, subexpressions.intersect([1, 2, 3], [1, 2]), [1, 2])
	expect(testResult, subexpressions.intersect([1, 2, 3], [2, 3]), [2, 3])
	expect(testResult, subexpressions.intersect([1, 2, 3], [1, 3]), [1, 3])
	expect(testResult, subexpressions.intersect([1, 2, 3], [1, 2, 3]), [1, 2, 3])

	return testResult.result_str()

def test_count_subexprs():
	testResult = TestResult(0, 0)

	expect(testResult, subexpressions.count_subexprs({}, [[ ], [ ]], 1), {})
	expect(testResult, subexpressions.count_subexprs({}, [[1], [ ]], 1), {})
	expect(testResult, subexpressions.count_subexprs({}, [[ ], [1]], 1), {})
	expect(testResult, subexpressions.count_subexprs({}, [[1], [1]], 1), {frozenset({1}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1], [1]], 2), {})
	expect(testResult, subexpressions.count_subexprs({}, [[1], [2]], 1), {})
	expect(testResult, subexpressions.count_subexprs({}, [[1], [1, 2]], 1), {frozenset({1}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1], [1, 2]], 2), {})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2], [1]], 1), {frozenset({1}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2], [1]], 2), {})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2], [1, 2]], 1), {frozenset({1, 2}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2], [1, 2]], 2), {frozenset({1, 2}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1]], 1), {frozenset({1}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1]], 2), {})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [2]], 1), {frozenset({2}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [2]], 2), {})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1, 2]], 1), {frozenset({1, 2}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1, 2]], 2), {frozenset({1, 2}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [2, 3]], 1), {frozenset({2, 3}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [2, 3]], 2), {frozenset({2, 3}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1, 3]], 1), {frozenset({1, 3}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1, 3]], 2), {frozenset({1, 3}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1, 2, 3]], 1), {frozenset({1, 2, 3}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1, 2, 3]], 2), {frozenset({1, 2, 3}): 2})
	expect(testResult, subexpressions.count_subexprs({}, [[1, 2, 3], [1, 2, 3]], 3), {frozenset({1, 2, 3}): 2})

	return testResult.result_str()

if __name__ == "__main__":
	print("test_intersect: " + test_intersect())
	print("test_count_subexprs: " + test_count_subexprs())