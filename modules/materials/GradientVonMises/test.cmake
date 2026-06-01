# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for GradientVonMises material
add_marmot_test("TestGradientVonMises" "${CURR_TEST_SOURCE_DIR}/test.cpp")
