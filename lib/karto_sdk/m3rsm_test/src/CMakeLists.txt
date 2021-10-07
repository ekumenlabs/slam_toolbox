# Include paths.
include_directories(
	../include
)

# Test sources.
set (GTEST_SOURCES
	M3RScanMatcher_TEST.cpp
)

m3rsm_tests_build(${GTEST_SOURCES})
