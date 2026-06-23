# math
project(math)

set(MATH_HEADERS
	Math.h
	MathType.h)
set(MATH_SOURCES
	Math.cpp
	MathType.cpp)

convert_filenames_to_full_paths(MATH_HEADERS)
convert_filenames_to_full_paths(MATH_SOURCES)

add_library(math ${MATH_HEADERS} ${MATH_SOURCES})
