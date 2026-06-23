# core
project(core)

set(CORE_HEADERS
	IdGenerator.h
	Object.h)
set(CORE_SOURCES
	Object.cpp)

convert_filenames_to_full_paths(CORE_HEADERS)
convert_filenames_to_full_paths(CORE_SOURCES)

add_library(core ${CORE_HEADERS} ${CORE_SOURCES})
