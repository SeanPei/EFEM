# fe
project(fe)

# Boost
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS filesystem system)
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
else ()
	message(FATAL_ERROR "Need boost library, but not found!")
endif()

set(CONSTRAINT_HEADERS
	Constraint.h
	ExternalConstraint.h
	ElementConstraint.h
	InertiaConstraint.h
	element/IntegralElement.h
	element/IntegralElement.inl)
set(CONSTRAINT_SOURCES
	Constraint.cpp
	ExternalConstraint.cpp
	ElementConstraint.cpp
	InertiaConstraint.cpp
	element/IntegralElement.cpp
	element/TetrahedraElement.cpp
	element/HexahedraElement.cpp)

convert_filenames_to_full_paths(CONSTRAINT_HEADERS)
convert_filenames_to_full_paths(CONSTRAINT_SOURCES)

add_library(fe ${CONSTRAINT_HEADERS} ${CONSTRAINT_SOURCES})
target_link_libraries(fe core math model io material ${Boost_LIBRARIES} ${EXTERNAL_LIBS})
