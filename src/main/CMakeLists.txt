# main
project(main)

set(MAIN_HEADERS)
set(MAIN_SOURCES
	main.cpp)

convert_filenames_to_full_paths(MAIN_HEADERS)
convert_filenames_to_full_paths(MAIN_SOURCES)

add_executable(main ${MAIN_HEADERS} ${MAIN_SOURCES})
target_link_libraries(main fe core io material math model simulation ${MKL_LIBRARIES})
set_target_properties(main PROPERTIES OUTPUT_NAME "EFEM")
