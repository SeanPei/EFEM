# io
project(io)

set(IO_HEADERS
	IO.h
	Config.h
	MeshParser.h
	Logger.h)
set(IO_SOURCES
	IO.cpp
	Config.cpp
	Logger.cpp)

convert_filenames_to_full_paths(IO_HEADERS)
convert_filenames_to_full_paths(IO_SOURCES)

add_library(io ${IO_HEADERS} ${IO_SOURCES})
target_link_libraries(io ${Boost_LIBRARIES})
