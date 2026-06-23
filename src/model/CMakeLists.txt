# model
project(model)

set(MODEL_HEADERS
	Model.h
	Body.h
	Element.h
	Node.h)
set(MODEL_SOURCES
	Model.cpp
	Body.cpp)

convert_filenames_to_full_paths(MODEL_HEADERS)
convert_filenames_to_full_paths(MODEL_SOURCES)

add_library(model ${MODEL_HEADERS} ${MODEL_SOURCES})
target_link_libraries(model core math io ${VTK_LIBRARIES})
