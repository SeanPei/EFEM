# material
project(material)

set(MATERIAL_HEADERS
	Material.h
)
set(MATERIAL_SOURCES
	NeoHookean.cpp
	MooneyRivlin.cpp
	Yeoh.cpp)

convert_filenames_to_full_paths(MATERIAL_HEADERS)
convert_filenames_to_full_paths(MATERIAL_SOURCES)

add_library(material ${MATERIAL_HEADERS} ${MATERIAL_SOURCES})
target_link_libraries(material math)
