# simulation
project(simulation)

set(SIMULATION_HEADERS
	Integrator.h
	Simulator.h
	nlsolver/LBFGS.h
	nlsolver/LBFGSMat.h
	nlsolver/LineSearchFEA.h
	nlsolver/Param.h)
set(SIMULATION_SOURCES
	Integrator.cpp
    Simulator.cpp)

convert_filenames_to_full_paths(SIMULATION_HEADERS)
convert_filenames_to_full_paths(SIMULATION_SOURCES)

add_library(simulation ${SIMULATION_HEADERS} ${SIMULATION_SOURCES})
target_link_libraries(simulation core math model io fe)