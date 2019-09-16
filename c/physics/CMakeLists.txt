set(SOURCES
	stickyFn.cpp
	integrators.cpp
	morse.cpp
	lennardJones.cpp
	sampling.cpp
	mcm.cpp)

add_library(physics ${SOURCES})
target_link_libraries(physics support nauty)
