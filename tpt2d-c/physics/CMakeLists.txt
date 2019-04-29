set(SOURCES
	integrators.cpp
	morse.cpp
	lennardJones.cpp
	sampling.cpp
	stickyFn.cpp)

add_library(physics ${SOURCES})
target_link_libraries(physics support nauty)