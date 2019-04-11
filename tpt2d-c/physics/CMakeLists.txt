set(SOURCES
	integrators.cpp
	morse.cpp
	sampling.cpp
	stickyFn.cpp)

add_library(physics ${SOURCES})
target_link_libraries(physics support nauty)