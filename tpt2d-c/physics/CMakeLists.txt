set(SOURCES
	adjMat.cpp
	euDist.cpp
	integrators.cpp
	morse.cpp
	sampling.cpp
	stickyFn.cpp)

add_library(physics ${SOURCES})
target_link_libraries(physics support)