set(SOURCES
	hydro.cpp)

add_library(hydro ${SOURCES})
target_link_libraries(hydro support physics nauty)
