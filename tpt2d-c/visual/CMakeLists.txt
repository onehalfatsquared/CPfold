set(SOURCES
	graph.cpp)

add_library(visual ${SOURCES})
target_link_libraries(visual support)