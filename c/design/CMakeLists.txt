set(SOURCES
	design.cpp
	)

add_library(tpt ${SOURCES})
target_link_libraries(design tpt nauty physics support)