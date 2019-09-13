set(SOURCES
	tpt.cpp
	)

add_library(tpt ${SOURCES})
target_link_libraries(tpt nauty physics support)