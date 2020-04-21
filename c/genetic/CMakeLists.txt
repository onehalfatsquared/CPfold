set(SOURCES
	genetics.cpp)

add_library(genetic ${SOURCES})
target_link_libraries(genetic design support)
