set(SOURCES
	genetics.cpp
	genetics_sampling.cpp)

add_library(genetic ${SOURCES})
target_link_libraries(genetic design support)
