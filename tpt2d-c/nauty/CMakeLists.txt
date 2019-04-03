set(SOURCES
	nauty.c
	nautil.c
	naugraph.c
	schreier.c
	naurng.c)

add_library(nauty ${SOURCES})
target_link_libraries(nauty physics support)