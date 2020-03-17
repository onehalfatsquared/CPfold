set(SOURCES
	design.cpp
	designTPT.cpp
	irreversibleDesign.cpp
	)

add_library(design ${SOURCES})
target_link_libraries(design tpt nauty physics support)