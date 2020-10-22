set(SOURCES
	protocol.cpp
	protocol_vec.cpp
	)

add_library(non_eq_protocol ${SOURCES})
target_link_libraries(non_eq_protocol dlib::dlib design tpt nauty physics support)
