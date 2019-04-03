add_executable(runBD runBD.cpp)
add_executable(runEstimator runEstimator.cpp)
add_executable(runTPT runTPT.cpp)

target_link_libraries(runBD physics support)
target_link_libraries(runEstimator physics support)
target_link_libraries(runTPT tpt nauty physics support)