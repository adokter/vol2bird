cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(mistnet)

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

add_library(mistnet SHARED libmistnet.cpp)
target_link_libraries(mistnet "${TORCH_LIBRARIES}")
set_property(TARGET mistnet PROPERTY CXX_STANDARD 14)
set_property(TARGET mistnet PROPERTY SUFFIX ".so")
install(TARGETS mistnet DESTINATION lib)
