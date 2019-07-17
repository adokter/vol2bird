cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(mistnet)

find_package(Torch REQUIRED)

add_library(mistnet STATIC libmistnet.cpp)
target_link_libraries(mistnet "${TORCH_LIBRARIES}" z)
set_property(TARGET mistnet PROPERTY CXX_STANDARD 11)
install(TARGETS mistnet DESTINATION ../lib)
