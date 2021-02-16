cmake_minimum_required(VERSION 3.16)
project(cupac CUDA CXX)

set(CMAKE_CUDA_STANDARD 14)

#add_executable(cupac main.cu cupac.cu cupac.cuh counter.cu counter.cuh grid.cu grid.cuh define.cuh)

#set_target_properties(
#        cupac
#        PROPERTIES
#        CUDA_SEPARABLE_COMPILATION ON)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_75")

add_library(cupac SHARED src/define.cuh src/grid.cu src/grid.cuh src/counter.cu src/counter.cuh src/cupac.cu src/cupac.cuh)
