cmake_minimum_required(VERSION 2.8)
project(comet-test)

set(MAKEMS2_DIR "../comet-ms-code/comet-ms/CometSearch/PeptideFragmentation/MakeMS2/")

INCLUDE_DIRECTORIES(${MAKEMS2_DIR})
INCLUDE_DIRECTORIES("/usr/include/eigen3")

aux_source_directory(. DIR_LIST)
aux_source_directory(${MAKEMS2_DIR} DIR_LIST)
aux_source_directory(${MAKEMS2_DIR}FragmentModels/ DIR_LIST)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include/eigen3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include -I/usr/local/include/eigen3 -I/usr/local/include/boost")

add_executable(${PROJECT_NAME} ${DIR_LIST})
