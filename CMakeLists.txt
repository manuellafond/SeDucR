cmake_minimum_required(VERSION 3.10)
project(segrec)



set(CMAKE_CXX_STANDARD 20)



SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -O3 -lz -lm -Wall -march=native")

set(CMAKE_BUILD_TYPE Release)



add_executable(ybrec main.cpp genespeciestreeutil.h genespeciestreeutil.cpp ybrec.h treeutils/define.h treeutils/newicklex.h treeutils/newicklex.cpp treeutils/node.h treeutils/treeutil.h treeutils/util.h treeutils/ewah/ewah.h)

