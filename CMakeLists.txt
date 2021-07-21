CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

PROJECT("Hodgkin-Huxley")

INCLUDE_DIRECTORIES( include /usr/local/Cellar/boost/1.76.0/include )

ADD_DEFINITIONS( --std=c++11 )

add_library(hhmod SHARED 
	include/HH_Model_class.hpp
	src/HH_Model_class.cpp
)

ADD_EXECUTABLE(hh-model HH-model )

target_link_libraries(hh-model hhmod)
