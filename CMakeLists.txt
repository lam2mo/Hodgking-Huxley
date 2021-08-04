CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

PROJECT("Hodgkin-Huxley")

INCLUDE_DIRECTORIES( include /shared/cadna_c-3.1.11/include )

ADD_DEFINITIONS( --std=c++11 )

LINK_DIRECTORIES( /shared/cadna_c-3.1.11/lib )
LINK_LIBRARIES( cadnaC )

ADD_EXECUTABLE(hh-model HH-model src/HH_Model_class)
