CADNAPATH=/shared/cadna_c-3.1.11
CPPFLAGS=-g -O2 --std=c++14 -Wall
LDFLAGS=-lm

TARGETS=hh-boost cadna-hh-boost

all: $(TARGETS)

hh-boost: hh-boost.cpp
	$(CXX) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

cadna-hh-boost: hh-boost.cpp
	$(CXX) -DUSE_CADNA $(CPPFLAGS) -I$(CADNAPATH)/include -o $@ $< $(LDFLAGS) -L$(CADNAPATH)/lib -lcadnaC

clean:
	rm -rf hh-boost cadna-hh-boost

