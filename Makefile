all: dep su2

CXX=g++ 
INCLUDE=
CXXFLAGS=-Wall -Wpedantic -std=c++11 -g -O3

MODULES = gaugeconfig su2 gauge_energy random_gauge_trafo get_staples sweep main

-include $(addsuffix .d,$(MODULES))

dep: $(addsuffix .d,$(MODULES))
	@ echo "...dependency files built"

$(addsuffix .d, $(MODULES)): %.d: %.cc Makefile
	@ $(CXX) -MM $(CXXFLAGS) ${INCLUDE} $< > $@

$(addsuffix .o, $(MODULES)): %.o: %.cc %.d Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

su2: $(addsuffix .o, $(MODULES)) Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(addsuffix .o, $(MODULES)) -o su2

clean:
	rm -f su2 *.o *~ core
