all: dep su2 su2hmc

CXX=g++ 
INCLUDE=
CXXFLAGS=-Wall -Wpedantic --std=c++11 -O3

MODULES = gaugeconfig su2 gauge_energy random_gauge_trafo get_staples sweep wilsonloop expsu2

-include $(addsuffix .d,$(MODULES) main hmc)

dep: $(addsuffix .d,$(MODULES) main hmc)
	@ echo "...dependency files built"

$(addsuffix .d, $(MODULES) main hmc): %.d: %.cc Makefile
	@ $(CXX) -MM $(CXXFLAGS) ${INCLUDE} $< > $@

$(addsuffix .o, $(MODULES) main hmc): %.o: %.cc %.d Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

su2: $(addsuffix .o, $(MODULES) main) Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(addsuffix .o, $(MODULES)) main.o -o su2

su2hmc: $(addsuffix .o, $(MODULES) hmc) Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(addsuffix .o, $(MODULES)) hmc.o -o su2hmc

clean:
	rm -f su2 su2hmc *.o *~ core
