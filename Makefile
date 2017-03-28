all: su2

CXX=g++ -std=c++11 -g
INCLUDE=-Iinclude/
CXXFLAGS=-Wpedantic

gaugeconfig.o: gaugeconfig.cc gaugeconfig.hh su2.hh random_su2.hh Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

su2.o: su2.cc su2.hh Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

gauge_energy.o: gauge_energy.cc gauge_energy.hh gaugeconfig.hh Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

random_gauge_trafo.o: random_gauge_trafo.cc random_gauge_trafo.hh su2.hh random_su2.hh Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

su2: main.cc gaugeconfig.o su2.o gauge_energy.o random_gauge_trafo.o include/geometry.hh Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) main.cc gaugeconfig.o gauge_energy.o random_gauge_trafo.o su2.o -o su2

clean:
	rm -f su2 *.o *~ core
