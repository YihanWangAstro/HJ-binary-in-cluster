CXX=g++
CXXFLAGS=-Wall -std=c++17 -O3 -pthread

all:flyby flyby-single

flyby: 
	$(CXX) $(CXXFLAGS) -o flyby scattering.cpp

flyby-single: 
	$(CXX) $(CXXFLAGS) -o flyby-single scattering-single.cpp

clean:
	rm -f flyby

