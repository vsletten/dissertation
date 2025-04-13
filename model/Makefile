CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -ffast-math
#CXXFLAGS = -std=c++11 -Wall -Wextra -ggdb
LDFLAGS = -lm
SOURCES = mckaol.cpp actions.cpp bfsearch.cpp envrn.cpp evtlist.cpp futil.cpp lattice.cpp myerr.cpp output.cpp ran2.cpp reactions.cpp rxnlist.cpp sim.cpp ucell.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = mckaol
HEADERS = $(wildcard *.hpp)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

.PHONY: all clean
