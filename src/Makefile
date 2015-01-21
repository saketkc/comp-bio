SOURCES = global_alignment.cpp utils/fasta_reader.cpp scoring_matrix.cpp timer.cpp
EXECUTABLES = global_alignment
INCLUDES = ./utils/
INCLUDEARGS = $(addprefix -I, $(INCLUDES))

# Replace suffix using macro
OBJS = $(SOURCES:.cpp=.o)
CXX = g++

# Turn on debug
DEBUG = -g

# Turn on warnings
CXXFLAGS = -Wall -O3 $(DEBUG) $(INCLUDEARGS)


all: global_alignment

## $<: Name of prerequisite (*.cpp)
## $@: Name of target (*.op)
#%.o: %.cpp %.hpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@

fasta_reader.o: utils/fasta_reader.cpp utils/fasta_reader.hpp
	$(CXX) $(CXXFLAGS) -c $<

scoring_matrix.o: scoring_matrix.cpp scoring_matrix.hpp
	$(CXX) $(CXXFLAGS) -c $<

timer.o: utils/timer.cpp utils/timer.hpp
	$(CXX) $(CXXFLAGS) -c $<

global_alignment.o: global_alignment.cpp
	$(CXX) $(CXXFLAGS) -c $<


global_alignment: global_alignment.o utils/fasta_reader.o scoring_matrix.o timer.o
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

install: $(EXECUTABLES)
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(EXECUTABLES) $(ROOT)/bin

test:
	global_alignment_test
clean:
	@-rm -f *.o *~
