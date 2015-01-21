SOURCES = global_alignment.cpp utils/fasta_reader.cpp scoring_matrix.cpp
EXECUTABLES = global_alignment
INCLUDES = ../include/ ../utils/
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
%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

global_alignment: fasta_reader.o scoring_matrix.o
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

install: $(EXECUTABLES)
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(EXECUTABLES) $(ROOT)/bin

test:
	global_alignment_test
clean:
	@-rm -f *.o *~
