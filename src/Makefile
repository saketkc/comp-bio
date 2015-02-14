GSOURCES =  global_alignment.cpp  utils/fasta_reader.cpp scoring_matrix.cpp timer.cpp
BSOURCES =  banded_global_alignment.cpp  utils/fasta_reader.cpp scoring_matrix.cpp timer.cpp
DSOURCES =  global_alignment_distancewise.cpp  utils/fasta_reader.cpp scoring_matrix.cpp timer.cpp
MSOURCES =  multiple_alignment.cpp utils/fasta_reader.cpp sequence_profile.cpp sequece_preprocessor.cpp
EXECUTABLES = banded_global_alignment global_alignment global_alignment_distancewise multiple_alignment
INCLUDES = ./utils/
INCLUDEARGS = $(addprefix -I, $(INCLUDES))

# Replace suffix using macro
GOBJS = $(GSOURCES:.cpp=.o)
BOBJS = $(BSOURCES:.cpp=.o)
DOBJS = $(DSOURCES:.cpp=.o)
MOBJS = $(MSOURCES:.cpp=.o)

CXX = g++

# Turn on debug
DEBUG = -g

# Turn on warnings
CXXFLAGS = -Wall -O3 $(DEBUG) $(INCLUDEARGS)


all: banded_global_alignment global_alignment global_alignment_distancewise multiple_alignment

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

banded_global_alignment.o: banded_global_alignment.cpp
	$(CXX) $(CXXFLAGS) -c $<

global_alignment.o: global_alignment.cpp
	$(CXX) $(CXXFLAGS) -c $<

multiple_alignment.o: multiple_alignment.cpp
	$(CXX) $(CXXFLAGS) -c $<

sequence_profile.o: sequence_profile.cpp sequence_profile.hpp
	$(CXX) $(CXXFLAGS) -c $<

sequence_preprocessor.o: sequence_preprocessor.cpp sequence_preprocessor.hpp
	$(CXX) $(CXXFLAGS) -c $<

global_alignment_distancewise.o: global_alignment_distancewise.cpp
	$(CXX) $(CXXFLAGS) -c $<


banded_global_alignment: banded_global_alignment.o utils/fasta_reader.o scoring_matrix.o timer.o
	$(CXX) $(CXXFLAGS) $(BOBJS) -o $@


global_alignment: global_alignment.o utils/fasta_reader.o scoring_matrix.o timer.o
	$(CXX) $(CXXFLAGS) $(GOBJS) -o $@

global_alignment_distancewise: global_alignment_distancewise.o utils/fasta_reader.o scoring_matrix.o timer.o
	$(CXX) $(CXXFLAGS) $(DOBJS) -o $@

multiple_alignment: multiple_alignment.o utils/fasta_reader.o sequence_profile.o sequence_preprocessor.o timer.o
	$(CXX) $(CXXFLAGS) $(GOBJS) -o $@

install: $(EXECUTABLES)
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(EXECUTABLES) $(ROOT)/bin

test:
	global_alignment_test
clean:
	@-rm -f *.o *~
	@-rm -f $(EXECUTABLES)
