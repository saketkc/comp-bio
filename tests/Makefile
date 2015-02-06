SOURCES = global_alignment_test.cpp test_fasta_reader.cpp
EXECUTABLES = global_alignment_test test_fasta_reader
INCLUDES =   ../src/ ../src/utils/
INCLUDEARGS = $(addprefix -I, $(INCLUDES))


# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall
# Replace suffix using macro
OBJS = $(SOURCES:.cpp=.o)
CXX = g++

# Turn on debug
DEBUG = -g

# Turn on warnings
CXXFLAGS = -Wall -O3 $(DEBUG) $(INCLUDEARGS)


all: $(EXECUTABLES)

## $<: Name of prerequisite (*.cpp)
## $@: Name of target (*.op)
#%.o: %.cpp %.hpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@
#%: %.cpp
#	$(CXX) $(CXXFLAGS) -o $@ $^

#$(EXECUTABLES): $(OBJS)
#	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

global_alignment_test.o: global_alignment_test.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

global_alignment_test: global_alignment_test.o
	$(CXX) $(CXXFLAGS) global_alignment_test.o -o $@ 

test_fasta_reader.o: test_fasta_reader.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test_fasta_reader: test_fasta_reader.o
	$(CXX) $(CXXFLAGS) test_fasta_reader.o -o $@ 


install: $(EXECUTABLES)
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(EXECUTABLES) $(ROOT)/bin
test:
	./global_alignment_test 
	./test_fasta_reader
clean:
	@-rm -f *.o *~
