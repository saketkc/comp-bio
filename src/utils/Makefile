SOURCES = fasta_reader.cpp timer.cpp
EXECUTABLES = fasta_reader
# Replace suffix using macro
OBJS = $(SOURCES:.cpp=.o)
CXX = g++

# Turn on debug
DEBUG = -g

# Turn on warnings
CXXFLAGS = -Wall $(DEBUG)
LFLAGS = -Wall $(DEBUG)




all: $(EXECUTABLES)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXECUTABLES): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

clean: 
	@-rm -f *.o *~
