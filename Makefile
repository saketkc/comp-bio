SOURCES = global_alignment.cpp
EXECUTABLES = global_alignment
# Replace suffix using macro
OBJS = $(SOURCES:.cpp=.o)
CXX = g++

# Turn on debug
DEBUG = -g

# Turn on warnings
CXXFLAGS = -Wall $(DEBUG)


all: $(EXECUTABLES)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXECUTABLES): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

clean:
	@-rm -f *.o *~
