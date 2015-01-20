SOURCES = global_alignment.cpp
##tests/global_alignment_test.cpp
EXECUTABLES = global_alignment global_alignment_test
INCLUDES = ./include
INCLUDEARGS = $(addprefix -I, $(INCLUDES))

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
%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EXECUTABLES): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

test:
	global_alignment_test
clean:
	@-rm -f *.o *~
