SOURCES = global_alignment.cpp
EXECUTABLES = global_alignment
INCLUDES = ../include
INCLUDEARGS = $(addprefix -I,$(INCLUDES))

# Replace suffix using macro
OBJS = $(SOURCES:.cpp=.o)
CXX = g++

# Turn on debug
DEBUG = -g

# Turn on warnings
CXXFLAGS = -Wall -O3 $(DEBUG)


all: $(EXECUTABLES)

## $<: Name of prerequisite (*.cpp)
## $@: Name of target (*.op)
%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@  $(INCLUDEARGS)

$(EXECUTABLES): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

clean:
	@-rm -f *.o *~
