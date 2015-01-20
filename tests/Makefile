SOURCES = global_alignment_test.cpp
EXECUTABLES = global_alignment_test
INCLUDES = ../include/cpptest ../utils/ ../src/
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
%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EXECUTABLES): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@


install: $(EXECUTABLES)
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(EXECUTABLES) $(ROOT)/bin
test:
	global_alignment_test
clean:
	@-rm -f *.o *~
