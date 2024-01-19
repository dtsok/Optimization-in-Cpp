# Makefile

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall

# Source and executable filenames
SRCS = common/common_file.cpp Nelder-Mead/nelder_mead.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = run

# Rule to build the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule to remove generated files
clean:
	rm -f $(OBJS) $(EXEC)
