CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native -DNDEBUG
TARGET = linear_algebra
SRCS = main.cpp Matrix.cpp Analysis.cpp Decomposer.cpp LinearSolver.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
