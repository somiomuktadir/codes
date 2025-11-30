CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native -DNDEBUG -Iinclude
LDFLAGS = -lm

# Directories
SRC_DIR = src
INCLUDE_DIR = include
OBJ_DIR = obj
BIN_DIR = bin

# Target executable
TARGET = $(BIN_DIR)/linear_algebra

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

# Library detection
LIBS =
DEFINES =

# Check for BLAS/LAPACK
ifneq ($(shell ldconfig -p | grep -q libblas && echo $$?),0)
    HAS_BLAS := $(shell echo "int main(){return 0;}" | $(CXX) -x c++ - -lblas -o /dev/null 2>/dev/null && echo yes || echo no)
    ifeq ($(HAS_BLAS),yes)
        LIBS += -lblas -llapack
        DEFINES += -DUSE_BLAS -DUSE_LAPACK
        $(info [INFO] BLAS/LAPACK libraries detected and will be used)
    else
        $(info [INFO] BLAS/LAPACK libraries not found, using native implementation)
    endif
endif

# Check for OpenMP
HAS_OPENMP := $(shell echo "int main(){return 0;}" | $(CXX) -x c++ -fopenmp - -o /dev/null 2>/dev/null && echo yes || echo no)
ifeq ($(HAS_OPENMP),yes)
    CXXFLAGS += -fopenmp
    $(info [INFO] OpenMP support enabled for parallel matrix multiplication)
else
    $(info [INFO] OpenMP not available)
endif

# Default target
all: $(TARGET)
	@echo "Build complete. Executable is at $(TARGET)"

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
	@echo "Cleaned build artifacts"

# Build without optional libraries
basic:
	@echo "Building basic version without BLAS/LAPACK/OpenMP..."
	$(MAKE) DEFINES="" LIBS="" CXXFLAGS="-std=c++17 -Wall -Wextra -O3 -march=native -DNDEBUG -Iinclude"

# Build with all optimizations (force)
full:
	@echo "Building full version with all optimizations..."
	$(MAKE) DEFINES="-DUSE_BLAS -DUSE_LAPACK" LIBS="-lblas -llapack" CXXFLAGS="-std=c++17 -Wall -Wextra -O3 -march=native -DNDEBUG -Iinclude -fopenmp"

.PHONY: all basic full clean
