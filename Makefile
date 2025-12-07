# Compiler
CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -march=native

# Target executable
TARGET = nbody2

# Source files
SRC = nbody2.cpp

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean build files
clean:
	rm -f $(TARGET) *.o

# Run the program
run: $(TARGET)
	./$(TARGET)
