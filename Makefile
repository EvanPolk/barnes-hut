CXX = g++
CXXFLAGS = -O3 -march=native -ffast-math
SRC = nbody2.cpp

# Output binaries
BIN_SERIAL = nbody2_serial
BIN_VEC    = nbody2_vec
BIN_OMP    = nbody2_omp

all: $(BIN_SERIAL) $(BIN_VEC) $(BIN_OMP)

# 1. Baseline serial
$(BIN_SERIAL): $(SRC)
	$(CXX) -O2 $(SRC) -o $(BIN_SERIAL)

# 2. Serial vectorized
$(BIN_VEC): $(SRC)
	$(CXX) $(CXXFLAGS) -DSERIAL_VEC $(SRC) -o $(BIN_VEC)

# 3. OpenMP parallel
$(BIN_OMP): $(SRC)
	$(CXX) $(CXXFLAGS) -fopenmp -DOPENMP $(SRC) -o $(BIN_OMP)

clean:
	rm -f $(BIN_SERIAL) $(BIN_VEC) $(BIN_OMP)
