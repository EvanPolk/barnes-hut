
CXX = g++
SRC = nbody2.cpp
CXXFLAGS = -O3 -march=native -ffast-math

# Output binaries
BIN_SERIAL   = nbody2_serial
BIN_OMP      = nbody2_omp
BIN_OMP_VEC  = nbody2_omp_vec

all: $(BIN_SERIAL) $(BIN_OMP) $(BIN_OMP_VEC)

# 1. Baseline serial
$(BIN_SERIAL): $(SRC)
	$(CXX) -O2 $(SRC) -o $(BIN_SERIAL)

# 2. OpenMP parallel (threads only)
$(BIN_OMP): $(SRC)
	$(CXX) -O2 -fopenmp -DOPENMP $(SRC) -o $(BIN_OMP)

# 3. Combined SIMD + OpenMP
$(BIN_OMP_VEC): $(SRC)
	$(CXX) $(CXXFLAGS) -fopenmp -DSERIAL_VEC -DOPENMP $(SRC) -o $(BIN_OMP_VEC)

clean:
	rm -f $(BIN_SERIAL) $(BIN_OMP) $(BIN_OMP_VEC)
