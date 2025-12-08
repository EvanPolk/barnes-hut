#!/bin/bash

echo "Running tests for Barnes Hut"
echo ""

echo "single thread"
for i in 128000 256000 512000 1024000; do
	./nbody2_serial -n $i | grep Average
done
echo ""

echo "OMP 2,4,8,16,32"
for i in 2 4 8 16 32; do
	export OMP_NUM_THREADS=$i
	for j in 128000 256000 512000 1024000; do
		./nbody2_omp_vec -n $j | grep Average
	done
done
