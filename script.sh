#!/bin/bash
# This is a script to measure perfomance
# Define matrices to be tested
mtx=(50 100 1000);
# Define block size for each matrix
blockSizes=(4 10 50);
# Define number of proccesses to be tested in v3 (MPI)
np=(1 2 4);

# v0
echo "matrix;1" > times.csv;
for m in "${mtx[@]}"; do
    echo -n "$m;" >> times.csv;
    ./v0 $m >> times.csv;
done

echo "" >> times.csv

# v1
echo "matrix;1" >> times.csv;
for m in "${mtx[@]}"; do
    echo -n "$m;" >> times.csv;
    ./v1 $m >> times.csv;
done

echo "" >> times.csv

# v2
echo "matrix;1" >> times.csv;
for m in "${!mtx[@]}"; do
    echo -n "$m;" >> times.csv;
    ./v2 ${mtx[m]} ${blockSizes[m]} >> times.csv;
done

echo "" >> times.csv;

# v3
echo -n "matrix" >> times.csv;
for n in "${np[@]}"; do
    echo -n ";"$n"" >> times.csv;
done
echo "" >> times.csv;

for m in "${!mtx[@]}"; do
    echo -n "${mtx[m]}" >> times.csv
    for n in "${np[@]}"; do
        echo -n ";" >> times.csv;
        mpirun -np "$n" ./v3 "${mtx[m]}" "${blockSizes[m]}" >> times.csv;
    done
    echo "" >> times.csv;
done
