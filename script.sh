#!/bin/bash
# This is a script to measure perfomance
matrices=(50 100 1000);
blockSizes=(4 10 50);
np=(1 2 4);

echo -n "" > times.csv;
echo -n "matrix" >> times.csv;

for n in "${np[@]}"; do
    echo -n ";"$n"" >> times.csv;
done
echo "" >> times.csv;

for m in "${!matrices[@]}"; do
    echo -n "${matrices[m]}" >> times.csv;

    for n in "${np[@]}"; do
        echo -n ";" >> times.csv;
        mpirun -np "$n" ./v0 "${matrices[m]}" "${blockSizes[m]}" >> times.csv;
    done

    echo "" >> times.csv;
done
