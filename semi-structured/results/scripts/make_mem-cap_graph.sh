#/bin/bash

cd ../../memory-capacity/
make
cd -
../../memory-capacity/memory-capacity > ../data/memory-cap_raw
sed -n 5,24p ../data/memory-cap_raw | cut -f1,5 > ../data/memory-cap_graph
gnuplot mem-cap_graph.gp
sed -n 4,24p ../data/memory-cap_raw | sed 's/\t/,/g'  > ../data/memory-cap_sphere.csv
sed -n 27,60p ../data/memory-cap_raw | sed 's/\t/,/g' > ../data/memory-cap_cube.csv
