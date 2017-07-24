#/bin/bash

cd ../../memory-capacity/
make
cd -
../../memory-capacity/memory-capacity > ../data/memory-cap_raw
sed -n 5,24p ../data/memory-cap_raw | cut -f1,5 > ../data/memory-cap_graph
gnuplot mem-cap_graph.gp
