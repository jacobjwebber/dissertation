set terminal tex
set output "../../../report/diagrams/mem-cap_sphere.tex"
set xlabel "Diameter (Grid Points)"
set autoscale
set ylabel "Percentage (%)"
set style data linespoints

plot "../data/memory-cap_graph" with lines lc rgb "blue" notitle, "../data/memory-cap_graph" lc rgb "blue" notitle
