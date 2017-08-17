set terminal epslatex input color
set output "../../../report/diagrams/mem-cap_sphere.tex"
set xlabel "Diameter (Grid Points)"
set autoscale
set ylabel "$\\frac{S_{ss}}{S_{s}}$"
set style data linespoints

plot "../data/memory-cap_graph" with lines lc rgb "blue" notitle, "../data/memory-cap_graph" lc rgb "blue" notitle, pi/6 lt 2 lc rgb "red" notitle
