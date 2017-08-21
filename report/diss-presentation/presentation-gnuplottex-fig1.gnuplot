set terminal epslatex color dashed
set output 'presentation-gnuplottex-fig1.tex'
set xlabel "Diameter (Grid Points)"
set autoscale
set ylabel "$\\frac{S_{ss}}{S_{s}}$"
#set style data linespoints
set size 0.8,0.8

plot "../../semi-structured/results/data/memory-cap_graph" \
with lines lt 1 lc rgb "blue" notitle, \
"../../semi-structured/results/data/memory-cap_graph" \
lc rgb "blue" notitle, \
pi/6 lt 2 lc rgb "red" notitle
