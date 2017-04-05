#set xrange [-1:1]
#set yrange [-1:1]
set ticslevel 0
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid ztics lt 0 lw 1 lc rgb "#bbbbbb"
#set grid
set border 4095

set hidden3d
set terminal pdf
set output 'locations.pdf'
splot "actual_locations.dat" using 1:2:3 t "Sensors", \
    "anchor_locations.dat" using 1:2:3 t "Anchors"
