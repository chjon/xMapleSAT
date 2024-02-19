set terminal push

#set terminal pdfcairo
set terminal pngcairo font "Helvetica,10" size 1100,900  # size: width hight
#set output outfile

set grid
#unset key
set key top left reverse Left


set xlabel "Original XMaple-LCM"
set ylabel "DIP disabling XMapleLCM"

#set ytics 0.1
#set logscale xy 1.005
set xrange [10:1810]
set yrange [10:1810]
#set trange [0:1800]

#set logscale xy 2
#set tics scale 0  # same as above?
#set xtics 5

#set xtics auto
#set ytics auto

set xtics border out scale 1,0.5 nomirror
set log x

set ytics border out scale 1,0.5 nomirror
set log y


set title "MapleLCM vs DIP-MapleLCM (SAT 2023 benchs)"
set output "plot_SAT_2034.png"

plot "times.txt" u 1:2 w point pt 7 ps 0.7 lc rgb "red" notitle, x dashtype 2 lc rgb "black", 2*x dashtype 2 lc rgb "blue", 0.5*x dashtype 2 lc rgb "blue"



#Restore the original file or channel name 
set output
#Restore the original terminal setting
set terminal pop


quit

