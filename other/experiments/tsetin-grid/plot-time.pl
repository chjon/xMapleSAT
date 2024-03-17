set terminal push

#set terminal pdfcairo
set terminal pngcairo font "Helvetica,12"  size 500,400
#set output outfile

set grid
#unset key
#set key top left reverse Left
set key top left

set style fill solid
#set boxwidth 0.95
set boxwidth 1.4

#set title "Average Result for 2 Seed"
set title "Tseiting on the grid results"

set xlabel "Problem size"
set ylabel "Running time (in seconds)

#set ytics 5
#set yrange [:20]

#set xtics 20
#set xrange [-100:102]
#set xrange [0:600]

#set xtics 0.5
#set xrange [-0.5:10.5]

#set xtics font ",17"

set output "out-time.png"
plot "data-time.txt" using 1:2 w lines title "Solver Time", 0.5*x**2 w lines title "0.5 n^2", 4*1.105**x w lines title "4*1.105^n"

#Restore the original file or channel name 
set output
#Restore the original terminal setting
set terminal pop


quit


