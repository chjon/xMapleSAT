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
set ylabel "Resolution steps"

#set ytics 5
#set yrange [:20]

#set xtics 20
#set xrange [-100:102]
#set xrange [0:600]

#set xtics 0.5
#set xrange [-0.5:10.5]

#set xtics font ",17"

set output "out-proof.png"
plot "data-proofs-complete.txt" using 1:4 w lines title "Resolution steps", 2000*x**2 title "2000*x^2", 100000*x title "100000*x"
#plot "data-proofs-complete.txt" using 1:3 w lines title "Lemmas in proofs", 230*x**2 title "230*x^2", 15000*x title "15000*x"

#Restore the original file or channel name 
set output
#Restore the original terminal setting
set terminal pop


quit


