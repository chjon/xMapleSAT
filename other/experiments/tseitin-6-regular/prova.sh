n=40
while [ $n -le 5000 ]; do
    echo "Size $n"
    echo "Size $n" >> outProva.txt
    cnfgen tseitin $n 6 > tseitin-prova-d6-$n.cnf
    xmaplelcm -dip-pair-min=20 -cpu-lim=2000 tseitin-prova-d6-$n.cnf | grep "CPU"
    n=$((n+5))
done
