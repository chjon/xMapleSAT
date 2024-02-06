for i in {7..100}; do
    cnfgen tseitin $i 6 > tseitin-d6-$i.cnf    

    echo "Size $i"
    
    echo "Size $i" >> kissat.txt
    ./kissat --time=300 tseitin-d6-$i.cnf >> kissat.txt

        echo "Size $i" >> maple.txt
    ../bin/xmaplelcm -cpu-lim=300 tseitin-d6-$i.cnf >> maple.txt

done
