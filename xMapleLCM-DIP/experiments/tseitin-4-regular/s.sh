for i in {10..100}; do
    cnfgen tseitin $i > tseitin-$i.cnf    

    echo "Size $i"
    
    echo "Size $i" >> kissat.txt
    ./kissat --time=300 tseitin-$i.cnf >> kissat.txt

        echo "Size $i" >> maple.txt
    ../bin/xmaplelcm -cpu-lim=300 tseitin-$i.cnf >> maple.txt

done
