for i in {87..100}; do
    rm -f proof.txt
    echo "Size $i" >> largeInstances.txt
    echo "Size $i" >> drat-output.txt
    echo "Size $i" 
    xmaplelcm -compute-dip -produce-proof first-grid-$i-$i.cnf >> largeInstances.txt
#    drat-trim first-grid-$i-$i.cnf proof.txt >> drat-output.txt
done
