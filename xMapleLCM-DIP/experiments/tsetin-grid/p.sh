for i in {41..100}; do
    rm -f proof.txt
    echo "Size $i" >> outNou.txt.txt
    echo "Size $i" outNou.txt.txt
    xmaplelcm -compute-dip -produce-proof first-grid-$i-$i.cnf >> outNou.txt.txt
    ../../drat-trim-master/drat-trim first-grid-$i-$i.cnf proof.txt
done
