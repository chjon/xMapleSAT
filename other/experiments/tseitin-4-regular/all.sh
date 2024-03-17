for file in tseitin*.cnf; do
    echo $file
    ~/recerca/systems/kissat --time=600 $file > $file".kissat"
    ../../bin/xmaplelcm -cpu-lim=600 $file > $file".DIP"
    cryptominisat5 --maxtime 600 $file > $file".crypto"
    /usr/bin/time -o "time.txt" -f "%e" --quiet timeout 600 bash run_svba.sh  $file . > $file".sbva"
    sbva=`cat time.txt`
    echo "$sbva" > $file".sbva"
    rm time.txt
done
