n=150
while [ $n -le 500 ]; do
    echo $n
    echo "KISSAT;CRYPTOMINISAT;XMAPLELCM"
    for i in {1..5}; do
	cnfgen randkxor 3 $n $n -T or 2 > s.cnf
	/usr/bin/time -o "time.txt" -f "%e" --quiet ~/recerca/systems/kissat-master-15-sept/build/kissat --time=60 s.cnf > kissat.txt
	kissat=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet cryptominisat5 --maxtime 60 s.cnf > cryptominisat.txt
	crypto=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet ../bin/xmaplelcm -compute-dip -cpu-lim=60 s.cnf > xmaple.txt
	maple=`cat time.txt`
	rm time.txt
	echo -n $kissat
	echo -n ";"
	echo -n $crypto
	echo -n ";"
	echo $maple
    done
    n=$((n+10))
done
