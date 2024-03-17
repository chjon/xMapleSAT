n=20
PATH_SBVA="/home/oliveras/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin"
cp albert_run_sbva_cadical $PATH_SBVA
echo "cd $PATH_SBVA && bash albert_run_sbva_cadical 40 20 s.cnf ." > run_sbva.sh
while [ $n -le 200 ]; do
    echo $n
    echo "KISSAT;CRYPTOMINISAT;XMAPLELCM;SBVA_CADICAL"
    for i in {1..5}; do
	cnfgen randkxor 4 $n $n -T or 3 > s.cnf
	/usr/bin/time -o "time.txt" -f "%e" --quiet ~/recerca/systems/kissat-master/build/kissat --time=5 s.cnf > kissat.txt
	kissat=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet cryptominisat5 --maxtime 5 s.cnf > cryptominisat.txt
	crypto=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet ../../bin/xmaplelcm -cpu-lim=5 s.cnf > xmaple.txt
	maple=`cat time.txt`	
	rm time.txt
	cp s.cnf ~/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin/.
	/usr/bin/time -o "time.txt" -f "%e" --quiet timeout 60 bash run_sbva.sh > sbva.txt
	sbva=`cat time.txt`	
	rm time.txt	
	echo -n $kissat
	echo -n ";"
	echo -n $crypto
	echo -n ";"
	echo -n $maple
	echo -n ";"
	echo $sbva
	echo
    done
    n=$((n+5))
done
