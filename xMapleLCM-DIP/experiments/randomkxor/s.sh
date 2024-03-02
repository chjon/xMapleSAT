echo "cnfgen randkxor 3 n n -T or 3"
n=150
#PATH_SBVA="/home/oliveras/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin"
PATH_SBVA="/home/oliveras/recerca/systems/SBVA/bin"
cp albert_run_sbva_cadical $PATH_SBVA
echo "cd $PATH_SBVA && bash albert_run_sbva_cadical 80 40 s.cnf ." > run_sbva.sh
while [ $n -le 500 ]; do
    echo $n
    echo "KISSAT;CRYPTOMINISAT;XMAPLELCM;SBVA_CADICAL"
    for i in {1..5}; do
	cnfgen randkxor 3 $n $n -T or 3 > s.cnf	
	/usr/bin/time -o "time.txt" -f "%e" --quiet ~/recerca/systems/kissat-master/build/kissat --time=300 s.cnf > kissat.txt
	kissat=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet cryptominisat5 --maxtime 300 s.cnf > cryptominisat.txt
	crypto=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet ../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=300 s.cnf > xmaple.txt
	maple=`cat time.txt`	
	rm time.txt
	cp s.cnf $PATH_SBVA/.
	/usr/bin/time -o "time.txt" -f "%e" --quiet timeout 300 bash run_sbva.sh > sbva.txt
	sbva=`cat time.txt`	
	rm time.txt	
	echo -n $kissat
	echo -n ";"
	echo -n $crypto
	echo -n ";"
	echo -n $maple
	echo -n ";"
	echo $sbva
    done
    n=$((n+10))
done

#######################################################

echo "cnfgen randkxor 3 n n -T xor 2"
n=150
#PATH_SBVA="/home/oliveras/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin"
PATH_SBVA="/home/oliveras/recerca/systems/SBVA/bin"
cp albert_run_sbva_cadical $PATH_SBVA
echo "cd $PATH_SBVA && bash albert_run_sbva_cadical 80 40 s.cnf ." > run_sbva.sh
while [ $n -le 500 ]; do
    echo $n
    echo "KISSAT;CRYPTOMINISAT;XMAPLELCM;SBVA_CADICAL"
    for i in {1..5}; do
	cnfgen randkxor 3 $n $n -T xor 2 > s.cnf	
	/usr/bin/time -o "time.txt" -f "%e" --quiet ~/recerca/systems/kissat-master/build/kissat --time=300 s.cnf > kissat.txt
	kissat=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet cryptominisat5 --maxtime 300 s.cnf > cryptominisat.txt
	crypto=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet ../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=300 s.cnf > xmaple.txt
	maple=`cat time.txt`	
	rm time.txt
	cp s.cnf $PATH_SBVA/.
	/usr/bin/time -o "time.txt" -f "%e" --quiet timeout 300 bash run_sbva.sh > sbva.txt
	sbva=`cat time.txt`	
	rm time.txt	
	echo -n $kissat
	echo -n ";"
	echo -n $crypto
	echo -n ";"
	echo -n $maple
	echo -n ";"
	echo $sbva
    done
    n=$((n+10))
done

#######################################################

echo "cnfgen randkxor 3 n n -T or 2"
n=150
#PATH_SBVA="/home/oliveras/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin"
PATH_SBVA="/home/oliveras/recerca/systems/SBVA/bin"
cp albert_run_sbva_cadical $PATH_SBVA
echo "cd $PATH_SBVA && bash albert_run_sbva_cadical 80 40 s.cnf ." > run_sbva.sh
while [ $n -le 800 ]; do
    echo $n
    echo "KISSAT;CRYPTOMINISAT;XMAPLELCM;SBVA_CADICAL"
    for i in {1..5}; do
	cnfgen randkxor 3 $n $n -T or 2 > s.cnf	
	/usr/bin/time -o "time.txt" -f "%e" --quiet ~/recerca/systems/kissat-master/build/kissat --time=300 s.cnf > kissat.txt
	kissat=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet cryptominisat5 --maxtime 300 s.cnf > cryptominisat.txt
	crypto=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet ../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=300 s.cnf > xmaple.txt
	maple=`cat time.txt`	
	rm time.txt
	cp s.cnf $PATH_SBVA/.
	/usr/bin/time -o "time.txt" -f "%e" --quiet timeout 300 bash run_sbva.sh > sbva.txt
	sbva=`cat time.txt`	
	rm time.txt	
	echo -n $kissat
	echo -n ";"
	echo -n $crypto
	echo -n ";"
	echo -n $maple
	echo -n ";"
	echo $sbva
    done
    n=$((n+10))
done


#######################################################

echo "cnfgen randkxor 3 n n -T xor 3"
n=20
#PATH_SBVA="/home/oliveras/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin"
PATH_SBVA="/home/oliveras/recerca/systems/SBVA/bin"
cp albert_run_sbva_cadical $PATH_SBVA
echo "cd $PATH_SBVA && bash albert_run_sbva_cadical 80 40 s.cnf ." > run_sbva.sh
while [ $n -le 300 ]; do
    echo $n
    echo "KISSAT;CRYPTOMINISAT;XMAPLELCM;SBVA_CADICAL"
    for i in {1..5}; do
	cnfgen randkxor 3 $n $n -T xor 3 > s.cnf	
	/usr/bin/time -o "time.txt" -f "%e" --quiet ~/recerca/systems/kissat-master/build/kissat --time=300 s.cnf > kissat.txt
	kissat=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet cryptominisat5 --maxtime 300 s.cnf > cryptominisat.txt
	crypto=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet ../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=300 s.cnf > xmaple.txt
	maple=`cat time.txt`	
	rm time.txt
	cp s.cnf $PATH_SBVA/.
	/usr/bin/time -o "time.txt" -f "%e" --quiet timeout 300 bash run_sbva.sh > sbva.txt
	sbva=`cat time.txt`	
	rm time.txt	
	echo -n $kissat
	echo -n ";"
	echo -n $crypto
	echo -n ";"
	echo -n $maple
	echo -n ";"
	echo $sbva
    done
    n=$((n+10))
done


#######################################################

echo "cnfgen randkxor 4 n n -T or 3"
n=20
#PATH_SBVA="/home/oliveras/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin"
PATH_SBVA="/home/oliveras/recerca/systems/SBVA/bin"
cp albert_run_sbva_cadical $PATH_SBVA
echo "cd $PATH_SBVA && bash albert_run_sbva_cadical 80 40 s.cnf ." > run_sbva.sh
while [ $n -le 300 ]; do
    echo $n
    echo "KISSAT;CRYPTOMINISAT;XMAPLELCM;SBVA_CADICAL"
    for i in {1..5}; do
	cnfgen randkxor 4 $n $n -T or 3 > s.cnf	
	/usr/bin/time -o "time.txt" -f "%e" --quiet ~/recerca/systems/kissat-master/build/kissat --time=300 s.cnf > kissat.txt
	kissat=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet cryptominisat5 --maxtime 300 s.cnf > cryptominisat.txt
	crypto=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet ../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=300 s.cnf > xmaple.txt
	maple=`cat time.txt`	
	rm time.txt
	cp s.cnf $PATH_SBVA/.
	/usr/bin/time -o "time.txt" -f "%e" --quiet timeout 300 bash run_sbva.sh > sbva.txt
	sbva=`cat time.txt`	
	rm time.txt	
	echo -n $kissat
	echo -n ";"
	echo -n $crypto
	echo -n ";"
	echo -n $maple
	echo -n ";"
	echo $sbva
    done
    n=$((n+10))
done

#######################################################

echo "cnfgen randkxor 4 n n -T xor 3"
n=20
#PATH_SBVA="/home/oliveras/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin"
PATH_SBVA="/home/oliveras/recerca/systems/SBVA/bin"
cp albert_run_sbva_cadical $PATH_SBVA
echo "cd $PATH_SBVA && bash albert_run_sbva_cadical 80 40 s.cnf ." > run_sbva.sh
while [ $n -le 300 ]; do
    echo $n
    echo "KISSAT;CRYPTOMINISAT;XMAPLELCM;SBVA_CADICAL"
    for i in {1..5}; do
	cnfgen randkxor 4 $n $n -T xor 3 > s.cnf	
	/usr/bin/time -o "time.txt" -f "%e" --quiet ~/recerca/systems/kissat-master/build/kissat --time=300 s.cnf > kissat.txt
	kissat=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet cryptominisat5 --maxtime 300 s.cnf > cryptominisat.txt
	crypto=`cat time.txt`
	rm time.txt
	/usr/bin/time -o "time.txt" -f "%e" --quiet ../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=300 s.cnf > xmaple.txt
	maple=`cat time.txt`	
	rm time.txt
	cp s.cnf $PATH_SBVA/.
	/usr/bin/time -o "time.txt" -f "%e" --quiet timeout 300 bash run_sbva.sh > sbva.txt
	sbva=`cat time.txt`	
	rm time.txt	
	echo -n $kissat
	echo -n ";"
	echo -n $crypto
	echo -n ";"
	echo -n $maple
	echo -n ";"
	echo $sbva
    done
    n=$((n+10))
done
