for i in {7..40}; do
    for file in tseitin-d6-$i.cnf; do
	echo $file
	# ~/recerca/systems/kissat --time=600 $file > $file".kissat"
	#../../bin/xmaplelcm -cpu-lim=5000 $file > $file".DIP"
	timeout 600 glucoser $file > $file".real-GLUCOSER" 2>&1
	# cryptominisat5 --maxtime 600 $file > $file".crypto"
	# /usr/bin/time -o "time.txt" -f "%e" --quiet timeout 600 bash run_svba.sh  $file . > $file".sbva"
	# sbva=`cat time.txt`
	# echo "$sbva" > $file".sbva"
	# rm time.txt
    done
done
