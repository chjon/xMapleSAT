extract_time () {
    file=$1
    t=`grep "CPU time" $file | awk '{print $5}' #print 4th column`
    echo $t
}

extract_conflicts () {
    file=$1
    t=`grep "conflicts" $file | head -1 | awk '{print $4}' #print 4th column`
    echo $t
}

extract_decisions () {
    file=$1
    t=`grep "decisions" $file | head -1 | awk '{print $4}' #print 4th column`
    echo $t
}

extract_result () {
    file=$1
    t=`egrep "SATISF|UNKNO" $file | awk '{print $2}' #print 4th column`
    echo $t
}


for file in $1/*.res.maplelcm-2; do
    time=$(extract_time $file)
    decs=$(extract_decisions $file)
    confs=$(extract_conflicts $file)
    stat=$(extract_result $file)
    name=`basename $file`
    echo "$name;$stat;$time;$confs;$decs"
done
