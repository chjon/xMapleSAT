extract_time () {
    file=$1
    t=`grep "CPU time" $file | awk '{print $4}' #print 4th column`
    echo $t
}

extract_conflicts () {
    file=$1
    t=`grep "conflicts" $file | head -1 | awk '{print $3}' #print 4th column`
    echo $t
}

extract_dip_conflicts () {
    file=$1
    t=`grep "dip conflicts" $file | awk '{print $4}' #print 4th column`
    echo $t
}

extract_decisions () {
    file=$1
    t=`grep "decisions" $file | head -1 | awk '{print $3}' #print 4th column`
    echo $t
}

extract_decisions_on_ext () {
    file=$1
    t=`grep "decisions on ext" $file | awk '{print $6}' #print 4th column`
    echo $t
}

extract_result () {
    file=$1
    t=`egrep "SATISF|INDETER|UNKNOWN" $file | awk '{print $1}' #print 4th column`
    echo $t
}

extract_perc_conflict_with_DIP ( ){
    file=$1
    t=`egrep "DIP found" $file | awk '{print $7}' #print 4th column`
    echo $t | cut -c2-   
}

extract_perc_conflict_with_DIP_learning ( ){
    file=$1
    t=`egrep "DIP-learning" $file | awk '{print $6}' #print 4th column`
    echo $t | cut -c2-      
}

extract_perc_decs_on_extended ( ){
    file=$1
    t=`egrep "on ext vars" $file | awk '{print $7}' #print 4th column`
    echo $t  | cut -c2-     
}

extract_perc_DIP_overhead ( ){
    file=$1
    t=`egrep "DIP computation" $file | awk '{print $7}' #print 4th column`
    echo $t | cut -c2-    
}



for file in $1/*.res.xmaplelcm-common20*; do
    time=$(extract_time $file)
    decs=$(extract_decisions $file)
    decs_on_ext=$(extract_decisions_on_ext $file)
    confs=$(extract_conflicts $file)
    stat=$(extract_result $file)
    perc_decs_on_ext=$(extract_perc_decs_on_extended $file)
    perc_confs_w_dip=$(extract_perc_conflict_with_DIP $file)
    perc_confs_w_dip_learn=$(extract_perc_conflict_with_DIP_learning $file)
    perc_dip_overhead=$(extract_perc_DIP_overhead $file)
    
    name=`basename $file`

    #    echo "$name;$stat;$time;$confs;$dip_confs;$decs"
    echo "$name;$stat;$time;$confs;$decs;$decs_on_ext;$perc_decs_on_ext;$perc_confs_w_dip;$perc_confs_w_dip_learn;$perc_dip_overhead"
done
