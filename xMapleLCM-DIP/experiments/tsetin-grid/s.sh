for i in {31..100}; do
    cnfgen tseitin first grid $i $i > first-grid-$i-$i.cnf
done
