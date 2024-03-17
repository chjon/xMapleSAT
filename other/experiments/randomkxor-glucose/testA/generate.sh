for i in {1..100}; do
    echo "SEED $i"
    bash gen_and_test.sh 280 200 30 $i
done
   
