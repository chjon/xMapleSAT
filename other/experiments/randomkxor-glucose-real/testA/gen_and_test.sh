n=$1
seed=$4

echo "kxor-3-xor-3 with size $1"
cnfgen randkxor 3 $1 $1 -T xor 3 > kxor-3-$1-$1-xor-3-seed-$seed.cnf
../../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=1800 kxor-3-$1-$1-xor-3-seed-$seed.cnf


# echo "kxor-4-or-3 with size $2"
# cnfgen randkxor 4 $2 $2 -T or 3 > kxor-4-$2-$2-or-3-seed-$seed.cnf
#  ../../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=1800 kxor-4-$2-$2-or-3-seed-$seed.cnf

# echo "kxor-4-xor-3 with size $3"
# cnfgen randkxor 4 $3 $3 -T xor 3 > kxor-4-$3-$3-xor-3-seed-$seed.cnf
#  ../../../bin/xmaplelcm -dip-pair-min=20 -cpu-lim=1800 kxor-4-$3-$3-xor-3-seed-$seed.cnf

 
 # 	cnfgen randkxor 3 $n $n -T or 3 > s.cnf	 --> FINISHED IN THE CLUSTER (A)
# 	cnfgen randkxor 3 $n $n -T xor 2 > s.cnf --> RUNNIG IN THE CLUSTER (B)
# 	cnfgen randkxor 3 $n $n -T or 2 > s.cnf	--> RUNNIG IN THE CLUSTER (C)
# 	cnfgen randkxor 3 $n $n -T xor 3 > s.cnf	
#       cnfgen randkxor 4 n n -T or 3
#       cnfgen randkxor 4 n n -T xor 3"


