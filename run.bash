#!/bin/bash  

gcc -Wall -O3 gota.c -lm -o gota.out 

for L in 240
	do

for R in 50
	do

for a in  3 8
	do
	
for h in 8
	do
	
for w in 5
	do

for fo in 0.1 1 10
	do

for CI in 2
	do
# Create the directory
dir_name="dados_${CI}_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}_fo_${fo}"
mkdir -p $dir_name

# Move the executable to the directory
cp gota.out $dir_name

# Change to the directory
cd $dir_name


# ==============================================================================
# =                               Rodando simulação                            = 
#===============================================================================

./gota.out -L ${L} -R ${R} -a ${a} -h ${h} -w ${w} -fo ${fo} -CI ${CI} -s  12345567
#sbatch ./sbatch.sh -L ${L} -R ${R} -a ${a} -h ${h} -w ${w} -fo ${fo} -CI ${CI}

# Change back to the original directory
cd -

done
done
done
done
done
done
done



# ==============================================================================

# ==============================================================================

