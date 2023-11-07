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

for fo in 0.00
	do

for CI in 1 2 
	do

# ==============================================================================
# =                               Rodando simulação                            = 
#===============================================================================

./gota.out -L ${L} -R ${R} -a ${a} -h ${h} -w ${w} -fo ${fo} -CI ${CI} -s  12345567 

done
done
done
done
done
done
done



# ==============================================================================

# ==============================================================================

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

for fo in 0.00
	do

for CI in 1 2 
	do

# ==============================================================================
# =                        Preparando script imagem final                      = 
#===============================================================================
cat << EOF > gota.gp
set term post eps enhan color font ",30" 
set out "gota.eps"
unset key

LL=${L} 
b=${a}   
h=${h} 
w=${w}   
BASE=3

EOF

cat teste.gp >> gota.gp

# ==============================================================================
# =                            Movendo os aqrquivos CB                         = 
#===============================================================================
if [ ${CI} -eq 1 ] 
then

#==============================================================================
#=                            Fazendo a imagem final                          = 
#===============================================================================

cp CB_gota_3d_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}_fo_${fo}_LAST.dsf file.dsf
awk -v file='file.dsf' -f separa_conf.awk file.dsf
gnuplot gota.gp

rm gota.gp
rm file.dsf

# ==============================================================================
# =                             Fazendo os graficos                            = 
#===============================================================================

cp CB_gota_3d_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}_fo_${fo}.dsf file.dsf

python3 graficos.py

mv volume_total.eps CB_volume_total_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv energia.eps CB_energia_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv B.eps CB_B_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv raio.eps CB_R_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv theta.eps CB_theta_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv v_baixo.eps CB_vpil_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv v_per.eps CB_vper_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps

fi

# ==============================================================================
# =                           Movendo os arquivos WE                           = 
#===============================================================================
if [ ${CI} -eq 2 ] 
then

# ==============================================================================
# =                            Fazendo a imagem final                          = 
#===============================================================================

cp WE_gota_3d_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}_fo_${fo}_LAST.dsf file.dsf
awk -v file='file.dsf' -f separa_conf.awk file.dsf
gnuplot gota.gp 

rm gota.gp
rm file.dsf

# ==============================================================================
# =                             Fazendo os graficos                            = 
#===============================================================================

cp WE_gota_3d_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}_fo_${fo}.dsf file.dsf

python graficos.py

mv volume_total.eps WE_volume_total_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv energia.eps WE_energia_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv B.eps WE_B_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv raio.eps WE_R_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv theta.eps WE_theta_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv v_baixo.eps WE_vpil_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps
mv v_per.eps WE_vper_L_${L}_R_${R}_a_${a}_h_${h}_w_${w}.eps

fi 
#===============================================================================


done
done
done
done
done
done
done
