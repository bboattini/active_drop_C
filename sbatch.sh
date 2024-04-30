#!/bin/bash 
#SBATCH -n 1 # Numero de CPU cores a serem alocados 
##SBATCH -N 1 # Numero de nodes a serem alocados
##SBATCH -t 0-00:01 # Tempo limite de execucao (D-HH:MM)
#SBATCH -p long # Particao (fila) a ser submetido
#SBATCH --qos qos_long # QOS 
##SBATCH -w node106
# Comandos de execução do seu programa:
L=$1
R=$2
a=$3
h=$4
w=$5
fo=$6
CI=$7
./gota.out -L ${L} -R ${R} -a ${a} -h ${h} -w ${w} -fo ${fo} -CI ${CI} -s  12345567
