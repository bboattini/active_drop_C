#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=10GB
#SBATCH --output=slurm_%j.out
#SBATCH --job-name=BD_water_oil

module purge
# Comandos de execução do seu programa:
L=$1
R=$2
a=$3
h=$4
w=$5
fo=$6
CI=$7
./gota.out -L ${L} -R ${R} -a ${a} -h ${h} -w ${w} -fo ${fo} -CI ${CI} -s  12345567
