#!/bin/bash
#SBATCH -N 1
#SBATCH -t 4:00:00
#SBATCH -A aosc-hi
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

time ./maooam
