#!/bin/bash

#SBATCH -A m2404
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --qos regular

EXE="/pscratch/sd/p/psivaram/spack/spack-install/linux-amzn2-x86_64_v3/gcc-7.3.1/salmon-1.10.2-6fawyjo5l35b3yvf62goig7c4qgks67x/bin/salmon"
SRUN="srun -n1 -c 64 --cpu-bind=cores"
SALMON="${EXE} quant -i sapien_index -l A -p 64"

for ((i = 1; i <= 20; i++))
do
  INPUT1="data/${i}_1.fq.gz"
  INPUT2="data/${i}_2.fq.gz"
  OUTPUT_VBEM="quants_vbem/sample_${i}"
  OUTPUT_EM="quants_em/sample_${i}"

  $SRUN $SALMON -1 $INPUT1 -2 $INPUT2 -o $OUTPUT_VBEM
  $SRUN $SALMON -1 $INPUT1 -2 $INPUT2 --useEM -o $OUTPUT_EM
done
