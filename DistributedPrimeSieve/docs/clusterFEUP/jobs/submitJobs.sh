#!/bin/bash
qsub -l nodes=1:ppn=1 -N DPS_Alg09_28bits -v RANGE_IN_BITS=28 DPS_Alg09.sh
qsub -l nodes=1:ppn=1 -N DPS_Alg09_32bits -v RANGE_IN_BITS=32 DPS_Alg09.sh
qsub -l nodes=1:ppn=1 -N DPS_Alg09_36bits -v RANGE_IN_BITS=36 DPS_Alg09.sh
qsub -l nodes=1:ppn=1 -N DPS_Alg09_40bits -v RANGE_IN_BITS=40 DPS_Alg09.sh
qsub -l nodes=1:ppn=1 -N DPS_Alg09_44bits -v RANGE_IN_BITS=44 DPS_Alg09.sh
qsub -l nodes=1:ppn=1 -N DPS_Alg09_48bits -v RANGE_IN_BITS=48 DPS_Alg09.sh

qsub -l nodes=1:ppn=16 -N DPS_Alg14_28bits -v RANGE_IN_BITS=28 DPS_Alg14.sh
qsub -l nodes=1:ppn=16 -N DPS_Alg14_32bits -v RANGE_IN_BITS=32 DPS_Alg14.sh
qsub -l nodes=1:ppn=16 -N DPS_Alg14_36bits -v RANGE_IN_BITS=36 DPS_Alg14.sh
qsub -l nodes=1:ppn=16 -N DPS_Alg14_40bits -v RANGE_IN_BITS=40 DPS_Alg14.sh
qsub -l nodes=1:ppn=16 -N DPS_Alg14_44bits -v RANGE_IN_BITS=44 DPS_Alg14.sh
qsub -l nodes=1:ppn=16 -N DPS_Alg14_48bits -v RANGE_IN_BITS=48 DPS_Alg14.sh
qsub -l nodes=1:ppn=16 -N DPS_Alg14_64bits -v RANGE_IN_BITS=64 DPS_Alg14.sh

qsub -l nodes=16:ppn=16 -N DPS_Alg19_28bits -v RANGE_IN_BITS=28,NUMBER_PROCESSES=33,NUMBER_SLOTS_PER_NODE=2 DPS_Alg19_.sh
qsub -l nodes=16:ppn=16 -N DPS_Alg19_32bits -v RANGE_IN_BITS=32,NUMBER_PROCESSES=33,NUMBER_SLOTS_PER_NODE=2 DPS_Alg19_.sh
qsub -l nodes=16:ppn=16 -N DPS_Alg19_36bits -v RANGE_IN_BITS=36,NUMBER_PROCESSES=33,NUMBER_SLOTS_PER_NODE=2 DPS_Alg19_.sh
qsub -l nodes=16:ppn=16 -N DPS_Alg19_40bits -v RANGE_IN_BITS=40,NUMBER_PROCESSES=33,NUMBER_SLOTS_PER_NODE=2 DPS_Alg19_.sh
qsub -l nodes=16:ppn=16 -N DPS_Alg19_44bits -v RANGE_IN_BITS=44,NUMBER_PROCESSES=33,NUMBER_SLOTS_PER_NODE=2 DPS_Alg19_.sh
qsub -l nodes=16:ppn=16 -N DPS_Alg19_48bits -v RANGE_IN_BITS=48,NUMBER_PROCESSES=33,NUMBER_SLOTS_PER_NODE=2 DPS_Alg19_.sh
qsub -l nodes=16:ppn=16 -N DPS_Alg19_64bits -v RANGE_IN_BITS=64,NUMBER_PROCESSES=33,NUMBER_SLOTS_PER_NODE=2 DPS_Alg19_.sh
