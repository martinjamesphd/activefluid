#!/bin/bash
#$ -V
#$ -N mj_ibm
#$ -cwd
#$ -o /scratch07/james/test/out.dat
#$ -e /scratch07/james/test/err.dat
#$ -pe openmp-fulla01 24
echo "got $NSLOTS slots."
echo "Start time is `date`"
./ns2d input_flow.dat 0
echo "End time is `date`"
exit 0
