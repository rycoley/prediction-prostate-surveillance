#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/jags/usr/include/JAGS
export LD_LIBRARY_PATH=/home/bst/student/rcoley/jags/usr/lib64
echo $LD_LIBRARY_PATH

cd /home/bst/student/rcoley/inhealth/prediction-model/for-git/simulations/IOP-BX

module load R/3.1

#$ -cwd

#set environment variables

#output files
#$ -o sim-iop-bx.out
#$ -e sim-iop-bx.err

#request memory
#$ -l h_vmem=12G

#request nodes

#send ids
#$ -t 1-200
#$ -V

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH  /home/bst/student/rcoley/inhealth/prediction-model/for-git/simulations/IOP-BX/sim-IOP-BX.R sim-IOP-BX.Rout

#clean up temp files?

