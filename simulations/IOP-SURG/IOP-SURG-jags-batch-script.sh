#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/jags/usr/include/JAGS
export LD_LIBRARY_PATH=/home/bst/student/rcoley/jags/usr/lib64
echo $LD_LIBRARY_PATH

cd /home/bst/student/rcoley/inhealth/prediction-model/for-git/simulations/IOP-SURG

module load R/3.1

#$ -cwd

#set environment variables

#output files
#$ -o sim-iop-surg.out
#$ -e sim-iop-surg.err

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
R CMD BATCH  /home/bst/student/rcoley/inhealth/prediction-model/for-git/simulations/IOP-SURG/sim-IOP-SURG.R sim-IOP-SURG.Rout

#clean up temp files?

