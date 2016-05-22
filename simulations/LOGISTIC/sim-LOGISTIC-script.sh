#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/jags/usr/include/JAGS
export LD_LIBRARY_PATH=/home/bst/student/rcoley/jags/usr/lib64
echo $LD_LIBRARY_PATH

cd /home/bst/student/rcoley/inhealth/prediction-model/for-git/simulations/LOGISTIC


#$ -cwd

#set environment variables

#output files
#$ -o sim-logistic.out
#$ -e sim-logistic.err

#request memory
#$ -l h_vmem=5G

#request nodes


# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH  /home/bst/student/rcoley/inhealth/prediction-model/for-git/simulations/LOGISTIC/sim-LOGISTIC.R sim-logistic.Rout

#clean up temp files?

