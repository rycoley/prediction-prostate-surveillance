#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/jags/usr/include/JAGS
export LD_LIBRARY_PATH=/home/bst/student/rcoley/jags/usr/lib64
echo $LD_LIBRARY_PATH

cd /home/bst/student/rcoley/inhealth/prediction-model/for-git/simulations/IOP-SURG


#$ -cwd

#set environment variables

#output files
#$ -o is-iop-surg.out
#$ -e is-iop-surg.err

#request memory
#$ -l h_vmem=5G

#request nodes


# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH  /home/bst/student/rcoley/inhealth/prediction-model/for-git/simulations/IOP-SURG/IS-for-eta-known-IOP-SURG.R IS-IOP-SURG.Rout

#clean up temp files?

