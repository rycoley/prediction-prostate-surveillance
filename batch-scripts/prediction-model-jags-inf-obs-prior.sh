#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/jags/usr/include/JAGS
export LD_LIBRARY_PATH=/home/bst/student/rcoley/jags/usr/lib64
echo $LD_LIBRARY_PATH

cd /home/bst/student/rcoley/inhealth/prediction-model/inf-obs

module load R/3.1

#$ -cwd

#set environment variables

#output files
#$ -o pred-mod-prior.out
#$ -e pred-mod-prior.err

#request memory
#$ -l h_vmem=10G

#request nodes

#send ids
#$ -t 1-5
#$ -V

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH  /home/bst/student/rcoley/inhealth/prediction-model/inf-obs/prediction-model-jags-inf-obs-prior.R pred-mod-prior.Rout

#clean up temp files?

