#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/jags/usr/include/JAGS
export LD_LIBRARY_PATH=/home/bst/student/rcoley/jags/usr/lib64
echo $LD_LIBRARY_PATH

cd /home/bst/student/rcoley/inhealth/prediction-model/for-git

module load R/3.1

#$ -cwd

#set environment variables

#output files
#$ -o unadj-pred-mod.out
#$ -e unadj-pred-mod.err

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
R CMD BATCH  /home/bst/student/rcoley/inhealth/prediction-model/for-git/jhas-analysis/UNADJ/UNADJ-call-jags.R unadj-pred-mod.Rout

#clean up temp files?

