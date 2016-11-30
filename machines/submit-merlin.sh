#!/bin/bash
# xhpl_pe_openmpi_using_openmpi-1.5.4-gcc-4.6.1-mkl-10.sge
#
# This script uses the parallel environment (PE) "openmpi" with an explicit machinefile.
# This script must be used with qsub command - do NOT run it as a stand-alone
# script unless NSLOTS and TMPDIR/machines are properly set for the MPI command MPICMD
# and PE_HOSTFILE is set to an empty file.

# Define your job name, parallel environment with the number of slots, and run time:
#$ -cwd
#$ -N Structure-Factor
#$ -pe openmpi 2
#$ -o script.out 
#$ -e script.err
#$ -m ea 
#$ -M dominik.schildknecht@psi.ch
#$ -l s_rt=00:05:00
#$ -l h_rt=00:05:30

###################################################
# Fix the SGE environment-handling bug (bash):
source /usr/share/Modules/init/sh
export -n -f module

# No modules are loaded explicitly - set the environment below.
###################################################
# Set the environment variables:
#OPENMPI=/opt/mpi/openmpi-1.5.4-gcc-4.6.1
GCC_DIR=/opt/gcc/gcc-4.6.1
MKL_DIR=/opt/intel/intel-12.0/composerxe-2011.3.174/mkl
MPIEXEC=/opt/psi/Compiler/openmpi/1.8.4/gcc/4.7.4/bin/mpirun
LD_LIBRARY_PATH=$OPENMPI/lib:$GCC_DIR/lib64
LD_LIBRARY_PATH=$MKL_DIR/lib/intel64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
export OMP_NUM_THREADS=1

source /gpfs/home/schildknecht_d/.bashrc

##############
# BEGIN DEBUG
# Print the SGE environment on master host:
echo "================================================================"
echo "=== SGE job  JOB_NAME=$JOB_NAME  JOB_ID=$JOB_ID"
echo "================================================================"
echo date=`date`
echo hostname=`hostname`
echo working directory=`pwd`
echo "NSLOTS=$NSLOTS"
echo "PE_HOSTFILE=$PE_HOSTFILE"
#cat $PE_HOSTFILE
#echo "Machinefile created by openmpi PE $TMPDIR/machines:"
#cat $TMPDIR/machines
#echo "================================================================"
#echo "Running environment:"
#env
#echo "================================================================"
#echo "Loaded environment modules:"
#module list 2>&1
#echo
# END DEBUG
##############

###################################################
# The command to run with mpiexec:
#CMD=$HOME/bin/mc++
CMD="/gpfs/home/schildknecht_d/src/mcpp/bin/mc++ --mpi"
ARGS='parm.in.xml'

###############
## BEGIN DEBUG
## Check that the libraries are available (on the master host):
#echo "ldd $CMD"
#ldd $CMD
#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
## Check the number of threads used by OpenMP:
#echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
## END DEBUG
###############

# The MPI command to run:
MPICMD="$MPIEXEC -x LD_LIBRARY_PATH -x OMP_NUM_THREADS -np $NSLOTS -machinefile $TMPDIR/machines $CMD $ARGS"
#MPICMD="$MPIEXEC --prefix $OPENMPI -x LD_LIBRARY_PATH -x OMP_NUM_THREADS -np $NSLOTS -machinefile $TMPDIR/machines $CMD $ARGS"
#echo "Command to run:"
#echo "$CMD $ARGS"
#echo
#echo "MPI Command to run:"
#echo "$MPICMD"

$MPICMD

################################################################################
