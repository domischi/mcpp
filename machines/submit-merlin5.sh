#!/bin/bash
#SBATCH --partition=merlin                      # name of slurm partition to submit
#SBATCH --time=23:55:00                         # time limit 
#SBATCH --ntasks=20                             # number of tasks
#SBATCH --nodes=2                               # number of nodes 
#SBATCH --ntasks-per-node=10                    # how to distribute the tasks on the different nodes
#SBATCH --cores-per-socket=8                    # How many cores each processor should have (this prevents execution on merlinc[<=16], which ends with illegal instruction error
#SBATCH --job-name="mc++"                       # Job name
#SBATCH --output=script.out                     # Output File
#SBATCH --open-mode=append                      # The script.out should not be overwritten 
#SBATCH --mail-type=ALL                         # When to write mails
#SBATCH --mail-user=dominik.schildknecht@psi.ch # Where to write mails

CMD="mc++ --mpi "
ARGS="parm.in.xml"
FULL_CMD="mpirun -np $SLURM_NTASKS $CMD $ARGS"

# ===================START==================== 
echo -e "========================================"
echo -e "Time and date:\n\t$(date)\n"
echo -e "Current working directory:\n\t$(pwd)\n"
echo -e "Job ID:\n\t$SLURM_JOB_ID\n"
echo -e "Number of Nodes:\n\t$SLURM_JOB_NUM_NODES\n"
echo -e "On nodes:\n\t$SLURM_JOB_NODELIST\n"
echo -e "Number of Tasks:\n\t$SLURM_NTASKS\n"
echo -e "CPUs on Node:\n\t$SLURM_CPUS_ON_NODE\n"
echo -e "Module List:\n\t$(module list 2>&1)\n"
echo -e "Host Name:\n\t$(mpirun hostname)\n"
echo -e "Command to execute:\n\t$FULL_CMD\n"
echo -e "========================================"
echo -e "\n\n"
echo -e "Run alpspython:"
alpspython single-spin-flip.py
echo -e "\nRun command:"
$FULL_CMD
# ====================END===================== 
