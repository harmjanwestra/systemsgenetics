#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH -J JOBNAME
#SBATCH -o OUTLOGFILE
#SBATCH -e ERRLOGFILE

# change above to match your cluster environment
# the JOBNAME, OUTLOGFILE and ERRLOGFILE variables will be replaced by the software

# this may not be needed on all clusters
set -e
set -u

# change with module load if module system available on your HPC
PATH=$PATH:/path/to/bwaBinary/

# this string below is replaced with the required commands
CMD