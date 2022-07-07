#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Author Nick Wong, Monash Bioinformatics Platform
#     nick.wong@monash.edu

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=bcl-convert


# To set a project account for credit charging, 
#SBATCH --account=tn26


# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --profile=task

# Memory usage (MB)
#SBATCH --mem-per-cpu=6000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-4:00:00
#SBATCH --partition=genomics
#SBATCH --qos=genomics

# To receive an email when job completes or fails
#SBATCH --mail-user=nick.wong@monash.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Command to run a serial job
module load singularity


BCLCONVERT=/home/nwong/bin/bcl-convert
OUTPUT=/home/nwong/tn26_scratch/gemma.hartley/results
mkdir -p ${OUTPUT}/logs

singularity exec --bind ${OUTPUT}:/var/log/bcl-convert --bind ${OUTPUT}:/results \
  --bind /home/nwong/tn26_scratch/gemma.hartley/2022-01-28-runFiles/20220125_GH_LB_nvSPXPv15_sideB/Files:/input  ${BCLCONVERT}/bcl-convert-3.6.3.sif bcl-convert \
  --sample-sheet /results/newSheet.csv \
  -f \
  --output-directory /results \
  --bcl-input-directory /input
