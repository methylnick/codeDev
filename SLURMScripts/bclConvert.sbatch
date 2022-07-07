#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Author Nick Wong, Monash Bioinformatics Platform
#     nick.wong@monash.edu

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=bcl-convert


# To set a project account for credit charging, 
#SBATCH --account=ls25


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


=/home/nick-wong/bin/bcl-convert
OUTPUT=/home/nick-wong/projects/david.curtis/2021-06-21-ASPREE-SP-2.4-3.1/fastqs
mkdir -p ${OUTPUT}/logs

singularity exec --bind ${OUTPUT}/logs:/var/log/bcl-convert ${}/bcl-convert-3.6.3.sif bcl-convert \
  --sample-sheet /home/nick-wong/projects/david.curtis/2021-06-21-ASPREE-SP-2.4-3.1/20210618_AL_ASPREE_rhampseq_nvSPXPv15_sideB_120pM/Files/SampleSheet.csv \
  -f \
  --output-directory ${OUTPUT} \
  --bcl-input-directory /home/nick-wong/projects/david.curtis/2021-06-21-ASPREE-SP-2.4-3.1/20210618_AL_ASPREE_rhampseq_nvSPXPv15_sideB_120pM/Files
