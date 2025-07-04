#!/usr/bin/bash

#SBATCH --job-name=bam_to_bw         ## Name of the job for the scheduler
#SBATCH --account=siwase99         ## Your PI's uniqname plus 0,99, or other number
#SBATCH --partition=largemem         ## Choose: standard, largemem, gpu, spgpu, debug            
##SBATCH --gpus=1  			   ## if partition=gpu, number of GPUS needed
#SBATCH --nodes=1                      ## number of nodes you are requesting
#SBATCH --ntasks=1                     ## how many task spaces do you want to reserve
#SBATCH --cpus-per-task=16            ## how many cores do you want to use per task
#SBATCH --time=24:00:00                 ## Maximum length of time you are reserving the 
#SBATCH --mem=240G                      ## Memory limit  
#SBATCH --mail-user=siwase@umich.edu  ## send email notifications to umich email listed
#SBATCH --mail-type=ALL                ## when to send email (standard values are:
                                       ## NONE,BEGIN,END,FAIL,REQUEUE,ALL.
   ## (See documentation for others)
#SBATCH --output=./%x-%j               ## send output and error info to the file listed
                                       ##(optional: different name format than default) 

# activate conda enviroment to use snakemakes
source ~/miniconda3/etc/profile.d/conda.sh # change if needed
conda activate /home/siwase/miniconda3/envs/snakemake # change


# load Bioinformatics module
module load Bioinformatics
module load star


# run snakemake to map the bowtie2 output to bigwig
snakemake --cores 16 -s code/bam_to_bigwig.smk --rerun-incomplete 