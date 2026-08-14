#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=cat_short_reads

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# WORKING_FOLDER is the core folder where this pipeline is being run.
TMP=/gpfs3tmp/pi/jcnunez

# Make folders
cd $TMP
if [ -d "All_shortreads" ]
then echo "Working All_shortreads folder exist"; echo "Let's move on."; date
else echo "Working All_shortreads folder doesnt exist. Let's fix that."; mkdir $TMP/All_shortreads; date
fi

#--------------------------------------------------------------------------------

# Change directory
cd $TMP/All_shortreads

# Copy reads from lane L006
scp $TMP/Population_genomics/Un_DTSA1002/Project_ELEL_Nova1179P_Longman/*fastq.gz .

# Copy reads from lane L008
scp $TMP/Population_genomics/Un_DTSA1030/Project_ELEL_Nova1179P_Longman/*fastq.gz .