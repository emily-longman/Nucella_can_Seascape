#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=glms_bio-oracle_rerun1

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=3-00:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Request CPU
#SBATCH --cpus-per-task=8

# Submit job array
#SBATCH --array=4,14,16,17,19,23,27,37,39,44,45,53,56,67,69,74,76,95,120,123,193,216,273,523,526,527,528,529,530,531,533,534,535,537,541,542,543,544,545,548,549,550,551,552,553,554,555,556,558,560,561,562,564,568,569,570,571,572,574,577,578,579,580,582,583,586,590,591,592,594,595,596,598,599,600,601,604,605,610,611,612,614,615,618,619,622,628,629,633,634,635,636,637,640,641,642,645,659,662,664,667,669,670,673,675,676,681,683,684,685,689,690,693,698,700,709,710,712,715,718,721,724,725,726,727,728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,791,792,793,794,795,796,797,798,799,800%50

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 02_glms_bio-oracle.R script. 
# This will run glms and look at the association of the allele frequencies of the outlier SNPs and the Bio-oracle environmental data.
# Note: Several glms required more than the 30 hours originally allotted. Thus, those array IDs were rerun in two sets.

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Script folder.
SCRIPT_FOLDER=$WORKING_FOLDER/src/02_GEA/02_glms

#--------------------------------------------------------------------------------

# Guide file 
guide_file=$WORKING_FOLDER/guide_files/scaffold_names_guide_file_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 18919, 2; 999 partitions each with 19 scaffold names)
# Scaffold name         # Partition/array
# Backbone_10001              1
# Backbone_10003              1
# Backbone_10004              1
# Backbone_10005              1
# ....

#--------------------------------------------------------------------------------

# Determine partition to process 

# Change directory
cd $WORKING_FOLDER/data/processed/GEA/glms

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $guide_file | awk '{print $1}' > scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

# List scaffold names
cat scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/glms

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_chunk_analysis" ]
then echo "Working glms_chunk_analysis folder exist"; echo "Let's move on."; date
else echo "Working glms_chunk_analysis folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_chunk_analysis; date
fi

#--------------------------------------------------------------------------------

# Run R script

Rscript $SCRIPT_FOLDER/02_glms_bio-oracle.R "${SLURM_ARRAY_TASK_ID}"

#--------------------------------------------------------------------------------

rm scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

# Say done
echo "done"