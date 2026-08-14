#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --array=1-19

#--------------------------------------------------------------------------------
# My additions:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "merge_lanes" ]
then echo "Working merge_lanes folder exist"; echo "Let's move on."; date
else echo "Working merge_lanes folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/merge_lanes; date
fi

# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_samtools=https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0

# Guide File
GUIDE_FILE=$WORKING_FOLDER/data/processed/repadapt/guide_files/Merge_bams.txt
# Extract sample names/files
i=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Population i:" ${i}

# Make temporary linefile with list of input BAM files
ls $WORKING_FOLDER/data/processed/repadapt/add_RG/${i}_*_sorted_dedup_RG.bam > ${i}.guide.txt

# Merge the 2 sequencing lanes
apptainer run $repadapt_samtools samtools merge \
-b ${i}.guide.txt \
$WORKING_FOLDER/data/processed/repadapt/merge_lanes/${i}_sorted_dedup_RG_lanes_merged.bam

# Housekeeping: Remove the temporary guide file
rm ${i}.guide.txt