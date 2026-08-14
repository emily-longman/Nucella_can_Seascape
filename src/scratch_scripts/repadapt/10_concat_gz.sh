#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=200G

#--------------------------------------------------------------------------------

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "vcf" ]
then echo "Working vcf folder exist"; echo "Let's move on."; date
else echo "Working vcf folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/vcf; date
fi

###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF BCFTOOLS IN YOUR MACHINE/SERVER  = bcftools v. 1.16
# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_bcftools=https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1
repadapt_samtools=https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0

### Concatenate all the chromsome vcfs produced in script 09. 
### list.txt is a list of the 14 (in this case) vcfs produced in script 09.
ls $WORKING_FOLDER/data/processed/repadapt/mpileup_gz/*vcf.gz > $WORKING_FOLDER/data/processed/repadapt/vcf.list.txt
#ls $WORKING_FOLDER/data/processed/repadapt/mpileup/*vcf.gz > $WORKING_FOLDER/data/processed/repadapt/vcf.list.txt

### Here we concatenate them in a single vcf

apptainer run $repadapt_bcftools bcftools concat -f $WORKING_FOLDER/data/processed/repadapt/vcf.list.txt -Oz > $WORKING_FOLDER/data/processed/repadapt/vcf/ncanaliculata.vcf.gz
cd $WORKING_FOLDER/data/processed/repadapt/vcf
apptainer run $repadapt_samtools tabix -p vcf ncanaliculata.vcf.gz