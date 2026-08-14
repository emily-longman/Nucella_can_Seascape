#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=50G

#--------------------------------------------------------------------------------

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape


###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF BCFTOOLS and VCFTOOLS IN YOUR MACHINE/SERVER  = bcftools v. 1.16   vcftools 0.1.16
# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_vcftools=https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321hdcf5f25_9
repadapt_bcftools=https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1
repadapt_samtools=https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0
#module load StdEnv/2020 intel/2020.1.217 tabix/0.2.6

########### Only filtering the VCF to exclude sites with QUAL < 30 and invariant ALT/ALT sites (AC = AN) ###########


apptainer run $repadapt_bcftools bcftools filter -e 'AC=AN || MQ < 30' $WORKING_FOLDER/data/processed/repadapt/vcf/ncanaliculata.vcf.gz -Ov > $WORKING_FOLDER/data/processed/repadapt/vcf/ncanaliculata_filtered.vcf

### As an alternative to the bcftools command above, vcftools can also be used to filter by QUAL:
### vcftools --gzvcf bplaty.vcf.gz --minQ 30 --recode --recode-INFO-all --stdout > bplaty_filtered.vcf

# Change directory
cd $WORKING_FOLDER/data/processed/repadapt/vcf

apptainer run $repadapt_samtools bgzip ncanaliculata_filtered.vcf
apptainer run $repadapt_samtools tabix -p vcf ncanaliculata_filtered.vcf.gz