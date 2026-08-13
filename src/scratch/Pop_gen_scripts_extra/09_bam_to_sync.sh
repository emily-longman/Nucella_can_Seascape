#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=make_sync

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Request CPU
#SBATCH --cpus-per-task=6

# Submit job array
#SBATCH --array=1-19

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Bam_sync.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will convert the bam files to gsync files (genome-wide extension of the SYNC file format) and will mask regions based on minimum and max depth thresholds. 

# Load modules 
module load python3.10-anaconda/2023.03-1
module load openjdk/1.8.0
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam

### program dependencies
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
gatk3=/netfiles/nunezlab/Shared_Resources/Software/gatk3/gatk/GenomeAnalysisTK.jar
#tabix=/netfiles/nunezlab/Shared_Resources/Software/htslib/tabix # loaded with samtools
#bgzip=/netfiles/nunezlab/Shared_Resources/Software/htslib/bgzip # loaded with samtools

# Load Mpileup2sync and MaskSYNC_snape_complete
Mpileup2Sync=/netfiles/nunezlab/Shared_Resources/Software/DESTv2/mappingPipeline/scripts/Mpileup2Sync.py

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked.fa

# This is the path to the gff file.
GFF=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/N.can.gff

# This is the path to the pickled genome file.
PICKLED=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/N.canaliculata_assembly.fasta.softmasked_pickled.ref

#--------------------------------------------------------------------------------

# Define parameters
JAVAMEM=59G
THREADS=5
echo ${SLURM_ARRAY_TASK_ID}

# GATK parameters
base_quality_threshold=25
illumina_quality_coding=1.8
minIndel=5

#--------------------------------------------------------------------------------

# Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.
GUIDE_FILE=$WORKING_FOLDER/guide_files/Merge_bams.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
## Population      Merged_name 1      Merged_name 2 
##   ARA	       ARA_S168_L006	  ARA_S13_L008
##   BMR	       BMR_S156_L006	  BMR_S1_L008
##   CBL	       CBL_S169_L006	  CBL_S14_L008
##   ...
##   VD	           VD_S161_L006	      VD_S6_L008

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo ${i}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "syncfiles" ]
then echo "Working syncfiles folder exist"; echo "Let's move on."; date
else echo "Working syncfiles folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/syncfiles; date
fi

# Change directory
cd syncfiles

if [ -d "${i}" ]
then echo "Working ${i} folder exist"; echo "Let's move on."; date
else echo "Working ${i} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/syncfiles/${i}; date
fi

#--------------------------------------------------------------------------------

# Use RealignerTargetCreator to identify and create a target intervals list
#java -jar $gatk3 -T RealignerTargetCreator \
#-nt $THREADS \
#-R $REFERENCE \
#-I $WORKING_FOLDER/RGSM_final_bams/${i}.bam \
#-o $WORKING_FOLDER/syncfiles/${i}/${i}.hologenome.intervals

# Perform local realignment for the target intervals using IndelRealigner 
#java -jar $gatk3 -T IndelRealigner \
#-R $REFERENCE \
#-I $WORKING_FOLDER/RGSM_final_bams/${i}.bam \
#-targetIntervals $WORKING_FOLDER/syncfiles/${i}/${i}.hologenome.intervals \
#-o $WORKING_FOLDER/syncfiles/${i}/${i}.contaminated_realigned.bam

#--------------------------------------------------------------------------------

# Run mpileup step 
echo "Run mpileup step"

samtools mpileup \
$WORKING_FOLDER/syncfiles/${i}/${i}.contaminated_realigned.bam \
-B \
-Q ${base_quality_threshold} \
-f $REFERENCE > $WORKING_FOLDER/syncfiles/${i}/${i}_mpileup.txt

#--------------------------------------------------------------------------------

# Transform Mpileup to Sync
echo "Transform Mpileup to Sync"

python3 $Mpileup2Sync \
--mpileup $WORKING_FOLDER/syncfiles/${i}/${i}_mpileup.txt \
--ref $PICKLED \
--output $WORKING_FOLDER/syncfiles/${i}/${i} \
--base-quality-threshold $base_quality_threshold \
--coding $illumina_quality_coding \
--minIndel $minIndel

#--------------------------------------------------------------------------------

# Housekeeping

# Rename files
mv $WORKING_FOLDER/syncfiles/${i}/${i}_masked.sync.gz $WORKING_FOLDER/syncfiles/${i}/${i}.masked.sync.gz
# Unzip sync files
gunzip $WORKING_FOLDER/syncfiles/${i}/${i}.masked.sync.gz

#--------------------------------------------------------------------------------

# Compress files so compatible with gzip (samtools program)
# This allows indexes to be built against the compressed file and you can still retrieve portions of the data without having to decompress.
bgzip $WORKING_FOLDER/syncfiles/${i}/${i}.masked.sync

# Index tab-deliminated genome position files (samtools program)
tabix -s 1 -b 2 -e 2 $WORKING_FOLDER/syncfiles/${i}/${i}.masked.sync

#--------------------------------------------------------------------------------

date
echo "done"