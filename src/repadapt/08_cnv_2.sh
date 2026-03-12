#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1 #1-19

#--------------------------------------------------------------------------------
# My additions/modifications:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#### calculate depth statistics for 19 samples -- see array number above and change accordingly ####

# Load needed modules
# Use apptainer to use the correct versions of programs (samtools v.1.16.1   BEDtools v.2.27.1)
module load apptainer/1.3.4
repadapt_samtools=https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0
repadapt_bedtools=https://depot.galaxyproject.org/singularity/bedtools:2.27.1--0

### Keep the lists below with the same order
#INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list1.txt)  ### list of realigned bam files
#OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt)  ### list of output names (just extract from input names removing bam suffix)

# Guide File
GUIDE_FILE=$WORKING_FOLDER/data/processed/repadapt/guide_files/Merge_bams.txt
# Extract sample names/files
i=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Population i:" ${i}

# dump depth of coverage at every position in the genome
apptainer run $repadapt_samtools samtools depth -aa $WORKING_FOLDER/data/processed/repadapt/realign_indel/${i}_sorted_dedup_RG_lanes_merged_realigned.bam > $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned.depth

# gene depth analysis
echo \n">>> Computing depth of each gene for ${i} <<<"\n
awk '{print $1"\t"$2"\t"$2"\t"$3}' $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned.depth | apptainer run $repadapt_bedtools bedtools map -a $WORKING_FOLDER/data/processed/repadapt/cnv/genes.bed -b stdin -c 4 -o mean -null 0 -g $WORKING_FOLDER/data/processed/repadapt/cnv/genome.bed | awk -F "\t" '{print $1":"$2"-"$3"\t"$4}' | sort -k1,1 > $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned-genes.tsv

# sort gene depth results based on input bed file
join -a 1 -e 0 -o '1.1 2.2' -t $'\t' genes.list $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned-genes.tsv > $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned-genes.sorted.tsv

# window depth analysis
echo \n">>> Computing depth of each window for ${i} <<<"\n
awk '{print $1"\t"$2"\t"$2"\t"$3}' $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned.depth | apptainer run $repadapt_bedtools bedtools map -a $WORKING_FOLDER/data/processed/repadapt/cnv/windows.bed -b stdin -c 4 -o mean -null 0 -g $WORKING_FOLDER/data/processed/repadapt/cnv/genome.bed | awk -F "\t" '{print $1":"$2"-"$3"\t"$4}' | sort -k1,1  > $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned-windows.tsv

# sort window depth results based on input bed file
join -a 1 -e 0 -o '1.1 2.2' -t $'\t' $WORKING_FOLDER/data/processed/repadapt/cnv/windows.list $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned-windows.tsv > $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned-windows.sorted.tsv

# overall genome depth
echo \n">>> Computing depth of whole genome for ${i} <<<"\n
awk '{sum += $3; count++} END {if (count > 0) print sum/count; else print "No data"}' $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned.depth > $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned-wg.txt

echo " >>> Cleaning a bit...
"
rm -rf $WORKING_FOLDER/data/processed/repadapt/cnv/${i}_sorted_dedup_RG_lanes_merged_realigned.depth
echo "
DONE! Check your files"