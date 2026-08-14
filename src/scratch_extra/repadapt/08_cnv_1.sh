#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=20G

#--------------------------------------------------------------------------------
# My additions/modifications:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "cnv" ]
then echo "Working cnv folder exist"; echo "Let's move on."; date
else echo "Working cnv folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/cnv; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/repadapt/genome

# make a genomefile for bedtools
awk '{print $1"\t"$2}' $WORKING_FOLDER/data/processed/repadapt/genome/*.fai > $WORKING_FOLDER/data/processed/repadapt/cnv/genome.bed

# make a BED file of 5000 bp windows from FASTA index
awk -v w=5000 '{chr = $1; chr_len = $2;
    for (start = 0; start < chr_len; start += w) {
        end = ((start + w) < chr_len ? (start + w) : chr_len);
        print chr "\t" start "\t" end;
    }
}' $WORKING_FOLDER/data/processed/repadapt/genome/*.fai > $WORKING_FOLDER/data/processed/repadapt/cnv/windows.bed

# now make location list from window bedfile and sort for join
awk -F "\t" '{print $1":"$2"-"$3}' $WORKING_FOLDER/data/processed/repadapt/cnv/windows.bed | sort -k1,1 > $WORKING_FOLDER/data/processed/repadapt/cnv/windows.list

# and a bed file of each gene (Note: must move gff file to genome directory)
awk '$3 == "gene" {print $1"\t"$4"\t"$5}' $WORKING_FOLDER/data/processed/repadapt/genome/*gff | uniq > $WORKING_FOLDER/data/processed/repadapt/cnv/genes.bed

# we also sort this file based on the order of the reference index
cut -f1 *.fai | while read chr; do awk -v chr=$chr '$1 == chr {print $0}' $WORKING_FOLDER/data/processed/repadapt/cnv/genes.bed | sort -k2,2n; done > $WORKING_FOLDER/data/processed/repadapt/cnv/genes.sorted.bed
mv $WORKING_FOLDER/data/processed/repadapt/cnv/genes.sorted.bed $WORKING_FOLDER/data/processed/repadapt/cnv/genes.bed

# now make location list from sorted gene bedfile and sort for join
awk -F "\t" '{print $1":"$2"-"$3}' $WORKING_FOLDER/data/processed/repadapt/cnv/genes.bed | sort -k1,1 > $WORKING_FOLDER/data/processed/repadapt/cnv/genes.list
