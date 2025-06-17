#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Blast_blastp_array_run2

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G

# Request CPU
#SBATCH --cpus-per-task=5

# Submit job array
#SBATCH --array=1-934%35

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out  # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run blastp on the N. canaliculata protein file.

#Load modules
ncbi=/gpfs1/home/e/l/elongman/software/ncbi-blast-2.16.0+/bin

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# PROTEIN_FILE is the N. canaliculata protein file generated from the gff
PROTEIN_FILE=$WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/Nucella_canaliculata_protein.fa

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology

if [ -d "blastp_array" ]
then echo "Working blastp_array folder exist"; echo "Let's move on."; date
else echo "Working blastp_array folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array

if [ -d "Run_2" ]
then echo "Working Run_2 folder exist"; echo "Let's move on."; date
else echo "Working Run_2 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2

if [ -d "Partition_${SLURM_ARRAY_TASK_ID}" ]
then echo "Working Partition_${SLURM_ARRAY_TASK_ID} folder exist"; echo "Let's move on."; date
else echo "Working Partition_${SLURM_ARRAY_TASK_ID} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}; date
fi

#--------------------------------------------------------------------------------

## Import master partition file 
guide_file=$WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/guide_files/protein_file_names_array_run_2.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 40634, 5; 903 partitions)
# Protein header                                   # Partition
#>g19552.t1 gene=g19552 seq_id=ntLink_0 type=cds         1
#>g19553.t1 gene=g19553 seq_id=ntLink_0 type=cds         1
#>g19554.t1 gene=g19554 seq_id=ntLink_0 type=cds         1

#--------------------------------------------------------------------------------

# Determine partition to process 

# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the headers of the protein file associated with a given partition/slurm arrray task ID
awk '$5=='${SLURM_ARRAY_TASK_ID}'' $guide_file | awk '{print $1 " " $2 " " $3 " " $4}' > protein.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}

# Cat file of protein headers and start while loop
cat $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/protein.names.${SLURM_ARRAY_TASK_ID}.txt | \
while read protein 
do 
echo ${protein}

# Extract protein name - give it the full header and extract just the protein name 
i=$(echo ${protein} | awk -e '{print($1)}' | sed 's/>//g')
echo ${i}

# Break up the protein file into each protein
grep -EA 1 "${protein}" ${PROTEIN_FILE} > \
$WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}/protein.${i}.partition.${SLURM_ARRAY_TASK_ID}.fa

# Use the blastp command to compare the Nucella protein file with the uniprot database
$ncbi/blastp -query $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}/protein.${i}.partition.${SLURM_ARRAY_TASK_ID}.fa \
-db $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/uniprot/uniref90 \
-out $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}/blastp_vs_uniref90.outfmt6_protein.${i} \
-outfmt 6 \
-evalue 1e-5 \
-max_target_seqs 1

# query <File_In>: Input file name
# db: BLAST database name
# outfmt 6: alignment view options: Tabular
# evalue: Expectation value (E) threshold for saving hits. Default = 10 
# max_target_seqs: Maximum number of aligned sequences to keep 

done

#--------------------------------------------------------------------------------

# Merge results into one txt file
cat $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}/blastp_vs_uniref90.outfmt6_protein.* \
> $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/blastp_vs_uniref90.outfmt6_partition_${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Housekeeping - remove intermediate files
rm $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}/protein.*.partition.${SLURM_ARRAY_TASK_ID}.fa
rm $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}/blastp_vs_uniref90.outfmt6_protein.*
rm $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/protein.names.${SLURM_ARRAY_TASK_ID}.txt
rm -R $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/blastp_array/Run_2/Partition_${SLURM_ARRAY_TASK_ID}
