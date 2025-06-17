#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Blast_format_db

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will format the blast uniprot database

#Load modules
ncbi=/gpfs1/home/e/l/elongman/software/ncbi-blast-2.16.0+/bin

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/outlier_analyses

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "gene_ontology" ]
then echo "Working gene_ontology folder exist"; echo "Let's move on."; date
else echo "Working gene_ontology folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology; date
fi

cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology

if [ -d "uniprot" ]
then echo "Working uniprot folder exist"; echo "Let's move on."; date
else echo "Working uniprot folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/uniprot; date
fi

#--------------------------------------------------------------------------------

# Prepare the Uniprot database
# Note: UniRef90 clusters are generated from the UniRef100 seed sequences with a 90% sequence identity threshold using the MMseqs2 algorithm.

# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/uniprot

# Download and unzip uniprot database in $WORKING_FOLDER/gene_ontology directory
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz 
gunzip -v uniref90.fasta.gz                                                         

# Format uniprot database
$ncbi/makeblastdb -in uniref90.fasta -dbtype prot -out uniref90 -title uniref90

# The uniref90.fasta file contains all UniRef90 entries in FASTA format.
# The format is as follows: UniqueIdentifier ClusterName n=Members Tax=Taxon RepID=RepresentativeMember where:
## UniqueIdentifier is the primary accession number of the UniRef cluster.
## ClusterName is the name of the UniRef cluster.
## Members is the number of UniRef cluster members.
## Taxon is the scientific name of the lowest common taxon shared by all UniRef cluster members.
## RepresentativeMember is the entry name of the representative member of the UniRef cluster.

