# Interactive sesssion to transform genome


#--------------------------------------------------------------------------------

# Load modules 
module load python3.10-anaconda/2023.03-1

pickle=/netfiles/nunezlab/Shared_Resources/Software/DESTv2/mappingPipeline/scripts/PickleRef.py

#--------------------------------------------------------------------------------

# Define important file locations

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked

#--------------------------------------------------------------------------------

# Change directory
cd /netfiles/pespenilab_share/Nucella/processed/Base_Genome

# Pickle genome
python3 $pickle \
--ref $REFERENCE \
--output N.canaliculata_assembly.fasta.softmasked_pickled > N.canaliculata_assembly.fasta.softmasked_pickled_ref.out
