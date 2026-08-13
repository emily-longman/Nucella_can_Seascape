# Revealing the abiotic and biotic drivers of past and future local adaptation in a coastal dogwhelk

## Project Summary

We used population genomic data along with abiotic and biotic datasets to identify the ecological drivers most strongly associated with patterns of adaptive genomic variation in the channeled dogwhelk, *Nucella canaliculata*. Our ecological datasets include climatic data on oceanic conditions ([Bio-Oracle](https://www.bio-oracle.org/)), biotic data on the abundance of N. canaliculata prey, competitors, and predators from a large-scale monitoring program on Pacific Coast rocky shores (i.e., [MARINe](https://marine.ucsc.edu/): Multi-Agency Rocky Intertidal Network), and a new dataset on the shell traits of the rocky shore foundation mussel species, *M. californianus*. Leveraging these signals of local adaptation, we used genomic analyses and population genetic simulations to predict the vulnerability of populations to ocean acidification. 

The bioinformatics pipeline was completed on the Vermont Advanced Computing Center ([VACC](https://www.uvm.edu/vacc)).

### Research Questions

1) What are the abiotic and biotic drivers of geographic patterns of adaptive genomic variation in a low dispersing marine species?
2) What are the genomic bases underlying patterns of local adaptation to varying selection pressures? 
3) Will future environmental change erode existing patterns of local adaptation, or can populations maintain adaptation through evolutionary responses? 

## File Structure

The files in this project are organized in the following structure. All files for arrays and metadata are in "guide_files". All subsequent directories will be generated therein. The VCF for the pool-seq data can be found on [Zenodo](https://zenodo.org/records/18623551).
 - data/
     - raw/
     - processed/
 - guide_files
 - src/
    - 01_ecological_data/
    - 02_demography/
    - 03_GEA_glms/
    - 04_baypass/
 - output/
     - figures/
     - tables/

## Part 1 - Assemble the Draft Genome

### 01 - Prepare the Raw Data 