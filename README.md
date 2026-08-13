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
    - 04_GEA_baypass/
    - 05_outlier_analyses/
    - 06_genomic_offest/
    - 07_SLiM/
 - output/
     - figures/
     - tables/

## Part 1 - Compile ecological data

To identify which ecological selective forces are driving patterns of adaptive variation in *N. canaliculata*, we used several ecologically abiotic and biotic variables.  

### 01 - Abiotic

We extracted decadal environmental data for the 19 *N. canaliculata* field sites from [Bio-Oracle](https://www.bio-oracle.org/). Our variables of interest were: sea surface temperature (mean, max, min, range, 2010-2019), pH (mean, min, 2010-2018), O2 (mean, 2010-2018), salinity (mean, 2010-2019), and chlorophyll (mean, 2010-2018). For each variable we also extracted future data (2090-2100 under the SSP5-8.5 scenario) The multiple iterations of the three focal scripts were for graphical purposes. 

01_extract_bio-oracle.R - Extract present and future data for variables for large geographic region along the eastern Pacific coastline.

02_format_bio-oracle.R - Filter each variable for the 19 focal *N. canaliculata* field sites.

03_graph_bio-oracle_raster.R - Graph abiotic variables.

### 02 - Biotic

There were two sets of biotic data. The first was abundance/density data of interacting species (prey, competitors, and predators) from [MARINe](https://marine.ucsc.edu/). 

01_summarize_graph_MARINe.R - Filter abundance/density data of interacting species to focal sites and summarize.

02_summarize_Mcalifornianus_shell.R - Filter shell morphology data and summarize.

## Part 2 - Demography and population structure

Filter the Pool-Seq dataset from [Longman et al. (2026)](https://royalsocietypublishing.org/rspb/article/293/2070/20253148/481573/Geographic-divergence-in-population-genomics-and) and perform population structure analyses. 

01_pca_poolfstat.R - Use poolfstat to fiter VCF and perform a PCA on the Pool-seq dataset.

## Part 3 - Genotype environment association analyses using GLM framework


## Part 4 - Identify outlier loci for top abiotic and biotic variable using Baypass


## Part 5 - Genomic offset analyses


## Part 6 - Population genetic simulation using SLiM

- 01_mean_ph:

- 02_Mcali:






