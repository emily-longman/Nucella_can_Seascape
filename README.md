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

GLMs...

01_run_glms_abiotic:

- 01_chunks_pt1.txt - Interactive session to make list of scaffold names of genome
- 01_chunks_pt2.R - Create guide file of scaffold chunks for GLMs.
- 02_glms_bio-oracle_pt1.sh
- 02_glms_bio-oracle_pt2.sh
- 02_glms_bio-oracle_pt3.sh
- 02_glms_bio-oracle.R
- 03_extract_glms.R and 03_extract_glms_launch.R
- 04_merge_glms.R and 04_merge_glms_launch.R

01_run_glms_biotic:

- 01_MARINe
- 02_Mcali

02_model_enrichment:


## Part 4 - Identify outlier loci for top abiotic and biotic variable using Baypass

To identify outlier loci associated with mean pH and mussel shell thickness we performed genotype-environment association scans while simultaneously accounting for demographic history using [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public). 

01_format_baypass.R - 

02_generate_omega.sh - 

03_baypass_Mcali_Thk.sh - 

03_baypass_ph_mean.sh - 

04_PODs.R - 

05_GEA_PODs:  

- 01_baypass_PODs_ph_mean.sh -
- 01_baypass_Mcali_Thk.sh -
- 02_calibrate_baypass.R -

06_summarize_graph_baypass.R -

07_annotate_SNPs_Mcali_thk.R - 

07_annotate_SNPs_ph_mean.R - 

08_window_rnp_100kb: 

09_post_hoc_outlier_analyses:



## Part 5 - Genomic offset analyses

Calculate genome-wide geometric genomic offset using [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) based on scaled future projections of mean pH from [Bio-Oracle](https://www.bio-oracle.org/).

01_scale_mean_ph.R - Scale mean pH for current and future data.

02_Baypass_GO_mean_ph_scaled.R - Perform genomic offset analyses for the 19 populations based on changes in mean pH.

## Part 6 - Population genetic simulation using SLiM

Perform population genetic simulations using [SLiM](https://github.com/MesserLab/SLiM). First, we created single locus population genetic simulations using a Wright-Fisher model for both mean pH and mussel shell thickness. In the models, we calculated population-level selection using a sigmoidal function, with parameters that dictate the structure of the curves. We used Approximate Bayesian Computation [abc](https://cran.r-project.org/web/packages/abc/index.html) to performs multivariate parameter estimation between the simulations and real data.  Second, we used the structure of the mean pH model to determine how genetic diversity would change with declining pH, as a result of ocean acidification over the time period 2020-2100.

01_mean_ph:
 
- 01_guide_file.R - Create a guide file of the parameters. 
- 02_ph.slim and 02_ph_launch_slim.sh - For the parameters in the guide file, use SLiM to perform simulations for an array of parameters.
- 03_ABC_SimData.R and 03_ABC_SimData_launch.sh - Calculate the summary statistics for ABC for the simulations.
- 04_ABC_realData.R - Calculate the summary statistics for ABC for the real data.
- 05_ABC_analysis_LOclinEst.R - Perform ABC analysis. 
- 06_format_ph_future.R - Use Bio-Oracle mean pH data for 2020-2100 to calculate population regression equations.
- 07_ph_future.slim and 07_ph_future_slim_launch.sh - Perform population genetic simulation for 2020-2100. Populations were started at the empirical allele frequencies and the pH regression equations were used to calculate yearly population-level selection coefficients.
- 08_analyze_ph_future.R - Compile simulations and graph output. 
- 09_ph_future_vary_m.slim and 09_ph_future_vary_m_launch.slim - Perform simulations with varying migration rates. 
- 10_analyze_ph_future_vary_m.R - Analyze and graph the output of the previous population genetic simulations.

02_Mcali:

- 01_guide_file.R - Create a guide file of the parameters. 
- 02_ph.slim and 02_ph_launch_slim.sh - For the parameters in the guide file, use SLiM to perform simulations for an array of parameters.
- 03_ABC_SimData.R and 03_ABC_SimData_launch.sh - Calculate the summary statistics for ABC for the simulations.
- 04_ABC_realData.R - Calculate the summary statistics for ABC for the real data.
- 05_ABC_analysis_LOclinEst.R - Perform ABC analysis. 

03_graph_selection_ph_Mcali.R - Graph the selection curves for the abiotic and biotic data based on the best fit parameters.




