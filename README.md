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

There were two sets of biotic data. The first was abundance/density data of interacting species (prey, competitors, and predators) from [MARINe](https://marine.ucsc.edu/), 16 of which match the *N. canaliculata* populations. The second is data on the shell morphology of the mussel prey species *Mytilus californianus* which is a foundation species on rocky shores along the west coast of North America. This latter dataset includes 18 field sites, which match the *N. canaliculata* populations.

01_summarize_graph_MARINe.R - Filter abundance/density data of interacting species to focal sites and summarize.

02_summarize_Mcalifornianus_shell.R - Filter the shell morphology data and summarize.

## Part 2 - Demography and population structure

Filter the Pool-Seq dataset from [Longman et al. (2026)](https://royalsocietypublishing.org/rspb/article/293/2070/20253148/481573/Geographic-divergence-in-population-genomics-and) and perform population structure analyses. 

01_pca_poolfstat.R - Use [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html) to fiter VCF, creating a poolobject with ~8M SNPs, and perform a PCA on the Pool-seq dataset.

## Part 3 - Genotype environment association analyses using GLM framework

We used generalized linear models (GLMs) to identify associations between *N. canaliculata* allele frequencies and abiotic and biotic ecological variables. We used null expectations, calculated from 100 permutations, to assess the relative importance of the ecological variables.

01_run_glms_abiotic:

- 01_chunks_pt1.txt - Interactive session to make list of scaffold names of genome
- 01_chunks_pt2.R - Create guide file of scaffold chunks for GLMs.
- 02_glms_bio-oracle_launch.sh and 02_glms_bio-oracle.R - Use GLMs to assess associations between allele frequencies and abiotic [Bio-Oracle](https://www.bio-oracle.org/) variables.
- 03_extract_glms.R and 03_extract_glms_launch.R - For the GLM output extract the data for each abiotic variable.
- 04_merge_glms.R and 04_merge_glms_launch.R - Merge all of the GLM outputs for each abiotic variable.

01_run_glms_biotic:

- 01_MARINe

    - 01_glms_MARINe_launch.sh and 01_glms_MARINe_*.R - Use GLMs to assess associations between allele frequencies and biotic [MARINe](https://marine.ucsc.edu/) variables. For efficiency, there are 5 R scripts, one that permforms GLMs with the real data and four scripts that peforms GLMs with the permuted variables. Note, the different walltime for the real data versus permutations.
    - 02_extract_glms.R and 02_extract_glms_launch.R - For the GLM output extract the data for each MARINe variable.
    - 03_merge_glms.R and 03_merge_glms_launch.R - Merge all of the GLM outputs for each MARINe variable.

- 02_Mcali

    - 01_glms_Mcali_launch.sh and 01_glms_Mcali_*.R - Use GLMs to assess associations between allele frequencies and biotic mussel morphology variables. For efficiency, there are 5 R scripts, one that permforms GLMs with the real data and four scripts that peforms GLMs with the permuted variables. Note, the different walltime for the real data versus permutations.
    - 02_extract_glms.R and 02_extract_glms_launch.R - For the GLM output extract the data for each mussel morphology variable.
    - 03_merge_glms.R and 03_merge_glms_launch.R - Merge all of the GLM outputs for each mussel morphology variable.


02_model_enrichment:

- 01_filt_glms_launch.sh and 01_filt_glm.R - For all of the ecological variables, filter the data so only the SNPs in the poolobject rather than the entire VCF.
- 02_model_enrichment_launch.sh and 02_model_enrichment.R - For each variable, summarize the permutations and perform model enrichment.
- 03_graph_glms.R - Graph the GLM results. 

## Part 4 - Identify outlier loci for top abiotic and biotic variable using Baypass

To identify outlier loci associated with mean pH and mussel shell thickness we performed genotype-environment association scans while simultaneously accounting for demographic history using [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public). 

01_format_baypass.R - Use [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html) to generate the Baypass input files for the entire 19 population dataset, and the 18 populations that we included in the mussel shell thickness analyses.

02_generate_omega.sh and  02_generate_omega_subset18pop.sh - Use [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) to generate the omega relatedness matrix for the appropriate number of populations.

03_baypass_Mcali_Thk.sh - Perform 5 runs of [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) using the mussel cross-sectional thickness data as a covariate. 

03_baypass_ph_mean.sh - Perform 5 runs of [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) using mean pH data as a covariate. 

04_PODs.R and 04_PODs_subset_18pop.R - Use [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html) create 10 pseudo-observed dataset (PODs).

05_GEA_PODs:  

- 01_baypass_PODs_ph_mean.sh - Perform 10 runs of BayPass using the POD as the input and mean pH as a covariate.
- 01_baypass_Mcali_Thk.sh - Perform 10 runs of BayPass using the POD as the input and mussel shell thickness as a covariate.
- 02_calibrate_baypass.R - Summarize the POD output and generate 95%, 99% and 99.9% thresholds for both ecological variables.

06_summarize_graph_baypass.R - Analyze the BayPass data for both ecological variables.

07_annotate_SNPs_Mcali.R - Annotate the outlier SNPs associated with mussel cross-sectional shell thickness using [SnpEff](https://pcingola.github.io/SnpEff/).

07_annotate_SNPs_ph_mean.R - Annotate the outlier SNPs associated with mean pH using [SnpEff](https://pcingola.github.io/SnpEff/).

08_window_rnp_100kb: 

- 01_generate_window.R - Create 100kb windows for genome scan.
- 02_baypass_window_rnp_Mcali.sh and 02_baypass_window_rnp_Mcali.R - Peform window scan to identify regions with the genome that are enriched with outlier SNPs associated with mussel shell thickness.
- 02_baypass_window_rnp_ph_mean.sh and 02_baypass_window_rnp_ph_mean.R - Peform window scan to identify regions with the genome that are enriched with outlier SNPs associated with mean pH.
- 03_graph_win_Mcali.R - Graph window scan output for mussel shell thickness.
- 03_graph_win_ph_mean.R - Graph window scan output for mean pH.
- 04_graph_ph_mean_outliers_win.R - Analyze and graph the Bayes Factors and XTX for the SNPs within the top outlier window. Additionally, graph allele frequencies for the SNPs with the highest association values.

09_post_hoc_outlier_analyses:

- 01_annotate_all_launch.sh and 01_annotate_all.R - Annotate the entire VCF, using a chunked 500-array approach.
- 02_merge_annotate.R - Merge the output from the previous scripts, and filter for the ~8M SNP dataset.
- 03_fishers_exact_test.R - Perform Fisher’s Exact Tests to assess if the mean pH model and shell thickness model were enriched with specific genetic variants.
- 04_baypass_corrected_AF_Mcali_launch.sh and 04_baypass_corrected_AF_Mcali.R - Summarize the baypass corrected allele frequencies (i.e., the mean of the posterior distribution of the alpha ij parameter) for the 5 shell thickness runs.
- 04_baypass_corrected_AF_ph_launch.sh and 04_baypass_corrected_AF_ph.R - Summarize the baypass corrected allele frequencies (i.e., the mean of the posterior distribution of the alpha ij parameter) for the 5 mean pH runs.

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




