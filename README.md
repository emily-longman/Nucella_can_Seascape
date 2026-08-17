# Revealing the abiotic and biotic drivers of past and future local adaptation in a coastal dogwhelk

## Project Summary

We integrate population genomic data with multiple ecological datasets to determine the ecological stressors most strongly associated with adaptive genomic variation in the channeled dogwhelk, *Nucella canaliculata*, as well as uncover their genetic bases. Our ecological datasets include climatic data on oceanic conditions ([Bio-Oracle](https://www.bio-oracle.org/)), biotic data on the abundance of *N. canaliculata* prey, competitors, and predators from a large-scale monitoring program on Pacific Coast rocky shores (i.e., [MARINe](https://marine.ucsc.edu/): Multi-Agency Rocky Intertidal Network), and a new dataset on the shell traits of the rocky shore mussel *Mytilus californianus*. Leveraging the signals of local adaptation identified, we used genomic offset analyses and population genetic simulations to predict the vulnerability of populations to ocean acidification. 

The bioinformatics pipeline was completed on the Vermont Advanced Computing Center ([VACC](https://www.uvm.edu/vacc)).

### Research Questions

1) What are the abiotic and biotic drivers of geographic patterns of adaptive genomic variation in a low dispersing marine species?
2) What are the genomic bases underlying patterns of local adaptation to varying selection pressures? 
3) Will future environmental change erode existing patterns of local adaptation, or can populations maintain adaptation through evolutionary responses? 

## File Structure

The files for this project are organized in the file structure presented below. All metadata and array guide files are in the "guide_files" directory. The raw data associated with this project are: (1) a previously published VCF of pool-seq data from 19 *N. canaliculata* populations, which can be found on [Zenodo](https://zenodo.org/records/18623551), (2) the *N. canaliculata* reference genome, which can be obtained from [NCBI](https://www.ncbi.nlm.nih.gov/), (3) data on the abundance/density of species that interact with *N. canaliculata* from [MARINe](https://marine.ucsc.edu/), and (4) *M. californianus* shell morphology data, which was generated as part of this project and can be found on Zenodo. All subsequent directories are generated therein. 

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

## Part 1 - Compile ecological data (01_ecological_data)

To identify which ecological selective forces are driving patterns of adaptive variation in *N. canaliculata*, we used several ecologically-relevant abiotic and biotic variables.  

### 01 - Abiotic

We extracted decadal environmental data for the 19 *N. canaliculata* field sites from [Bio-Oracle](https://www.bio-oracle.org/). Our variables of interest were: sea surface temperature (mean, max, min, range, 2010-2019), pH (mean, min, 2010-2018), O2 (mean, 2010-2018), salinity (mean, 2010-2019), and chlorophyll (mean, 2010-2018). For each variable we also extracted future data (2020-2100 under the SSP5-8.5 scenario). Note: the multiple iterations of the three focal scripts were for graphical purposes. 

01_extract_bio-oracle.R - Extract present and future environmental data for a large geographic region along the eastern Pacific coastline.

02_format_bio-oracle.R - Filter each variable for the 19 focal *N. canaliculata* field sites.

03_graph_bio-oracle_raster.R - Graph abiotic variables.

### 02 - Biotic

There were two sets of biotic data. The first was abundance/density data of interacting species (prey, competitors, and predators) from [MARINe](https://marine.ucsc.edu/), 16 of which field sites match the *N. canaliculata* populations. The second is data on the shell morphology of the mussel prey species *Mytilus californianus*, which is a foundation species on rocky shores along the west coast of North America. This latter dataset includes 18 field sites. Mussels for these collections were performed in the same wave-exposed location that the *N. canaliculata* were collected.

01_summarize_graph_MARINe.R - Filter the abundance/density data to the focal *N. canaliculata* field sites and summarize the data.

02_summarize_Mcalifornianus_shell.R - Filter the shell morphology data and summarize.

## Part 2 - Demography and population structure (02_demography)

Filter the Pool-Seq dataset from [Longman et al. (2026)](https://royalsocietypublishing.org/rspb/article/293/2070/20253148/481573/Geographic-divergence-in-population-genomics-and) and perform population structure analyses. 

01_pca_poolfstat.R - Use [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html) to fiter the VCF, creating a poolobject with ~8M SNPs. Subsequently, perform a PCA on the Pool-seq dataset.

## Part 3 - Genotype environment association analyses using GLM framework (03_GEA_glms)

We used generalized linear models (GLMs) to identify associations between *N. canaliculata* allele frequencies and abiotic and biotic ecological variables. We used null expectations, calculated from 100 permutations, to assess the relative importance of the ecological variables.

01_run_glms_abiotic:

- 01_chunks_pt1.txt - Interactive session to make a list of the scaffold names present in the genome.
- 01_chunks_pt2.R - Create a guide file, which chunks the list of scaffold names.
- 02_glms_bio-oracle_launch.sh and 02_glms_bio-oracle.R - Use GLMs to assess associations between allele frequencies (chunked based on the scaffold lists created in the previous scripts) and abiotic [Bio-Oracle](https://www.bio-oracle.org/) variables 
- 03_extract_glms.R and 03_extract_glms_launch.R - Extract the data for each abiotic variable from the GLM outputs.
- 04_merge_glms.R and 04_merge_glms_launch.R - Merge all of the GLM outputs for each abiotic variable.

01_run_glms_biotic:

- 01_MARINe

    - 01_glms_MARINe_launch.sh and 01_glms_MARINe_*.R - Use GLMs to assess associations between allele frequencies (chunked based on the scaffold guide file) and biotic [MARINe](https://marine.ucsc.edu/) variables. For efficiency, there are 5 R scripts, one that performs GLMs with the real ecological variables and four scripts that performs GLMs with the permuted variables. Note the different walltime for the real data versus permutations.
    - 02_extract_glms.R and 02_extract_glms_launch.R - Extract the data for each MARINe variable from the GLM outputs.
    - 03_merge_glms.R and 03_merge_glms_launch.R - Merge all of the GLM outputs for each MARINe variable.

- 02_Mcali

    - 01_glms_Mcali_launch.sh and 01_glms_Mcali_*.R - Use GLMs to assess associations between allele frequencies (chunked based on the scaffold guide file) and the mussel morphology variables. For efficiency, there are 5 R scripts, one that performs GLMs with the real data and four scripts that performs GLMs with the permuted variables. Note the different walltime for the real data versus permutations.
    - 02_extract_glms.R and 02_extract_glms_launch.R - Extract the data for each mussel morphology variable from the GLM outputs.
    - 03_merge_glms.R and 03_merge_glms_launch.R - Merge all of the GLM outputs for each mussel morphology variable.


02_model_enrichment:

- 01_filt_glms_launch.sh and 01_filt_glm.R - For all of the ecological variables, filter the data to the ~8M SNP list present in the poolobject.
- 02_model_enrichment_launch.sh and 02_model_enrichment.R - For each variable, summarize the permutations and perform model enrichment.
- 03_graph_glms.R - Graph the GLM results. 

## Part 4 - Identify outlier loci for top abiotic and biotic variable using Baypass (04_baypass)

To identify outlier loci associated with mean pH and mussel shell thickness we performed genotype-environment association scans while simultaneously accounting for demographic history using [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public). 

01_format_baypass.R - Use [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html) to generate the Baypass input files for the entire 19 population dataset, as well as the 18 population subset used in the mussel shell thickness analyses.

02_generate_omega.sh and 02_generate_omega_subset18pop.sh - Use [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) to generate the omega relatedness matrix for the appropriate number of populations.

03_baypass_Mcali_Thk.sh - Perform 5 runs of [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) using the mussel cross-sectional thickness data as a covariate. 

03_baypass_ph_mean.sh - Perform 5 runs of [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) using mean pH data as a covariate. 

04_PODs.R and 04_PODs_subset_18pop.R - Use [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html) to create 10 pseudo-observed dataset (PODs).

05_GEA_PODs:  

- 01_baypass_PODs_ph_mean.sh - Perform 10 runs of BayPass using the PODs as the input and mean pH as a covariate.
- 01_baypass_Mcali_Thk.sh - Perform 10 runs of BayPass using the PODs as the input and mussel shell thickness as a covariate.
- 02_calibrate_baypass.R - Summarize the POD outputs and generate 95%, 99% and 99.9% significance thresholds for both ecological variables.

06_summarize_graph_baypass.R - Analyze the BayPass data for both ecological variables.

07_annotate_SNPs_Mcali.R - Annotate the outlier SNPs associated with mussel cross-sectional shell thickness using [SnpEff](https://pcingola.github.io/SnpEff/).

07_annotate_SNPs_ph_mean.R - Annotate the outlier SNPs associated with mean pH using [SnpEff](https://pcingola.github.io/SnpEff/).

08_window_rnp_100kb: 

- 01_generate_window.R - Create 100kb windows for the genome scan.
- 02_baypass_window_rnp_Mcali_launch.sh and 02_baypass_window_rnp_Mcali.R - Perform a window scan to identify regions with the genome that are enriched with outlier SNPs associated with mussel shell thickness.
- 02_baypass_window_rnp_ph_mean_launch.sh and 02_baypass_window_rnp_ph_mean.R - Perform a window scan to identify regions with the genome that are enriched with outlier SNPs associated with mean pH.
- 03_graph_win_Mcali.R - Graph the window scan output for mussel shell thickness.
- 03_graph_win_ph_mean.R - Graph the window scan output for mean pH.
- 04_graph_ph_mean_outliers_win.R - For the top outlier window that is enriched with SNPs associated with mean pH, analyze and graph Bayes Factors and genetic differentiation (XTX). Additionally, graph allele frequencies for the SNPs with the highest association values.

09_post_hoc_outlier_analyses:

- 01_annotate_all_launch.sh and 01_annotate_all.R - Annotate the entire VCF, using a chunked 500-array approach.
- 02_merge_annotate.R - Merge the output from the previous script, and filter for the ~8M SNP dataset.
- 03_fishers_exact_test.R - Perform Fisher’s Exact Tests to assess if the mean pH model and shell thickness model were enriched with specific genetic variants.
- 04_baypass_corrected_AF_Mcali_launch.sh and 04_baypass_corrected_AF_Mcali.R - Summarize the baypass corrected allele frequencies (i.e., the mean of the posterior distribution of the alpha ij parameter) for the 5 shell thickness runs.
- 04_baypass_corrected_AF_ph_launch.sh and 04_baypass_corrected_AF_ph.R - Summarize the baypass corrected allele frequencies (i.e., the mean of the posterior distribution of the alpha ij parameter) for the 5 mean pH runs.
- 05_Mcali_AF_cor.R - Calculate Pearson’s Correlation between the residual structure of the top outlier loci associated with shell thickness and a random sample of 1,000 SNPs.
- 05_ph_AF_cor.R - Calculate Pearson’s Correlation between the residual structure of the top 51 outlier loci associated with mean pH and a random sample of 1,000 SNPs.
- 06_effect_sizes.R - Calculate local effect size (Cohen's F2) for individual SNPs for both ecological models. 


## Part 5 - Genomic offset analyses (05_genomic_offset)

Calculate genome-wide geometric genomic offset using [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) based on scaled future projections of mean pH from [Bio-Oracle](https://www.bio-oracle.org/).

01_scale_mean_ph.R - Scale mean pH for current and future data.

02_Baypass_GO_mean_ph_scaled.R - Perform genomic offset analyses for the 19 populations based on changes in mean pH.

## Part 6 - Population genetic simulation using SLiM (06_SLiM)

We performed population genetic simulations using [SLiM](https://github.com/MesserLab/SLiM) to contextualize the patterns of local adaptation. To do so, we created single locus population genetic simulations using a Wright-Fisher model based on the top outlier loci associated with both mean pH and mussel shell thickness. In the models, we calculated population-level selection using a sigmoidal function, with multiple parameters that dictate the structure of the curves. We used Approximate Bayesian Computation [ABC](https://cran.r-project.org/web/packages/abc/index.html) to perform multivariate parameter estimation between the simulations and real data.  

To determine the adaptive potential of populations to ocean acidification, we assessed how genetic diversity would change with declining pH.

01_mean_ph:
 
- 01_guide_file.R - Create a guide file of the parameters. 
- 02_ph.slim and 02_ph_launch_slim.sh - For the parameters in the guide file, use SLiM to perform population genetic simulations for an array of parameters.
- 03_ABC_SimData_luanch.R and 03_ABC_SimData.R - Calculate the summary statistics for ABC for the simulated data.
- 04_ABC_realData.R - Calculate the summary statistics for ABC for the real data.
- 05_ABC_analysis_LOclinEst.R - Perform ABC analysis. 
- 06_format_ph_future.R - Use Bio-Oracle mean pH data for 2020-2100 to generate population regression equations.
- 07_ph_future_slim_launch.sh and 07_ph_future.slim - Perform population genetic simulation for 2020-2100. Populations were started at the empirical allele frequencies and the pH regression equations were used to calculate yearly population-level selection coefficients.
- 08_analyze_ph_future.R - Compile simulations and graph output. 
- 09_ph_future_vary_m_launch.sh and 09_ph_future_vary_m.slim - Perform simulations with varying migration rates. 
- 10_analyze_ph_future_vary_m.R - Analyze and graph the output of the population genetic simulations generated in the previous script.

02_Mcali:

- 01_guide_file.R - Create a guide file of the parameters. 
- 02_ph_launch_slim.sh and 02_ph.slim - For the parameters in the guide file, use SLiM to perform population genetic simulations for an array of parameters.
- 03_ABC_SimData_launch.sh and 03_ABC_SimData.R - Calculate the summary statistics for ABC for the simulated data.
- 04_ABC_realData.R - Calculate the summary statistics for ABC for the real data.
- 05_ABC_analysis_LOclinEst.R - Perform ABC analysis. 

03_graph_selection_ph_Mcali.R - Graph the selection curves for the abiotic and biotic data based on the best fit parameters.




