Versions:

Packages.txt - Lists Python, MATLAB, and R versions and relevant packages with corresponding versions.

----------
Files:

atlasMvec2.csv - Upper triangular matrix vectorization indices.

colenetlabels.txt - Network labels retrieved from: https://github.com/ColeLab/ColeAnticevicNetPartition

arrow3.m - Script for 3D arrow plotting in MATLAB retrieved from: https://github.com/juanitacabral/LEiDA_Psilocybin

cluster_goodness.m - Script for plotting elbow plot in MATLAB retrieved from: https://github.com/trendscenter/gift

fit_L_tocurve_area.m - Script for plotting elbow plot in MATLAB retrieved from: https://github.com/trendscenter/gift

elbow_k.m - Script for plotting elbow plot in MATLAB which uses cluster_goodness.m and fit_L_tocurve_area.m.

----------
Analysis:

1. Extract ICA-FIX resting-state fMRI data and motion.

    param_maker.py - Read HCP1200 subjects with resting-state fMRI scans based on the directory structure extracted from DataLad.

    fmri_get.txt - For each resting-state subject ID, extract the associated 3T REST1-2 ICA-FIX surface-based fMRI data files from the DataLad repository.

    voxdemean.py - For each resting-state scan for each subject ID, subtract mean BOLD timeseries voxelwise.

    meants_do.txt - For each resting-state subject ID, use the 360 cortical Glasser atlas to average demeaned voxelwise fMRI data.

    motion_get.txt - For each resting-state subject ID, extract the associated 3T REST1-2 fMRI relative FD at each data point from the DataLad repository.

2. Select subjects.

    subselector.py - Read in all participants with ICA-FIX 3T resting-state data, then isolate participants who have all four scans. Use the HCP1200 data dictionary to isolate participants who have all 10 cognitive test scores of interest. Retrieve the mean relative FD for each participant for each run and exclude participants with higher than 0.2 in any run. Find the intersection across subjects.

3. Do LEiDA, extract reconfiguration metrics, and assess ICC.

    LE_dFC.m - For the mean ROI fMRI timeseries for each run and each participant, conduct LEiDA.

    LE_group_statecalc.m - For the LE(t) for each run and each participant, stack them and conduct k-medians clustering with a specified k for one random iteration of 500. 

    LE_group_iter_compare.py - Find the best iteration of k-medians clustering for each k. Values with lower within-cluster sum of distance (WCS) scores are better.

    LE_group_stateanalyze.m - For the best clustering, generate 2D matrices for LE(t) centroids and dFC(t) centroids by averaging according to the clustering. 

    LE_group_statesort.py - For the set of LE or dFC centroids across a range of k, print out the 2D images and generate a table which sorts them. Then, visually sort them and edit the table. Run this code again with the edited table to visualize the sorted plots.

    LE_group_stateclass.py - For the selected k, use the sorting matrix to translate the original clustering to the universal indices across k.

    LE_group_statestat.m - For the specified k, find the frequency metrics for each run separately, then average across runs. This includes occurrence, dwell time, transition number, and transition probability.

    LE_group_statejump.py - For the selected k and subject ID, find the average transition distance between LE(t) for given clustering transitions within each run, and average across runs.

    LE_group_statejump_gather.py - For the selected k, gather the transition distances across all subjects across runs and within run.

    LE_group_statesim.py - For the selected k and subject ID, find the average idiosyncrasy of each LE(t) from group LE(t) across runs and within run.

    LE_group_statesim_gather.py - For the selected k, gather the idiosyncrasies across all subjects across runs and within run.

    LE_group_staterecon_ICC.R - For the specified k, find ICC(1,1) between the session averages across runs for each network reconfiguration metric.
    
4. Extract cognitive scores and confounds.

    predstat.R - Gather cognitive tests and confounds into one file.

5. Do PLSC.

    r_LE_group_staterecon_PLSC.R - For the specified k, conduct PLSC where X = network reconfiguration metrics and Y = 10 cognitive tests of interest. Do PLSC to find direction of loadings to flip, conduct permutation testing on singular values, investigate bivariate normality and outliers for significant latent variable pairs and correlate, conduct bootstrap resampling on loadings, and investigate reproducibility of singular values and loadings.

    r_LE_group_staterecon_PLSC_distplots.R - For the specified k, color the distribution of reconfiguration values across subjects in box plots by the PLSC loading with different colors for stable values. 

6. Do univariate analyses.

    cogscorer_PCA.py - Do PCA to extract the first PC score and loadings to use as gPCA. 

    subselector_FA.py - Relax the fMRI requirement to increase sample size to improve cognitive score extraction.  Get cognitive score data for expanded sample and add gPCA to it.

    cogscorer_FA.R - For the expanded sample, conduct EFA and CFA and retrieve scores and loadings.

    predstat_expanded.R - Gather cognitive tests and confounds into one file, with the addition of gPCA, gEFA, and gCFA.

    cogcorr.R - For the specified k, find the correlations between gPLSC, gPCA, gEFA, gCFA, and the 10 input cognitive variables.

    LE_group_staterecon_normtest.R - For the specified k, for pairs of cognitive and network reconfiguration variables, conduct bivariate normality tests, plot scatterplots, and plot QQ plots.

    LE_group_state_unicorr.R - For the specified k and cognitive variable, find the correlations with the frequency metrics and the corresponding uncorrected p-values.

    LE_group_statejump_unicorr.R - For the specified k and cognitive variable, find the correlations with the transition distance metrics and the corresponding uncorrected p-values.

    LE_group_statesim_unicorr.R - For the specified k and cognitive variable, find the correlations with the idiosyncrasy metrics and the corresponding uncorrected p-values.

    LE_group_staterecon_unicorr_FDR.m - For the specified k and cognitive variable, extract the frequency, transition distance, and idiosyncrasy p-values. FDR-correct to generate new p-values and threshold correlations with them.

    LE_group_staterecon_unicorr_FDR.R - # For the specified k and cognitive variable, show box plots for the distributions of reconfiguration metrics across subjects. Additionally, explore coloring by significant correlations for box plots and bar plots.

7. Do summary score co-existence analysis.

    LE_group_staterecon_summary_score.R - For the specified k, average for the specific set of selected network reconfiguration metrics.

    LE_group_staterecon_summary_normtest.R - For the specified k, for pairs of summary scores relating to g or processing speed, conduct bivariate normality tests, plot scatterplots, and plot QQ plots.

    LE_group_staterecon_summary_corr.R - For the specified k, for pairs of summary scores relating to g or processing speed, find the correlations and associated p-values.

    LE_group_staterecon_summary_corr_FDR.m - For the specified k, FDR-correct the p-values for the correlation between pairs of summary scores relating to g or processing speed.

8. Visualize, for a specified k, the dFC and LE states, and average positive and negative LE scores in radar plots.

    r_LE_group_state_vis.py - For the selected k, plot the 2D dFC or LE states to display with and without limits.

    r_LE_group_state_quantradar.R - For the specified k, separately average the positive and negative LE values for each region for each network to produce a score. 

9. Visualize, for a range of k, the elbow plot, and positive LE scores on the brain.

    LE_group_state_elbow.m - Use the group clustering and the centroids for each k to produce the cluster validity indices to produce the elbow plot, where lower values are better.

    LE_group_state_elbow_vis.R - Plot the cluster validity index for each k to produce the elbow plot.

    LE_group_state_brainvis.R - Plot the ROIs with positive LE values for each state across the ks of interest on the brain.

10. Visualize ICNs on the brain.

    colenetwork_vis.R - For each network, delineate it in brain ROIs by color.

11. Visualize FC strength and variability.

    LE_group_state_FC_strvar.m - For the specified k and participant, extract dFC(t) across and within states and within runs, calculate mean for FC strength and std for FC variability, average across runs.

    LE_group_state_FC_strvar_collect.m - For the specified k, collect the dFC FC strength and variability across all states and within states from all participants. Generate averaged 2D matrices from vectorized forms.

    LE_group_state_FC_strvar_vis.py - For the selected k, calculate the correlation between total FC strength and variability across states, plot the FC strength and variability within each state as 2D matrices with limits, and plot the FC strength and variability within each state as histograms.

12. Visualize methods figure items.

    methods.m - Produce sample dFC(t) and LE(t) to visualize for the methods description.

    methods.py - Visualize sample fMRI timeseries, dFC(t), and LE(t) for the methods description.

    hilbert_plot.m - Visualize Hilbert transform for sample fMRI timeseries for the methods description.
