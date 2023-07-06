Versions:

Packages.txt - Lists Python, MATLAB, and R versions and relevant packages with corresponding versions.



Files:

atlasMvec2.csv - Upper triangular matrix vectorization indices.

colenetlabels.txt - Network labels retrieved from: https://github.com/ColeLab/ColeAnticevicNetPartition

arrow3.m - Script for 3D arrow plotting in MATLAB.

cluster_goodness.m - Script for plotting elbow plot in MATLAB.

elbow_k.m - Script for plotting elbow plot in MATLAB.

fit_L_tocurve_area.m - Script for plotting elbow plot in MATLAB.



Analysis:
Step01_r_param_maker.py - Read HCP1200 subjects with resting-state fMRI scans based on the directory structure extracted from DataLad.

Step02_r_fmri_get.txt - For each resting-state subject ID, extract the associated 3T REST1-2 ICA-FIX surface-based fMRI data files from the DataLad repository.

Step03_r_voxdemean.py - For each resting-state scan for each subject ID, subtract mean BOLD timeseries voxelwise.

Step04_r_meants_do.txt - For each resting-state subject ID, use the 360 cortical Glasser atlas to average demeaned voxelwise fMRI data.

Step05_r_motion_get.txt - For each resting-state subject ID, extract the associated 3T REST1-2 fMRI relative FD at each data point from the DataLad repository.

Step06_r_subselector.py - Read in all participants with ICA-FIX 3T resting-state data, then isolate participants who have all four scans. Use the HCP1200 data dictionary to isolate participants who have all 10 cognitive test scores of interest. Retrieve the mean relative FD for each participant for each run and exclude participants with higher than 0.2 in any run. Find the intersection across subjects and output the full list for analysis and a half list to examine parameters.

Step07_c_cogscorer.py - For the subject group, extract 10 cognitive test scores of interest for each subject, then find the first PC from PCA to use as g.

Step08_c_predstat.R - For the subject group, gather g, cognitive tests used to find g, and confounds into one file.

Step09_c_cogboot.R - For the subject group, do bootstrap resampling for the g PC1 loadings.

Step10_r_LE_dFC.m - For the mean ROI fMRI timeseries for each run and each participant, conduct LEiDA.

Step11_r_LE_group_statecalc.m - Read in the subject file for all subjects. For the LE(t) for each run and each participant, stack them and conduct k-medians clustering with a specified k for one random iteration of 500. 

Step12_r_LE_iter_compare.py - For the current subject group for the subject file, find the best iteration of k-medians clustering for each k. Values with lower within-cluster sum of distance (WCS) scores are better.

Step13_r_LE_group_stateanalyze.m - Read in the subject file for all subjects. For the best clustering, generate 2D matrices for LE(t) and dFC(t) centroids by averaging according to the clustering. Next, average while multiplying the dFC(t) by occurrence.

Step14_r_LE_group_statesort.py - For the set of LE or dFC centroids at each k, print out the 2D images and generate a table which sorts them. Then, visually sort them and edit the table. Run this code again with the edited table to visualize the sorted plots.

Step15_r_LE_group_stateclass.py - For the subject file group and selected k, use the sorting matrix to translate the original clustering to the universal indices across k.

Step16_r_LE_group_statestat.m - For the subject group and the specified k, find the dynamics metrics for  each run separately, then average across runs. This includes occurrence, dwell time, transition number, and transition probability.

Step17_r_LE_group_statejump.py - For the subject file group, selected k, and subject ID, find the average transition distance between LE(t) for given clustering transitions within each run, then average across runs.

Step18_r_LE_group_statejump_gather.py - For the subject file group and selected k, gather the transition distances across all subjects.

Step19_r_LE_group_statesim.py - For the subject file group, selected k, and subject ID, find the average idiosyncrasy of each LE(t) from group LE(t) across runs.

Step20_r_LE_group_statesim_gather.py - For the subject file group and selected k, gather the idiosyncrasies across all subjects.

Step21_r_LE_group_staterecon_normtest.R - For the subject group and the specified k, plot the Q-Q and scatter plots for the dynamics metrics and the cognitive variable.

Step22_r_LE_group_state_unicorr.R - For the subject group and the specified k, find the Pearson's correlations between the dynamics metrics and the cognitive variables and the corresponding uncorrected p-values.

Step23_r_LE_group_statejump_unicorr.R - For the subject group and the specified k, find the Pearson's correlations between the transition distance metrics and the cognitive variables and the corresponding uncorrected p-values.

Step24_r_LE_group_statesim_unicorr.R - For the subject group and the specified k, find the Pearson's correlations between the idiosyncrasy metrics and the cognitive variables and the corresponding uncorrected p-values.

Step25_r_LE_group_staterecon_unicorr_FDR.m - For the subject group and the specified k, extract the dynamics, transition distance, and idiosyncrasy p-values. FDR-correct to generate new p-values and threshold correlations with them.

Step26_r_LE_group_staterecon_unicorr_FDR.R - For the subject group and the specified k, show box plots for the distributions of reconfiguration metrics across subjects. Additionally, explore coloring by significant correlations for box plots and bar plots.

Step27_r_LE_group_staterecon_PLSC.R - For the subject group and the specified k, conduct PLSC where X = network reconfiguration metrics and Y = 10 cognitive tests of interest. Do PLSC to find direction of loadings to flip, conduct permutation testing on singular values, conduct bootstrap resampling on loadings, and investigate reproducibility of singular values and loadings.

Step28_r_LE_group_staterecon_PLSC_distplots.R - For the subject group and the specified k, color the distribution of reconfiguration values across subjects in box plots by the PLSC loading with different colors for stable values. 

Step29_r_methods.m - Produce sample dFC(t) and LE(t) to visualize for the methods description.

Step30_r_methods.py - Visualize sample fMRI timeseries, dFC(t), and LE(t) for the methods description.

Step31_r_hilbert_plot.m - Visualize Hilbert transform for sample fMRI timeseries for the methods description.

Step32_r_LE_group_state_vis.py - For the subject file group and selected k, plot the 2D dFC or LE states to display.

Step33_r_LE_group_state_quantradar.R - For the subject group and the specified k, separately average the positive and negative LE(t) values for each region for each network to produce a score.

Step34_r_LE_group_statestd_vis.py - For the subject file group and selected k, plot the dFC variability matrices for each state.

Step35_r_LE_group_state_strvar.m - For the subject group and the specified k, calculate the average dFC(t) strength and variability across edges for each state and across all states for each subject.

Step36_r_LE_group_state_strvar_collect.m - For the subject group and the specified k, gather the average dFC(t) strength and variability across edges for each state and across all states from all subjects.

Step37_r_LE_group_state_strvar.R - For the subject group and the specified k, find the Pearson's correlation between edge-averaged dFC(t) strength and variability across all states, and within each state and plot corresponding scatter plots.

Step38_r_LE_group_state_sFCdFC.m - For the subject group and the specified k, calculate the occurrence-weighted dFC(t) mean and Pearson's correlation static FC for each subject. Calculate the correlation and average across all subjects.

Step39_r_LE_group_state_sFCdFC_vis.py - For the subject file group and selected k, plot the static FC and occurrence-weighted dFC(t) matrices.

Step40_r_LE_group_state_dimtaker_LEfull.py - For the subject file group and selected k, do PCA on the LE(t) across subjects and retrieve the PCs.

Step41_r_LE_group_state_dimtaker_LEfull.R - For the subject group and the specified k, plot the 75% and 25% Q1-3 whiskers for the PCs for each state. For three random subjects, plot LE(t) as dots where colors are states along three PC dimensions.

Step42_r_LE_group_staterecon_highlow_score.R - For the subject group and the specified k, average for the specific set of selected network reconfiguration metrics and rank subjects according to it.

Step43_r_LE_group_staterecon_highlow_normtest.R - For the subject group and the specified k, plot the Q-Q and scatter plots for the sets of summary network reconfiguration metrics and state dFC(t) strength and variability.

Step44_r_LE_group_staterecon_highlow_corr.R - For the subject group and the specified k, get the average for the specific target set of network reconfiguration metrics and find the Pearson's correlations and associated p-values with edge-averaged dFC(t) strength and variability for each state. 

Step45_r_LE_group_staterecon_highlow_corr_FDR.m - For the subject group and the specified k, FDR-correct the p-values for the correlation between strength and variability with reconfiguration set average scores across all sets of interest.

Step46_r_LE_group_state_strvar_highlow_score.R - For the subject group and the specified k, read in the strength and variability data each state and rank subjects according to the average values.

Step47_r_LE_group_state_strvar_highlow_subview.py - For the subject file group and selected k, plot the dFC states for the top and bottom 10% of average strength and variability scores.

Step48_r_LE_group_statek.m - For the subject group, use the group clustering and the centroids for each k to produce the cluster validity indices to produce the elbow plot, where lower values are better.

Step49_r_LE_group_statek_vis.R - For the subject group, plot the cluster validity index for each k to produce the elbow plot.

Step50_r_LE_group_state_brainvis.R - For the subject group, plot the ROIs with positive LE values for each state across the ks of interest on the brain.

Step51_r_colenetwork_vis.R - For each network, delineate it in brain ROIs by color.

