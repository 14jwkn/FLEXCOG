%{
For the subject group and the specified k, calculate the
occurrence-weighted dFC(t) mean and Pearson's correlation static FC for
each subject. Calculate the correlation and average across all subjects.
Adapts code from: 
https://github.com/juanitacabral/LEiDA
https://github.com/trendscenter/gift
Output:
sub_mean_dFC.h5 Occurrence-weighted dFC(t) mean across all subjects, vectorized and matrix.
sub_mean_sFC.h5 Average ROI correlation static FC mean across all subjects, vectorized and matrix.
subcorr_sFCdFC.h5 Correlation between occurence-weighted dFC(t) mean and static FC matrix for each subject.
subcorr_sFCdFC.csv Correlation between occurence-weighted dFC(t) mean and static FC matrix averaged across all subjects.
%}

%Define command line arguments.
function [] = 37_r_LE_group_state_sFCdFC(subfile,k)
disp(append('Doing ',subfile,' ',k));

%Set up I/O. 
if strcmp(subfile,'r_full_submain.txt') 
    subgroup = 'full';  
elseif strcmp(subfile,'r_half_submain.txt') 
    subgroup = 'half'; 
end        
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                subgroup,'/',k,'/');
subdFC_file = append(outpath,'sub_mean_dFC.h5');
subsFC_file = append(outpath,'sub_mean_sFC.h5');
sFCdFC_file = append(outpath,'subcorr_sFCdFC.h5');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
nsubj = size(subjects,1);

%Read in the best iteration.
best_file = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/',...
                  'group/best_iter/',subgroup,'/best_iter.csv');
best_clust = readtable(best_file);
best_lab = k;
iteration = num2str(best_clust{:,best_lab});
disp(append('Best: ',iteration))

%Read in best clustering and centroids, then initialize values.
clustfile = append(outpath,'subclust_',iteration,'.h5');
centfile = append(outpath,'subcent_',iteration,'.h5');
inkey = append('/',subgroup);
IDX = h5read(clustfile,inkey);
V = h5read(centfile,inkey);
Tmax = size(IDX,1);
[N_Cl,N_ba] = size(V);

%Calculate occurrence.
F = zeros(1,N_Cl);
for jj = 1:N_Cl
    F(jj) = (sum(IDX == jj))/Tmax;
end

%Do dFC weighted mean.
nconn = (N_ba*(N_ba-1))/2;
all_dFC_mean = zeros(nsubj,nconn);
for subin = 1:nsubj

    %Initialize subject dFC.
    dFC_mean = zeros(1,nconn);
        
    %Read in subject data.
    disp(append('Doing dFC: ',subjects(subin)))
    inkey = '/LE_dFC';
    Lfile1 = append('../outputs/r_LE_dFC/REST1_LR/',subjects(subin),...
                    '/LE_dFC.h5');     
    Rfile1 = append('../outputs/r_LE_dFC/REST1_RL/',subjects(subin),...
                    '/LE_dFC.h5'); 
    Lfile2 = append('../outputs/r_LE_dFC/REST2_LR/',subjects(subin),...
                    '/LE_dFC.h5'); 
    Rfile2 = append('../outputs/r_LE_dFC/REST2_RL/',subjects(subin),...
                    '/LE_dFC.h5'); 
    dFCmatL1 = h5read(Lfile1,inkey);
    dFCmatR1 = h5read(Rfile1,inkey);
    dFCmatL2 = h5read(Lfile2,inkey);
    dFCmatR2 = h5read(Rfile2,inkey);
    dFCmat = vertcat(dFCmatL1,dFCmatR1,dFCmatL2,dFCmatR2);
    delvars = {'dFCmatL1','dFCmatR1','dFCmatL2','dFCmatR2'};
    clear(delvars{:})  
    nwin = size(dFCmat,1);

    %Set time indices and go through them.
    startin = (nwin*(subin-1)) + 1;
    endin = nwin*subin;
    subt = 1;
    disp(append('Subject t: ',num2str(startin),' to ',num2str(endin)));
    for t=startin:endin

        %Add to the running total weighted mean.
        idx_c = IDX(t);
        dFC_mean = dFC_mean + (dFCmat(subt,:)*F(idx_c));

        %Increment subject time.
        subt = subt + 1; 
    end
    delvars = {'dFCmat'};
    clear(delvars{:})
    
    %Append.
    all_dFC_mean(subin,:) = dFC_mean;
end

%Save.
outkey = append('/sub_1D');
try
    h5create(subdFC_file,outkey,size(all_dFC_mean));
catch
    disp('File and key created already.')
end
h5write(subdFC_file,outkey,all_dFC_mean);
disp('Subject dFC saved.')

%Find the dFC average, unvectorize and save it.
dFC_mean_2D = zeros(N_ba,N_ba);
dFC_mean_2D(triu(true(N_ba),1)) = mean(all_dFC_mean);
dFC_mean_2D = dFC_mean_2D + dFC_mean_2D';
outkey = append('/mean_2D');
try
    h5create(subdFC_file,outkey,size(dFC_mean_2D));
catch
    disp('File and key created already.')
end
h5write(subdFC_file,outkey,dFC_mean_2D);
disp('Group dFC saved.')

%Set 2D reverser.
I_sup_diag=find(triu(ones(N_ba),1));

%Do sFC.
all_sFC = zeros(nsubj,nconn);
for subin = 1:nsubj
    
    %Set subject indices.
    disp(append('Doing sFC: ',subjects(subin)))
    
    %Read in the subject file.
    inpath = append('../outputs/r_meants/',subjects(subin),'/');
    prefix = append('demean_rfMRI_');
    
    %Set file names.
    Lfile1 = append(inpath,prefix,'REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');
    Rfile1 = append(inpath,prefix,'REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');              
    Lfile2 = append(inpath,prefix,'REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');
    Rfile2 = append(inpath,prefix,'REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');  
    
    %Read and remove first and last time points.
    sub_tc_L1 = readmatrix(Lfile1)';
    sub_tc_L1(1200,:) = [];
	sub_tc_L1(1,:) = [];
    sub_tc_R1 = readmatrix(Rfile1)';
    sub_tc_R1(1200,:) = [];
	sub_tc_R1(1,:) = [];
    sub_tc_L2 = readmatrix(Lfile2)';
    sub_tc_L2(1200,:) = [];
	sub_tc_L2(1,:) = [];
    sub_tc_R2 = readmatrix(Rfile2)';
    sub_tc_R2(1200,:) = [];
	sub_tc_R2(1,:) = [];
    
    %Correlate.
    corr_L1 = corrcoef(sub_tc_L1);
    corr_L1 = corr_L1(I_sup_diag);
    corr_R1 = corrcoef(sub_tc_R1);
    corr_R1 = corr_R1(I_sup_diag);
    corr_L2 = corrcoef(sub_tc_L2);
    corr_L2 = corr_L2(I_sup_diag);
    corr_R2 = corrcoef(sub_tc_R2);
    corr_R2 = corr_R2(I_sup_diag);
    
    %Concatenate, average, and add.
    subcorr = horzcat(corr_L1,corr_R1,corr_L2,corr_R2);
    subcorr = mean(subcorr');
    all_sFC(subin,:) = subcorr;
    delvars = {'sub_tc','sub_tc_L1','sub_tc_R1','sub_tc_L2','sub_tc_R2'};
    clear(delvars{:})   
end 

%Save.
outkey = append('/sub_1D');
try
    h5create(subsFC_file,outkey,size(all_sFC));
catch
    disp('File and key created already.')
end
h5write(subsFC_file,outkey,all_sFC);
disp('Subject sFC saved.')

%Find the dFC average, unvectorize and save it.
sFC_2D = zeros(N_ba,N_ba);
sFC_2D(triu(true(N_ba),1)) = mean(all_sFC);
sFC_2D = sFC_2D + sFC_2D';
outkey = append('/mean_2D');
try
    h5create(subsFC_file,outkey,size(sFC_2D));
catch
    disp('File and key created already.')
end
h5write(subsFC_file,outkey,sFC_2D);
disp('Group sFC saved.')

%Find the correlation between mean dFC and sFC for each subject.
allcorr = zeros(nsubj,1);
for subin = 1:nsubj
    
    %Correlate.
    disp(append('Doing sFC-dFC correlation: ',subjects(subin)))
    [cc,~]=corrcoef(all_sFC(subin,:),all_dFC_mean(subin,:));
    allcorr(subin,1) = cc(1,2);
end

%Save.
outkey = append('/sub');
try
    h5create(sFCdFC_file,outkey,size(allcorr));
catch
    disp('File and key created already.')
end
h5write(sFCdFC_file,outkey,allcorr);
disp('Subject sFC-dFC correlation saved.')

%Find the average and save it.
corrmean = array2table(mean(allcorr));
outfile = append(outpath,'subcorr_sFCdFC.csv');
writetable(corrmean,outfile)
disp('Group sFC-dFC correlation saved.')
end
