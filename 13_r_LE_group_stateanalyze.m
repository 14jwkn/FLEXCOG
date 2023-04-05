%{
Read in the subject file for all subjects. For the best clustering, generate
2D matrices for LE(t) and dFC(t) centroids by averaging according to 
the clustering. Next, average while multiplying the dFC(t) by occurrence.
Finally, find the Pearson's correlation between ROI timeseries to find
static FC and find the correlation with the total occurrence-weighted dFC(t).
Adapts code from: 
https://github.com/juanitacabral/LEiDA
https://github.com/trendscenter/gift
Output:
k_LE.h5 2D centroids for LE(t).
k_dFC.h5 2D centroids for dFC(t).
mean_dFC.h5 Occurrence-weighted mean dFC(t).
sFC.h5 sFC matrix.
corr_sFC.csv Correlation between the occurrence-weighted dFC(t) and sFC.
%}

%Define command line arguments.
function [] = r_LE_group_stateanalyze(subfile,k)
% subfile = 'r_full_submain.txt';
% k = '6';
disp(append('Doing ',subfile,' ',k));

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
if strcmp(subfile,'r_half_submain.txt')
    subgroup = 'half';
elseif strcmp(subfile,'r_full_submain.txt') 
    subgroup = 'full';    
end        
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                subgroup,'/',k,'/');
mean_dFCfile = append(outpath,'mean_dFC.h5');
mean_sFCfile = append(outpath,'sFC.h5');
corrfile = append(outpath,'corr_sFC.csv');

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

%Do LE weighted mean.
disp('Doing LE')
tic
VVT_mean=zeros(N_ba);
outfile = append(outpath,'k_LE.h5');
for c=1:N_Cl
    
    %Produce 2D matrix.
    VVT=V(c,:)'*V(c,:); 
    
    %Save 2D matrix.
    outkey = append('/',subgroup,'_',num2str(c));
    try
        h5create(outfile,outkey,size(VVT));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,VVT);
end
disp('LE saved.')
toc

%Do dFC weighted mean.
tic
nconn = (N_ba*(N_ba-1))/2;
dFC_mean = zeros(1,nconn);
k_mean = zeros(N_Cl,nconn);
for subin = 1:nsubj
        
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

        %Add to the corresponding running k unweighted mean.
        idx_c = IDX(t);
        k_mean(idx_c,:) = k_mean(idx_c,:) + (dFCmat(subt,:));

        %Add to the running total weighted mean.
        dFC_mean = dFC_mean + (dFCmat(subt,:)*F(idx_c));

        %Increment subject time.
        subt = subt + 1; 
    end
    delvars = {'dFCmat'};
    clear(delvars{:})
end
k_mean = k_mean/Tmax;

%Unvectorize and save the k means.
outfile = append(outpath,'k_dFC.h5');
for c=1:N_Cl
    
    %Produce 2D matrix. 
    k_2D = zeros(N_ba,N_ba);
    k_2D(triu(true(N_ba),1)) = k_mean(c,:);
    k_2D = k_2D + k_2D.';
    
    %Save 2D matrix.
    outkey = append('/',subgroup,'_',num2str(c));
    try
        h5create(outfile,outkey,size(k_2D));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,k_2D);
end

%Unvectorize and save the dFC total mean.
dFC_mean_2D = zeros(N_ba,N_ba);
dFC_mean_2D(triu(true(N_ba),1)) = dFC_mean;
dFC_mean_2D = dFC_mean_2D + dFC_mean_2D';
outkey = append('/',subgroup);
try
    h5create(mean_dFCfile,outkey,size(dFC_mean_2D));
catch
    disp('File and key created already.')
end
h5write(mean_dFCfile,outkey,dFC_mean_2D);
disp('dFC saved.')
toc

%Do sFC.
tic
tc = zeros(Tmax,N_ba);
for subin = 1:nsubj
    
    %Set subject indices.
    disp(append('Doing sFC: ',subjects(subin)))
    startin = (nwin*(subin-1)) + 1;
    endin = nwin*subin;
    
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
    
    %Concatenate and add.
    sub_tc = vertcat(sub_tc_L1,sub_tc_R1,sub_tc_L2,sub_tc_R2);
    tc(startin:endin,:) = sub_tc;
    delvars = {'sub_tc','sub_tc_L1','sub_tc_R1','sub_tc_L2','sub_tc_R2'};
    clear(delvars{:})
end    
sFC = corrcoef(tc);

%Save the sFC file.
outkey = append('/',subgroup);
try
    h5create(mean_sFCfile,outkey,size(sFC));
catch
    disp('File and key created already.')
end
h5write(mean_sFCfile,outkey,sFC);

%Find the dFC weighted mean correlation with sFC.
I_sup_diag=find(triu(ones(N_ba),1));
[cc,~]=corrcoef(sFC(I_sup_diag),dFC_mean_2D(I_sup_diag));

%Save the correlation.
corrmat = zeros(2,1);
corrmat(1,1) = cc(1,2);
corrmat = array2table(corrmat,'RowNames',{'dFC'},'VariableNames',{'r'});
writetable(corrmat,corrfile,'WriteRowNames',1)
disp('sFC saved.')
toc
end
