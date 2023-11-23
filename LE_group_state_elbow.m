%{
Use the group clustering and the centroids for each
k to produce the cluster validity indices to produce the elbow plot, where
lower values are better.
Adapts code from: 
https://github.com/trendscenter/gift
Output:
metric_group.h5 Contains the CVI values for each k for the group.
%}

%Define command line arguments.
function [] = LE_group_state_elbow()

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O.
subgroup = 'full';
kkey = append('/',subgroup);
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                 subgroup,'/');
groupfile = append(outpath,'metric_group.h5');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
nsubj = size(subjects,1);

%Set ks.
klist = 2:12;
lenk = size(klist,2);

%Read in best iterations.
best_file = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/',...
                  'group/best_iter/',subgroup,'/best_iter.csv');
best_clust = readtable(best_file);

%Produce clusterings.
clusts = {};
for i = 1:size(klist,2)
    curr_k = num2str(klist(i));
    inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                    subgroup,'/',curr_k,'/');
    best_lab = curr_k;
    iteration = num2str(best_clust{:,best_lab}); 
    infile = append(inpath,'subclust_',iteration,'.h5');
    clust = h5read(infile,kkey);
    clusts{i} = clust;
end    

%Produce centroids.
cents = {};
for i = 1:lenk
    curr_k = num2str(klist(i));
    inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                    subgroup,'/',curr_k,'/');
    best_lab = curr_k;
    iteration = num2str(best_clust{:,best_lab});                
    centfile = append(inpath,'subcent_',iteration,'.h5');
    cents{i} = h5read(centfile,kkey);
end    

%Do elbow metric finder for WCS and CVI.
results = elbow_k(dFCfull,klist,clusts,cents);
cvi_groupscore = results.cvi;
cvi_groupk = results.c_k;

%Package the group scores.
group_out = [[cvi_groupk repelem(0,lenk-1)];...
            cvi_groupscore;...
            ];
disp('Group done.')
        
%Save.
outkey = kkey;
try
    h5create(groupfile,outkey,size(group_out));
catch
    disp('File and key created already.')
end

%Write out outputs.
h5write(groupfile,outkey,group_out)
disp('Saved.')
end

