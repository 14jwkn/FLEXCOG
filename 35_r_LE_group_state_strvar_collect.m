%{
For the subject group and the specified k, gather the average dFC(t) and
LE(t) strength and variability across edges for each state and across all
states from all subjects.
Output:
dFC_strvar.csv Average dFC(t) strength/variability values.
LE_strvar.csv Average LE(t) strength/variability values.
%}

%Define command line arguments.
function [] = 35_r_LE_group_state_strvar_collect(subfile,k)
disp(append('Doing: ',subfile,' ',k))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
if strcmp(subfile,'r_full_submain.txt') 
    subgroup = 'full';    
elseif strcmp(subfile,'r_half_submain.txt') 
    subgroup = 'half'; 
end                     
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                 subgroup,'/',k,'/');
file_dFC = append(outpath,'dFC_strvar.csv');
file_LE = append(outpath,'LE_strvar.csv');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
nsubj = size(subjects,1);

%Set labels.
klabs = (1:6);
inlabs = string({'full_mean' 'full_std'});
inlabs = [inlabs append('mean_',string(klabs))];
inlabs = [inlabs append('std_',string(klabs))];
nin = size(inlabs,2);

%For each state, save subject tables and group matrices.
nk = str2num(k);
tab_dFC = zeros(nsubj,nin);
tab_LE = zeros(nsubj,nin);
for sidx = 1:nsubj

    %Read in.
    inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                    clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/',subjects(sidx),'/');
    infile = append(inpath,'strvar.h5');
    inkey = append('/outvals');
    inmat = h5read(infile,inkey);

    %Extract and append.
    tab_dFC(sidx,:) = inmat(1,:);
    tab_LE(sidx,:) = inmat(2,:); 
end    

%Produce tables.
out_dFC = array2table(tab_dFC);
out_dFC.Properties.VariableNames = inlabs;
out_dFC.Properties.RowNames = subjects;
out_LE = array2table(tab_LE);
out_LE.Properties.VariableNames = inlabs;
out_LE.Properties.RowNames = subjects;

%Save the tables.
writetable(out_dFC,file_dFC,'WriteRowNames',true);
writetable(out_LE,file_LE,'WriteRowNames',true);
end
