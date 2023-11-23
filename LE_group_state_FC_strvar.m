%{
For the specified k and participant, extract dFC(t) across and within states 
and within runs, calculate mean for FC strength and std for FC variability, 
average across runs.
Output:
dFC_strvar.h5 Contains dFC FC strength and variability across all states and within states.
%}

%Define command line arguments.
function [] = LE_group_state_FC_strvar(k,subject)

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
subfile = 'r_full_submain.txt';
subgroup = 'full';           
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                 subgroup,'/',k,'/',subject,'/group_strvar/');
dFCfile = append(outpath,'dFC_strvar.h5');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 

%Get the subject index and produce the start and end indices.
subin = find(strcmp(subjects,subject));
nwin = 4792;
startin = (nwin*(subin-1)) + 1;
endin = nwin*subin;
disp(append('Index: ',num2str(startin),'-',num2str(endin)))

%Read in the clustering and set up parameters.
inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/');
infile = append(inpath,'uni_subcent.h5');
inkey = append('/',subgroup,'/table');
fullclust = double(h5read(infile,inkey).values);
klabs = unique(fullclust)';
nk = str2num(k);
fullclust = fullclust(startin:endin);

%Read in the dFC and LE data.
Lfile1 = append('../outputs/r_LE_dFC/',clean,'/REST1_LR/',subject,...
                '/LE_dFC.h5');     
Rfile1 = append('../outputs/r_LE_dFC/',clean,'/REST1_RL/',subject,...
                '/LE_dFC.h5'); 
Lfile2 = append('../outputs/r_LE_dFC/',clean,'/REST2_LR/',subject,...
                '/LE_dFC.h5'); 
Rfile2 = append('../outputs/r_LE_dFC/',clean,'/REST2_RL/',subject,...
                '/LE_dFC.h5'); 

%dFC data.
inkey = '/LE_dFC';
dFCmatL1 = h5read(Lfile1,inkey);
dFCmatR1 = h5read(Rfile1,inkey);
dFCmatL2 = h5read(Lfile2,inkey);
dFCmatR2 = h5read(Rfile2,inkey);
dFCmat = vertcat(dFCmatL1,dFCmatR1,dFCmatL2,dFCmatR2);
delvars = {'dFCmatL1','dFCmatR1','dFCmatL2','dFCmatR2'};
clear(delvars{:})  
disp('Read input matrices.')

%Set up runs.
nrun = 4;   
runwin = size(fullclust,1)/nrun;

%Set dimensions.
ndFC = size(dFCmat,2);

%Produce the average STR and STD across runs.
strvar_dFC = zeros(2,ndFC);
dFCstr = zeros(nk,ndFC);
dFCvar = zeros(nk,ndFC);
for ridx = 1:nrun

    %Extract.
    runstart = (runwin*(ridx-1)) + 1;
    runend = (runwin*ridx);
    runclust = fullclust(runstart:runend);
    dFCrun = dFCmat(runstart:runend,:);
    
    %Find STR and STD and append.
    strvar_dFC(1,:) = strvar_dFC(1,:) + mean(dFCrun);
    strvar_dFC(2,:) = strvar_dFC(2,:) + std(dFCrun);
    
    %For each k.
    for kidx = 1:nk
        kval = klabs(kidx);

        %Find all rows that belong to the state.
        dFCstate = dFCrun(runclust == kval,:);

        %Find STR and STD and append.
        dFCstr(kidx,:) = dFCstr(kidx,:) + mean(dFCstate);
        dFCvar(kidx,:) = dFCvar(kidx,:) + std(dFCstate);
    end
end       

%Divide by the number of runs.
strvar_dFC = strvar_dFC / nrun;
dFCstr = dFCstr / nrun;
dFCvar = dFCvar / nrun;

%Attach labels.
dFCstr = horzcat(klabs',dFCstr);
dFCvar = horzcat(klabs',dFCvar);

%Save as tables.
outkey_list = {'strvar_dFC','dFCstr','dFCvar'};
outmat_list = {strvar_dFC,dFCstr,dFCvar};
outfile_list = {dFCfile,dFCfile,dFCfile};
nout = length(outkey_list);
for oidx=1:nout
    
    %Create file.
    outkey = append('/',outkey_list{oidx});
    outfile = outfile_list{oidx};
    outmat = outmat_list{oidx};
    try
        h5create(outfile,outkey,size(outmat));
    catch
        disp('File and key created already.')
    end
    
    %Write to file.
    h5write(outfile,outkey,outmat);
    % testout = h5read(outfile,outkey);
end
disp('Saved.')
end
