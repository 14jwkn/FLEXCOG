%{
For the subject group and the specified k, calculate the average dFC(t) and
LE(t) strength and variability across edges for each state and across all
states for each subject.
Output:
strvar.h5 Average dFC(t) and LE(t) strength/variability values for the subject.
%}

%Define command line arguments.
function [] = r_LE_group_state_strvar(subfile,k,subject)
disp(append('Doing: ',subfile,' ',k,' ',subject))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
if strcmp(subfile,'r_full_submain.txt') 
    subgroup = 'full';    
elseif strcmp(subfile,'r_half_submain.txt') 
    subgroup = 'half'; 
end        
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                 '/',subgroup,'/',k,'/',subject,'/');
outfile = append(outpath,'strvar.h5');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 

%Get the subject index and produce the start and end indices.
subin = find(strcmp(subjects,subject));
if ~strcmp(order,'both')
    nwin = 2396;
else
    nwin = 4792;
end
startin = (nwin*(subin-1)) + 1;
endin = nwin*subin;
disp(append('Index: ',num2str(startin),'-',num2str(endin)))

%Read in the clustering and set up parameters.
inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                subgroup,'/',k,'/');
infile = append(inpath,'uni_subcent.h5');
inkey = append('/',subgroup,'/table');
fullclust = double(h5read(infile,inkey).values);
klabs = unique(fullclust)';
nk = str2num(k);
fullclust = fullclust(startin:endin);

%Read in the LE and dFC data for two runs.
Lfile1 = append('../outputs/r_LE_dFC/REST1_LR/',subject,...
                '/LE_dFC.h5');     
Rfile1 = append('../outputs/r_LE_dFC/REST1_RL/',subject,...
                '/LE_dFC.h5'); 
Lfile2 = append('../outputs/r_LE_dFC/REST2_LR/',subject,...
                '/LE_dFC.h5'); 
Rfile2 = append('../outputs/r_LE_dFC/REST2_RL/',subject,...
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

%LE data.
inkey = '/LE_LE';
LEmatL1 = h5read(Lfile1,inkey);
LEmatR1 = h5read(Rfile1,inkey);
LEmatL2 = h5read(Lfile2,inkey);
LEmatR2 = h5read(Rfile2,inkey);
LEmat = vertcat(LEmatL1,LEmatR1,LEmatL2,LEmatR2);
delvars = {'LEmatL1','LEmatR1','LEmatL2','LEmatR2'};
clear(delvars{:}) 
disp('Read input matrices.')

%Set up runs.
nrun = 4;
runwin = size(fullclust,1)/nrun;

%Produce the full average absolute strength.
fullavg_dFC = mean(mean(abs(dFCmat),1));
fullavg_LE = mean(mean(abs(LEmat),1));

%Produce the full average STD.
fullstd_dFC = mean(std(dFCmat,1));
fullstd_LE = mean(std(LEmat,1));

%Produce the average absolute strength and STD for each state.
stateavg_dFC = zeros(nk,nrun);
stateavg_LE = zeros(nk,nrun);
statestd_dFC = zeros(nk,nrun);
statestd_LE = zeros(nk,nrun);
for kidx = 1:nk
    
    %Extract.
    kval = klabs(kidx);

    %For each run.
    for ridx = 1:nrun
        
        %Extract.
        runstart = (runwin*(ridx-1)) + 1;
        runend = (runwin*ridx);
        runclust = fullclust(runstart:runend);
        dFCrun = dFCmat(runstart:runend,:);
        LErun = LEmat(runstart:runend,:);
        
        %Find all rows that belong to the state.
        dFCstate = dFCrun(runclust == kval,:);
        LEstate = LErun(runclust == kval,:);
        
        %Find the average absolute strength and append.
        stateavg_dFC(kidx,ridx) = mean(mean(abs(dFCstate),1));
        stateavg_LE(kidx,ridx) = mean(mean(abs(LEstate),1));
        
        %Find the average STD and append.
        statestd_dFC(kidx,ridx) = mean(std(dFCstate,1));
        statestd_LE(kidx,ridx) = mean(std(LEstate,1));
    end       
end    

%Collect values.
outlabs = string({'full_mean' 'full_std'});
outlabs = [outlabs append('mean_',string(klabs))];
outlabs = [outlabs append('std_',string(klabs))];
nout = size(outlabs,2);

%Collect dFC.
out_dFC = zeros(1,nout);
out_dFC(1,1) = fullavg_dFC;
out_dFC(1,2) = fullstd_dFC;
out_dFC(1,3:(3+nk-1)) = mean(stateavg_dFC,2);
out_dFC(1,(3+nk):(nk*2+2)) = mean(statestd_dFC,2);

%Collect LE.
out_LE = zeros(1,nout);
out_LE(1,1) = fullavg_LE;
out_LE(1,2) = fullstd_LE;
out_LE(1,3:(3+nk-1)) = mean(stateavg_LE,2);
out_LE(1,(3+nk):(nk*2+2)) = mean(statestd_LE,2);

%Attach labels.
outtab = zeros(2,nout);
outtab(1,:) = out_dFC;
outtab(2,:) = out_LE;

%Save as tables.
valkey = append('/outvals');
try
    h5create(outfile,valkey,size(outtab));
catch
    disp('File and key created already.')
end
h5write(outfile,valkey,outtab);
disp('Saved.')
end
