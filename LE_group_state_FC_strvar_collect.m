%{
For the specified k, collect the dFC FC strength and variability across all 
states and within states from all participants. Generate averaged 2D matrices from
vectorized forms.
Output:
dFC_strvar.h5 Contains dFC FC strength and variability across all states and within states from all participants and average 2D matrices.
%}

%Define command line arguments.
function [] = LE_group_state_FC_strvar_collect(k)

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
subfile = 'r_full_submain.txt';
subgroup = 'full';                
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                 subgroup,'/',k,'/group_strvar/');
out_dFCfile = append(outpath,'dFC_strvar.h5');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
nsubj = size(subjects,1);

%Read in one subject's data for parameters.
inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                subgroup,'/',k,'/',subjects(1),'/group_strvar/');
in_dFCfile = append(inpath,'dFC_strvar.h5');
inkey = append('/dFCstr');
dFCsub = h5read(in_dFCfile,inkey);
dFCdim = size(dFCsub,2) - 1;
klabs = dFCsub(:,1);

%Across states, produce subject tables. Set up table.
dFCstr = zeros(nsubj,dFCdim);
dFCvar = zeros(nsubj,dFCdim);

%For each subject.
for sidx = 1:nsubj
    disp('All');

    %Read in.
    inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                    subgroup,'/',k,'/',subjects(sidx),'/group_strvar/');
    in_dFCfile = append(inpath,'dFC_strvar.h5');
    inkey = append('/strvar_dFC');
    dFCsub = h5read(in_dFCfile,inkey);

    %Extract and append.
    dFCstr(sidx,:) = dFCsub(1,:);
    dFCvar(sidx,:) = dFCsub(2,:);
end    

%Save each subject table and averaged matrix.
outkey_list = {'dFCstr','dFCvar'};
outmat_list = {dFCstr,dFCvar};
nout = length(outkey_list);
for oidx=1:nout
 
    %Get matrix and average.
    cmat = outmat_list{oidx};
    avg_cmat = mean(cmat);
    
    %Make 2D and set outfile depending on dFC or LE.
    avg_2D = zeros(LEdim,LEdim);
    avg_2D(triu(true(LEdim),1)) = avg_cmat;
    avg_2D = avg_2D + avg_2D';
    outfile = out_dFCfile;
    
    %Save the subject table.
    ckey = append('/',outkey_list{oidx});
    outkey = ckey;
    outmat = cmat;
    try
        h5create(outfile,outkey,size(outmat));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,outmat);
    
    %Save the 2D matrix.
    outkey = append(ckey,'_2D');
    outmat = avg_2D;
    try
        h5create(outfile,outkey,size(outmat));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,outmat);
end

%For each state, repeat.
nk = str2num(k);
for kidx = 1:nk
    
    %Extract.
    disp(kidx);
    kval = klabs(kidx);
    
    %Set up table.
    dFCstr = zeros(nsubj,dFCdim);
    dFCvar = zeros(nsubj,dFCdim);

    %For each subject.
    for sidx = 1:nsubj

        %Read in and append.
        inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
        clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/',subjects(sidx),'/group_strvar/');
        in_dFCfile = append(inpath,'dFC_strvar.h5');
        
        %For dFC.
        inkey = append('/dFCstr');
        dFCsub = h5read(in_dFCfile,inkey);
        dFCstr(sidx,:) = dFCsub(dFCsub(:,1) == kval,2:(dFCdim+1));
        inkey = append('/dFCvar');
        dFCsub = h5read(in_dFCfile,inkey);
        dFCvar(sidx,:) = dFCsub(dFCsub(:,1) == kval,2:(dFCdim+1));
    end    

    %Save each subject table and averaged matrix.
    outkey_list = {'dFCstr','dFCvar'};
    outmat_list = {dFCstr,dFCvar};
    nout = length(outkey_list);
    for oidx=1:nout

        %Get matrix and average.
        cmat = outmat_list{oidx};
        avg_cmat = mean(cmat);

        %Make 2D and set outfile depending on dFC or LE.
        avg_2D = zeros(LEdim,LEdim);
        avg_2D(triu(true(LEdim),1)) = avg_cmat;
        avg_2D = avg_2D + avg_2D';
        outfile = out_dFCfile;

        %Save the subject table.
        ckey = append('/',outkey_list{oidx},'_k',num2str(kval));
        outkey = ckey;
        outmat = cmat;
        try
            h5create(outfile,outkey,size(outmat));
        catch
            disp('File and key created already.')
        end
        h5write(outfile,outkey,outmat);

        %Save the 2D matrix.
        outkey = append(ckey,'_2D');
        outmat = avg_2D;
        try
            h5create(outfile,outkey,size(outmat));
        catch
            disp('File and key created already.')
        end
        h5write(outfile,outkey,outmat);
    end
end
end
