%{
For the subject group and the specified k, FDR-correct the p-values for the
correlation between strength and variability with reconfiguration set 
average scores across all sets of interest.
Output:
reconvar_dFC_P_FDR.csv FDR-corrected p-values.
reconvar_dFC_P_FDR_thres.csv FDR-corrected p-values, thresholded p < 0.05.
reconvar_dFC_R_FDR_thres.csv Correlation, thresholded p < 0.05 for FDR-corrected p.
%}

%Define command line arguments.
function [] = 43_r_LE_group_staterecon_highlow_corr_FDR(subgroup,k)
disp(append('Doing: ',subgroup,' ',k))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
basepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                  clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/');
setpath = append(basepath,'allreconfig_subsets/');

%Set parameters.
nk = str2num(k);
nvals = nk*2;

%Extract labels.
targtypes = {'g_jumppos','g_auto','g_simneg',...
             'pr_126far','pr_auto','pr_simpos'};     
ntarg = size(targtypes,2);

%For each subset.
allxlabs = [];
allpvals = [];
allxvals = [];
for (tidx = 1:ntarg)
    
    %Extract.
    ctarg = string(targtypes(tidx));

    %Read in.
    targpath = append(setpath,ctarg,'/strvar/');
    infile = append(targpath,'reconvar_dFC_P.csv');
    pmat = readtable(infile);
    xlabs = string(table2array(pmat(5:end,1)));
    pvals = pmat(5:end,2:end);
    infile = append(targpath,'reconvar_dFC.csv');
    xvals = readtable(infile);
    xvals = xvals(5:end,2:end);
    
    %Append.
    allxlabs = vertcat(allxlabs,xlabs);
    allpvals = vertcat(allpvals,pvals);  
    allxvals = vertcat(allxvals,xvals);
end    

%Produce vector and FDR correct.
pvec = table2array(allpvals);
allfdr = mafdr(pvec,'BHFDR',1);
allfdr_thres = allfdr;
allfdr_thres(allfdr >= 0.05) = NaN;
allx_thres = table2array(allxvals);
allx_thres(allfdr >= 0.05) = NaN;

%For each subset.
for (tidx = 1:ntarg)
    
    %Extract.
    ctarg = string(targtypes(tidx));
    
    %Produce index.
    startidx = nvals*(tidx-1) + 1;
    endidx = nvals*tidx;

    %Save p.
    targpath = append(setpath,ctarg,'/strvar/');
    outfile = append(targpath,'reconvar_dFC_P_FDR.csv');
    outtab = array2table([allxlabs(startidx:endidx) allfdr(startidx:endidx)]);
    outtab.Properties.VariableNames = {'Label','Value'};
    writetable(outtab,outfile)
    
    %Save p thresholded.
    outfile = append(targpath,'reconvar_dFC_P_FDR_thres.csv');
    outtab = array2table([allxlabs(startidx:endidx) allfdr_thres(startidx:endidx)]);
    outtab.Properties.VariableNames = {'Label','Value'};
    writetable(outtab,outfile)
    
    %Save values thresholded.
    outfile = append(targpath,'reconvar_dFC_R_FDR_thres.csv');
    outtab = array2table([allxlabs(startidx:endidx) allx_thres(startidx:endidx)]);
    outtab.Properties.VariableNames = {'Label','Value'};  
    writetable(outtab,outfile)
end    
end    
