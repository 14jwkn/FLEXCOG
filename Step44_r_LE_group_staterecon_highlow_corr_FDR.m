%{
For the subject group and the specified k, FDR-correct the p-values for the
correlation between strength and variability with reconfiguration set 
average scores across all sets of interest.
Output:
recon_strvar_FDR_P.csv FDR-corrected p-values.
recon_strvar_FDR_P_thres.csv FDR-corrected p-values, thresholded p < 0.05.
recon_strvar_FDR_thres.csv Correlation, thresholded p < 0.05 for FDR-corrected p.
%}

%Define command line arguments.
function [] = Step44_r_LE_group_staterecon_highlow_corr_FDR(subgroup,k)
disp(append('Doing: ',subgroup,' ',k))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
basepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                  clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/');
setpath = append(basepath,'allreconfig_subsets/');

%Extract labels.
targtypes = {'g_jumppos','g_auto','g_simneg',...
             'pr_126far','pr_auto','pr_simpos'};     
ntarg = size(targtypes,2);

%For each subset.
allxlabs = [];
allpvals = [];
allxvals = [];
for tidx = 1:ntarg
    
    %Extract.
    ctarg = string(targtypes(tidx));

    %Read in.
    infile = append(setpath,'recon_strvar_P.csv');
    pmat = readtable(infile,'ReadRowNames',1);
    xlabs = string(pmat.Properties.RowNames);
    infile = append(setpath,'recon_strvar.csv');
    xmat = readtable(infile,'ReadRowNames',1);

    %Isolate based on label.
    if strcmp(ctarg,'g_jumppos')||strcmp(ctarg,'pr_126far')
        strvar_labs = string({'mean_'});
    elseif strcmp(ctarg,'g_auto')||strcmp(ctarg,'pr_auto')
        strvar_labs = string({'std_'}); 
    elseif strcmp(ctarg,'g_simneg')||strcmp(ctarg,'pr_simpos')
        strvar_labs = string({'mean_','std_'}); 
    end       
    xlabs = xlabs(contains(xlabs,strvar_labs));
    pvals = table2array(pmat(xlabs,ctarg));
    xvals = table2array(xmat(xlabs,ctarg));

    %Append.
    xlabs = append(ctarg,'_',xlabs);
    allxlabs = vertcat(allxlabs,xlabs);
    allpvals = vertcat(allpvals,pvals);  
    allxvals = vertcat(allxvals,xvals);
end    

%Produce vector and FDR correct.
pvec = allpvals;
allfdr = mafdr(pvec,'BHFDR',1);
allfdr_thres = allfdr;
allfdr_thres(allfdr >= 0.05) = NaN;
allx_thres = allxvals;
allx_thres(allfdr >= 0.05) = NaN;

%Save.
outpath = setpath;
outfile = append(outpath,'recon_strvar_FDR_P.csv');
outtab = array2table([allxlabs allfdr]);
outtab.Properties.VariableNames = {'Label','Value'};
writetable(outtab,outfile)
outfile = append(outpath,'recon_strvar_FDR_P_thres.csv');
outtab = array2table([allxlabs allfdr_thres]);
outtab.Properties.VariableNames = {'Label','Value'};
writetable(outtab,outfile)
outfile = append(outpath,sumlabel,'recon_strvar_FDR_thres.csv');
outtab = array2table([allxlabs allx_thres]);
outtab.Properties.VariableNames = {'Label','Value'};
writetable(outtab,outfile)   
end    
