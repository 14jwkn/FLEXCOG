%{
For the specified k, FDR-correct the p-values for the correlation between pairs of 
summary scores relating to g or processing speed.
Output:
pcorr_vec_fullfdr.csv FDR-corrected p-values.
%}

%Define command line arguments.
function [] = LE_group_staterecon_summary_corr_FDR(k)

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
basepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                 subgroup,'/',k,'/');
sumpath = append(basepath,'allreconfig_summary/');
setpath = append(sumpath,'/unicorr/');

%Extract labels.
targtypes = {'g_dynwithin_pos','g_dynwithin_neg',
             'g_auto','g_jumppos_between','g_jumpneg_between',
             'g_simneg',
             'pr_dyn2345','pr_dyn_exit1_pos','pr_dyn_exit1_neg','pr_dyn_exit6_pos',
             'pr_16auto',
             'pr_simpos'}
ntarg = size(targtypes,2);

%Produce vector and FDR correct.
infile = append(setpath,'pcorr_vec.csv');
pmat = readtable(infile,'ReadRowNames',1);
xlabs = string(pmat.Properties.VariableNames);
pvec = table2array(pmat);
allfdr = mafdr(pvec,'BHFDR',1);

%Save.
outfile = append(setpath,'pcorr_vec_fullfdr.csv');
outtab = array2table([xlabs' allfdr']);
outtab.Properties.VariableNames = {'Label','Value'};
writetable(outtab,outfile)
end    
