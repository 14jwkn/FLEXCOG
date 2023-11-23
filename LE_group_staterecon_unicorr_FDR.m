%{
For the specified k and cognitive variable, extract the frequency, transition
distance, and idiosyncrasy p-values. FDR-correct to generate new p-values
and threshold correlations with them.
Output:
pcorr.csv Original p-values from all metrics in one table.
pcorr_fullfdr.csv FDR-corrected p-values from all metrics in one table.
pcorr_fullfdr_thres.csv FDR-corrected p-values from all metrics thresholded p < 0.05.
corr.csv Original correlations from all metrics in one table.
corr_fullfdr_thres.csv  Original correlations thresholded p < 0.05.
%}

%Define command line arguments.
function [] = LE_group_staterecon_unicorr_FDR(k,ccog)

%Set base path and univariate path.
basepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                  subgroup,'/',k,'/');
dynpath = append(basepath,'dynstate/unicorr/',ccog,'/');
simpath = append(basepath,'statesim/unicorr/',ccog,'/');
jumppath = append(basepath,'statejump/unicorr/',ccog,'/');
outpath = append(basepath,'allreconfig/unicorr/',ccog,'/');
                
%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end               

%Set parameters.
nk = str2double(k);

%Read in dynamics p-values and correlations.
unipath = dynpath;
infile = append(unipath,'pcorr.csv');
pmat = readtable(infile);
xlabs = string(table2array(pmat(:,1)));
pmat = pmat(:,2:end);
pmat = pmat(:,ccog);
pmat.Properties.RowNames = xlabs;
dynmat = pmat;
infile = append(unipath,'corr.csv');
cmat = readtable(infile);
cmat = cmat(:,2:end);
cmat = cmat(:,ccog);
cmat.Properties.RowNames = xlabs;
cdyn = cmat;

%Read in similarity p-values and correlations.
unipath = simpath;
infile = append(unipath,'pcorr.csv');
pmat = readtable(infile);
xlabs = append('sim_',string(1:nk));
pmat = pmat(:,2:end);
pmat = pmat(:,ccog);
pmat.Properties.RowNames = xlabs;
simmat = pmat;
infile = append(unipath,'corr.csv');
cmat = readtable(infile);
cmat = cmat(:,2:end);
cmat = cmat(:,ccog);
cmat.Properties.RowNames = xlabs;
csim = cmat;

%Read in jump p-values and correlations.
unipath = jumppath;
infile = append(unipath,'pcorr.csv');
pmat = readtable(infile);
xlabs = string(table2array(pmat(:,1)));
% xlabs = strrep(xlabs,'sptrans_','jump_');
pmat = pmat(:,2:end);
pmat = pmat(:,ccog);
pmat.Properties.RowNames = xlabs;
jumpmat = pmat;
infile = append(unipath,'corr.csv');
cmat = readtable(infile);
cmat = cmat(:,2:end);
cmat = cmat(:,ccog);
cmat.Properties.RowNames = xlabs;
cjump = cmat;

%Concatenate.
allpmat = vertcat(dynmat,simmat,jumpmat);
xlabs = string(allpmat.Properties.RowNames);
pvec = table2array(allpmat);

%FDR correct for every test in the experiment.
allfdr = mafdr(pvec,'BHFDR',1);
allfdr_thres = allfdr;
allfdr_thres(allfdr >= 0.05) = NaN;

%Save.
outfile = append(outpath,'pcorr.csv');
outtab = array2table([xlabs pvec]);
outtab.Properties.VariableNames = {'Label',ccog};
writetable(outtab,outfile)
outfile = append(outpath,'pcorr_fullfdr.csv');
outtab = array2table([xlabs allfdr]);
outtab.Properties.VariableNames = {'Label',ccog};
writetable(outtab,outfile)
outfile = append(outpath,'pcorr_fullfdr_thres.csv');
outtab = array2table([xlabs allfdr_thres]);
outtab.Properties.VariableNames = {'Label',ccog};
writetable(outtab,outfile)

%Threshold original correlations.
allcmat = table2array(vertcat(cdyn,csim,cjump));
outfile = append(outpath,'corr.csv');
outtab = array2table([xlabs allcmat]);
outtab.Properties.VariableNames = {'Label',ccog};
writetable(outtab,outfile)
cmat_thres = allcmat;
cmat_thres(allfdr >= 0.05) = NaN;
outfile = append(outpath,'corr_fullfdr_thres.csv');
outtab = array2table([xlabs cmat_thres]);
outtab.Properties.VariableNames = {'Label',ccog};
writetable(outtab,outfile)
end
