# For the specified k, for pairs of summary scores relating to g or processing 
# speed, find the correlations and associated p-values. 
# Output: 
# olab,'.csv' Correlation, p-value, and p < 0.05 thresholded versions.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)

#Catches arguments.
args <- commandArgs(trailingOnly=T)
k <- args[1]

#Set base path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                  subgroup,'/',k,'/')
sumpath <-  paste0(basepath,'allreconfig_summary/')
outpath <- paste0(sumpath,'/unicorr/')
dir.create(outpath,recursive=T)

#Set parameters.
nk <- as.numeric(k)
targlist <- c('g_dynwithin_pos','g_dynwithin_neg',
              'g_auto','g_jumppos_between','g_jumpneg_between',
              'g_simneg',
              'pr_dyn2345','pr_dyn_exit1_pos','pr_dyn_exit1_neg','pr_dyn_exit6_pos',
              'pr_16auto',
              'pr_simpos')
ntarg <- length(targlist)

#Read in subjects and get length.
subjects <- scan('r_full_submain.txt',character(),quote='') 
nsub <- length(subjects)

#Read in reconfiguration summary values.
data_XY <- matrix(NA,nrow=nsub,ncol=ntarg)
for (taidx in 1:ntarg) {
  targtype <- targlist[taidx]
  infile <- paste0(sumpath,targtype,'_values.csv')
  cvals <- read_csv(infile,col_names=T) %>% 
    rename_at(1,~'Subject')
  data_XY[,taidx] <- unlist(cvals[,'Value'])
}
colnames(data_XY) <- targlist

#Correlate.
x1list <- targlist
x2list <- targlist
collector <- matrix(NA,nrow=length(x1list),ncol=length(x2list))
rownames(collector) <- x1list
colnames(collector) <- x2list
p_collector <- collector
lab_collector <- collector
for (x1 in x1list) {
  for (x2 in x2list) {
    collector[x1,x2] <- cor(data_XY[,x1],data_XY[,x2],method='spearman')
    p_collector[x1,x2] <- cor.test(unlist(data_XY[,x1]),unlist(data_XY[,x2]),method='spearman')$p.value
    lab_collector[x1,x2] <- paste0(x1,'_X_',x2)
  }
}

#Get thresholded and rounded values.
pthres <- p_collector >= 0.05
p_collector_thres <- p_collector
p_collector_thres[pthres] <- NA
collector_thres <- collector
collector_thres[pthres] <- NA
collector_thres <- round(collector_thres,2)

#Turn the p-values into a vector.
p_vec <- p_collector[upper.tri(p_collector,diag=F)]
names(p_vec) <- lab_collector[upper.tri(p_collector,diag=F)]

#Save.
outlabs <- c('corr','pcorr','corr_thres','pcorr_thres','pcorr_vec')
nout <- length(outlabs)
outvars <- list(collector,p_collector,collector_thres,p_collector_thres,p_vec)
for (oidx in 1:nout) {
  
  #Extract.
  olab <- outlabs[[oidx]]
  ovar <- outvars[[oidx]]
  
  #Save.
  outfile <- paste0(outpath,olab,'.csv')
  write.csv(ovar,outfile,row.names=T)
}
