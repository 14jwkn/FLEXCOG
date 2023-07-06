# For the subject group and the specified k, find the Pearson's correlations between
# the idiosyncrasy metrics and the cognitive variables and the corresponding 
# uncorrected p-values.
# Output:
# olab,'.csv' Pearson's correlation, p-value, and p < 0.05 thresholded versions.

#Set libraries.
library(tidyverse)

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]
k <- args[2]

#Set base path and univariate path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/statesim/')
unipath <- paste0(basepath,'unicorr/')

#Set parameters.
nk <- as.numeric(k)

#Get predictor data.
infile <- paste0(basepath,'meansim.csv')
predmat <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
nstate <- ncol(predmat) - 1
statelabs <- colnames(predmat)[2:ncol(predmat)]

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Isolate desired cognitive variables.
ourcog <- colnames(othercog)
matchlist <- c('Subject','gPCA')
ourcog <- ourcog[ourcog %in% matchlist]
cogmat <- othercog[,ourcog]
ncog <- ncol(cogmat) - 1
coglabs <- ourcog[2:length(ourcog)]

#Remove desired confounds from predictors.
con_list <- c('Age','Gender')
con_pred <- paste0(con_list,collapse='+')

#Isolate confounds.
con_mat <- othercog[,con_list]

#Define values.
val_pred <- colnames(predmat)[-1]

#Join confounds and values.
resmat <- cbind(predmat,con_mat)

#For each value.
for (val_lab in val_pred) {
  
  #Produce residuals.
  inform <- as.formula(paste0(val_lab,'~',con_pred))
  respred <- resid(lm(inform,data=resmat))
  
  #Replace.
  predmat[,val_lab] <- respred
}

#Remove desired confounds from cognitive variables.
con_list <- c('Age','Gender')
con_pred <- paste0(con_list,collapse='+')

#Isolate confounds.
con_mat <- othercog[,con_list]

#Define values.
val_pred <- colnames(cogmat)[-1]

#Join confounds and values.
resmat <- cbind(cogmat,con_mat)

#For each value.
for (val_lab in val_pred) {
  
  #Produce residuals.
  inform <- as.formula(paste0(val_lab,'~',con_pred))
  respred <- resid(lm(inform,data=resmat))
  
  #Replace.
  cogmat[,val_lab] <- respred
}

#Prepare data.
data_X <- data.frame(predmat)
rownames(data_X) <- data_X[,'Subject']
data_X <- within(data_X,rm('Subject'))
data_Y <- data.frame(cogmat)
rownames(data_Y) <- data_Y[,'Subject']
data_Y <- within(data_Y,rm('Subject'))

# Correlations ------------------------------------------------------------

#Produce path.
outpath <- paste0(unipath)
dir.create(outpath,recursive=T)

#Set up matrix.
corrmat <- matrix(NA,nrow=nstate,ncol=ncog)
rownames(corrmat) <- statelabs
colnames(corrmat) <- coglabs
p_corrmat <- corrmat

#For each cognitive score.
for (cidx in 1:ncog) {
  
  #Extract.
  ycog <- coglabs[cidx]
  
  #For each state metric.
  for (sidx in 1:nstate) {
    
    #Extract.
    xstate <- statelabs[sidx]
    
    #Get the correlation and p-value.
    corrmat[xstate,ycog] <- cor(data_Y[,ycog],data_X[,xstate],method='spearman')
    p_corrmat[xstate,ycog] <- cor.test(data_Y[,ycog],data_X[,xstate],method='spearman')$p.value
  }
}

#Threshold.
pthres <- p_corrmat >= 0.05
pthres_p_corrmat <- p_corrmat
pthres_p_corrmat[pthres] <- NA
pthres_corrmat <- corrmat
pthres_corrmat[pthres] <- NA

#Save each file.
outlabs <- c('corr','pcorr','corr_thres','pcorr_thres')
nout <- length(outlabs)
outvars <- list(corrmat,p_corrmat,pthres_corrmat,pthres_p_corrmat)
for (oidx in 1:nout) {
  
  #Extract.
  olab <- outlabs[[oidx]]
  ovar <- outvars[[oidx]]
  
  #Save.
  outfile <- paste0(outpath,olab,'.csv')
  write.csv(ovar,outfile,row.names=T)
}
