# For the specified k, find the correlations between gPLSC, gPCA, gEFA, gCFA,
# and the 10 input cognitive variables.
# Output:
# cogcorr.csv Correlations between cognitive variables.
# cogcorr_round.csv Correlations between cognitive variables, rounded.

#Set libraries.
library(tidyverse)

#Catches arguments.
args <- commandArgs(trailingOnly=T)
k <- args[1]

#Get gPLSC.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                  subgroup,'/',k,'/')
plscpath <-  paste0(basepath,'allreconfig/PLSC/')
inpath <- paste0(plscpath,'perm/')
infile <- paste0(inpath,'LVY_scores.csv')
lvy <- read.csv(infile,row.names=1)
gPLSC <- lvy[,'LVY1']
names(gPLSC) <- rownames(lvy)

#Set outpath.
outpath <- paste0('../outputs/c_cognition/',subgroup,'/')
dir.create(outpath,recursive=T)

#Read subjects.
subfile <- 'r_full_submain.txt'
subjects <- scan(subfile,character(),quote='') 

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) 
othercog <- data.frame(othercog)
rownames(othercog) <- othercog[,'Subject']

#Isolate desired cognitive variables and subjects.
matchlist <- c('gPCA','gCFA','gEFA',
               'PMAT24','PicVoc','ProcSpeed','CardSort','Flanker',
               'ListSort','PicSeq','ReadVoc','WordMem','LineOrient')
cogmat <- othercog[subjects,matchlist]
coglabs <- colnames(cogmat)
ncog <- length(coglabs)

#Remove desired confounds from cognitive variables to match with gPLSC.
cogmat_r <- cogmat

#Define confounds.
con_list <- c('Age','Gender')
con_pred <- paste0(con_list,collapse='+')

#Isolate confounds.
con_mat <- othercog[,con_list]

#Define values.
val_pred <- colnames(cogmat_r)[-1]

#Join confounds and values.
resmat <- cbind(cogmat_r,con_mat)

#For each value.
for (val_lab in val_pred) {

  #Produce residuals.
  inform <- as.formula(paste0(val_lab,'~',con_pred))
  respred <- resid(lm(inform,data=resmat))

  #Replace.
  cogmat[,val_lab] <- respred
}

#Attach gPLSC.
cogmat$gPLSC <- gPLSC
coglabs <- c('gPLSC',coglabs)
ncog <- length(coglabs)

#Generate correlation table.
corrmat <- matrix(NA,nrow=ncog,ncol=ncog)
rownames(corrmat) <- coglabs
colnames(corrmat) <- coglabs
corrmat_r <- corrmat
for (ccog1 in coglabs) {
  for (ccog2 in coglabs) {
    corrmat[ccog1,ccog2] <- cor(cogmat[,ccog1],cogmat[,ccog2],method='spearman')
  }
}
corrmat_round <- round(corrmat,2)

#Save.
outfile <- paste0(outpath,'cogcorr.csv')
write.csv(corrmat,outfile)
outfile <- paste0(outpath,'cogcorr_round.csv')
write.csv(corrmat_round,outfile)
