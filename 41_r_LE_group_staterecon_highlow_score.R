# For the subject group and the specified k, average for the specific set of
# selected network reconfiguration metrics and rank subjects according to it.
# Output:
# cset,'_values.csv' Average values for the reconfiguration metrics in this set.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(TExPosition)
library(PTCA4CATA)
library(data4PCCAR)

#Catches arguments.
subgroup <- args[1]
k <- args[2]

#Set base path and set path,
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/')
setpath <-  paste0(basepath,'allreconfig_subsets/')
dir.create(setpath,recursive=T)

#Set parameters.
nk <- as.numeric(k)
limval <- 1
if (k == '6') {
  subsets <- c('g_jumppos','g_auto','g_simneg','pr_auto','pr_126far','pr_simpos')
  nsets <- length(subsets)
}

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Isolate desired cognitive variables.
ourcog <- colnames(othercog)
matchlist <- paste0('Subject|CardSort|Flanker|ListSort|PicSeq|PV|PS|ReadEng|P24_CR|',
                    'IWRD_TOT|VSPLOT_TC')
removelist <- paste0('_Age')
ourcog <- str_subset(ourcog,matchlist)
ourcog <- {ourcog[!str_detect(ourcog,removelist)]}
cogmat <- othercog[,ourcog] 
ncog <- ncol(cogmat) - 1

#Add in gPCA.
cogmat[,'gPCA'] <- othercog[,'gPCA']

#Get predictor data from dynamics, state similarity, and state jump.
infile <- paste0(basepath,'avg_substat.csv')
dynmat <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
colnames(dynmat) <- gsub('sptrans_','trans_',colnames(dynmat))
infile <- paste0(basepath,'statesim/meansim.csv')
simmat <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
infile <- paste0(basepath,'statejump/transmean.csv')
colnames(simmat) <- gsub('State','sim',colnames(simmat))
jumpmat <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
colnames(jumpmat) <- c('Subject',paste0('jump_',colnames(jumpmat)[2:ncol(jumpmat)]))
predmat <- dynmat %>%
  inner_join(simmat) %>%
  inner_join(jumpmat)
colnames(predmat) <- gsub('-','t',colnames(predmat))

#Remove desired confounds from predictor variables.
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

#Generate full matrix after regressing out confounds.
fullmat <- predmat %>%
    mutate(Subject=as.character(Subject)) %>%
    mutate_if(is.character,as.factor)
fullmat <- cogmat %>%
    inner_join(fullmat) %>%
    rename('PMAT24'='P24_CR',
           'PicVoc'='PV',
           'ProcSpeed'='PS',
           'WordMem'='IWRD_TOT',
           'ReadVoc'='ReadEng',
           'LineOrient'='VSPLOT_TC')

#Remove gPCA.
cogmat <- within(cogmat,rm('gPCA'))

# Limit Subjects ----------------------------------------------------------

#Check mean and SD, standardize. 
standmat <- within(fullmat,rm('Subject'))
standmat <- scale(standmat)

#Find number to retain.
nsub <- nrow(fullmat)
nretain <- round(nsub*limval,0)

#Go through each subset.
for (setidx in 1:nsets) {
  
  #Extract.
  cset <- subsets[setidx]
  print(cset)
  if (k == '6') {
    if (cset == 'g_jumppos') {
      posvars <- c('jump_2t3','jump_2t4','jump_3t5','jump_4t2','jump_4t6','jump_5t3','jump_6t4')
    } else if (cset == 'g_auto') {
      posvars <- c('jump_2t2','jump_3t3','jump_4t4','jump_5t5')
    } else if (cset == 'g_simneg') {
      posvars <- c('sim_2','sim_3','sim_4','sim_5')
    } else if (cset == 'pr_auto') {
      posvars <- c('jump_1t1','jump_2t2','jump_4t4','jump_5t5','jump_6t6')
    } else if (cset == 'pr_126far') {
      posvars <- c('jump_1t6','jump_6t2')
    } else if (cset == 'pr_simpos') {
      posvars <- c('sim_1','sim_6')
    } 
  }
  
  #Set outpath.
  outpath <- paste0(setpath,cset,'/')
  print(outpath)
  dir.create(outpath,recursive=T)
  
  #Retrieve and average.
  cmat <- data.frame(standmat[,posvars])
  cmean <- apply(cmat,1,mean)

  #Attach and sort.
  cmean <- cbind(fullmat[,'Subject'],cmean)
  colnames(cmean) <- c('Subject','Value')
  
  #Save values.
  outfile <- paste0(outpath,cset,'_values.csv')
  write_csv(cmean,outfile)
}
