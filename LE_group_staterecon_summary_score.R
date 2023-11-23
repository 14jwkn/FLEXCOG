# For the specified k, average for the specific set of selected network reconfiguration metrics.
# Output:
# cset,'_values.csv' Average values for the reconfiguration metrics in this set.

#Set libraries.
library(tidyverse)

#Catches arguments.
k <- args[1]

#Set base path and output path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
outpath <-  paste0(basepath,'allreconfig_summary/')
dir.create(outpath,recursive=T)

#Set parameters.
nk <- as.numeric(k)
limval <- 1
if (k == '6') {
  subsets <- c('g_dynwithin_pos','g_dynwithin_neg',
               'g_auto','g_jumppos_between','g_jumpneg_between',
               'g_simneg',
               'pr_dyn2345','pr_dyn_exit1_pos','pr_dyn_exit1_neg','pr_dyn_exit6_pos',
               'pr_16auto',
               'pr_simpos')
  nsets <- length(subsets)
}

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

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

#Format full matrix.
fullmat <- predmat %>%
    mutate(Subject=as.character(Subject)) %>%
    mutate_if(is.character,as.factor)

# Extract Scores ----------------------------------------------------------

#Standardize. 
standmat <- within(fullmat,rm('Subject'))
standmat <- scale(standmat)

#Go through each subset.
for (setidx in 1:nsets) {
  
  #Extract.
  cset <- subsets[setidx]
  print(cset)
  if (k == '6') {
    if (cset == 'g_dynwithin_pos') {
      posvars <- c('dwell_2','dwell_3','trans_2t2','trans_3t3')
    } else if (cset == 'g_dynwithin_neg') {
      posvars <- c('trans_1t6','trans_3t6')
    } else if (cset == 'g_auto') {
      posvars <- c('jump_2t2','jump_3t3','jump_4t4','jump_5t5')
    } else if (cset == 'g_jumppos_between') {
      posvars <- c('jump_2t3','jump_2t4','jump_3t5','jump_4t2','jump_4t6','jump_5t3','jump_6t4')
    } else if (cset == 'g_jumpneg_between') {
      posvars <- c('jump_1t3','jump_1t4','jump_1t5','jump_3t1','jump_5t1',
                   'jump_2t5','jump_5t2')
    } else if (cset == 'g_simneg') {
      posvars <- c('sim_2','sim_3','sim_4','sim_5')
    } else if (cset == 'pr_dyn2345') {
      posvars <- c('occur_2','occur_3','occur_4','occur_5',
                   'dwell_3','dwell_5',
                   'trans_3t3','trans_5t5')
    } else if (cset == 'pr_dyn_exit1_pos') {
      posvars <- c('trans_1t3','trans_1t4','trans_1t5')
    } else if (cset == 'pr_dyn_exit1_neg') {
      posvars <- c('occur_1','dwell_1','trans_1t1',
                   'trans_2t1','trans_3t1','trans_4t1','trans_5t1','trans_6t1')
    } else if (cset == 'pr_dyn_exit6_pos') {
      posvars <- c('trans_6t3','trans_6t5')
    } else if (cset == 'pr_16auto') {
      posvars <- c('jump_1t1','jump_6t6')
    } else if (cset == 'pr_simpos') {
      posvars <- c('sim_1','sim_6')
    }
  }
  
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
