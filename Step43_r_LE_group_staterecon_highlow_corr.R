# For the subject group and the specified k, get the average for the specific target 
# set of network reconfiguration metrics and find the Pearson's correlations and
# associated p-values with edge-averaged dFC(t) strength and variability for 
# each state. 
# Output:
# recon_strvar.csv Correlation between strength and variability with the reconfiguration set average score.
# recon_strvar_P.csv P-value for the correlation between strength and variability with the reconfiguration set average score.
# recon_strvar_thres.csv Correlation between strength and variability with the reconfiguration set average score, thresholded p < 0.05.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]
k <- args[2]

#Set base path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/')
viewpath <- paste0(basepath,'allreconfig_subsets/')

#Set parameters.
nk <- as.numeric(k)
targlist <- c('g_jumppos','g_auto','g_simneg',
              'pr_126far','pr_auto','pr_simpos')
ntarg <- length(targlist)

#Read in strvar.
infile <- paste0(basepath,'dFC_strvar.csv')
mat_dFC <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject') 
strvar_labs <- c(colnames(mat_dFC)[2:ncol(mat_dFC)])

#Get confounds.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Define confounds.
con_list <- c('Age','Gender')
con_pred <- paste0(con_list,collapse='+')
con_mat <- othercog[,con_list]

#Join confounds and values.
val_pred <- strvar_labs
resmat <- cbind(mat_dFC[,val_pred],con_mat)
colnames(resmat) <- c(val_pred,con_list)

#For each value.
for (val_lab in val_pred) {
  
  #Produce residuals.
  inform <- as.formula(paste0(val_lab,'~',con_pred))
  respred <- resid(lm(inform,data=resmat))
  
  #Replace.
  mat_dFC[,val_lab] <- respred
}

#Read in reconfiguration summary values which have already been confound regressed.
nsub <- nrow(mat_dFC)
metvals <- matrix(NA,nrow=nsub,ncol=ntarg)
for (taidx in 1:ntarg) {
  targtype <- targlist[taidx]
  infile <- paste0(hlpath,targtype,'/',targtype,'_values.csv')
  cvals <- read_csv(infile,col_names=T) %>% 
    rename_at(1,~'Subject')
  metvals[,taidx] <- unlist(cvals[,'Value'])
}
colnames(metvals) <- targlist
metvals <- cbind(mat_dFC[,'Subject'],metvals)

#Combine data into one table.
total_mat <- metvals %>%
  inner_join(mat_dFC) 

#Correlate state values for dFC mean and std with reconfiguration summary values.
x1list <- targlist
x2list <- strvar_labs
collector <- matrix(NA,nrow=length(x1list),ncol=length(x2list))
rownames(collector) <- x1list
colnames(collector) <- x2list
p_collector <- collector
for (x1 in x1list) {
  for (x2 in x2list) {
    
    #Do specific ones.
    if ((((x1=='g_jumppos')|(x1=='pr_126far')) & grepl('mean_',x2,fixed=T)) |
        (((x1=='g_auto')|(x1=='pr_auto')) & grepl('std_',x2,fixed=T)) |
        (((x1=='g_simneg')|(x1=='pr_simpos')) & (grepl('mean_',x2,fixed=T)|grepl('std_',x2,fixed=T)))) {
      collector[x1,x2] <- cor(mat_dFC[,x1],mat_dFC[,x2],method='spearman')
      p_collector[x1,x2] <- cor.test(unlist(mat_dFC[,x1]),unlist(mat_dFC[,x2]),method='spearman')$p.value
    }
  }
}

#Get thresholded and rounded values.
collector_thres <- collector
collector_thres[p_collector >= 0.05] <- NA
collector_thres <- round(collector_thres,2)

#Save.
outpath <- paste0(viewpath,'strvar/')
dir.create(outpath,recursive=T)
outfile <- paste0(outpath,'recon_strvar.csv')
outtab <- t(as.data.frame(collector))
write.csv(outtab,outfile,row.names=T)
outfile <- paste0(outpath,'recon_strvar_P.csv')
outtab <- t(as.data.frame(p_collector))
write.csv(outtab,outfile,row.names=T)
outfile <- paste0(outpath,'recon_strvar_thres.csv')
outtab <- t(as.data.frame(collector_thres))
write.csv(outtab,outfile,row.names=T)
