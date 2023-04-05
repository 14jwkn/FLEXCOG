# For the subject group and the specified k, get the average for the specific target 
# set of network reconfiguration metrics and find the Pearson's correlations and
# associated p-values with edge-averaged dFC(t) strength and variability for 
# each state. 
# Output:
# reconvar_dFC.csv Correlation between strength and variability with the reconfiguration set average score.
# reconvar_dFC_P.csv P-value for the correlation between strength and variability with the reconfiguration set average score.
# reconvar_dFC_thres.csv Correlation between strength and variability with the reconfiguration set average score, thresholded p < 0.05.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(scales)

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]
k <- args[2]
targtype <- args[3]

#Set base path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/')
viewpath <-  paste0(basepath,'allreconfig_subsets/')

#Set parameters.
nk <- as.numeric(k)

#Read in data.
infile <- paste0(viewpath,targtype,'_values.csv')
metvals <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
met_labs <- c(colnames(metvals)[2:ncol(metvals)])
infile <- paste0(basepath,'dFC_strvar.csv')
mat_dFC <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject') 
infile <- paste0(basepath,'LE_strvar.csv')
mat_LE <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject') 
strvar_labs <- c(colnames(mat_dFC)[2:ncol(mat_dFC)])
mat_dFC <- mat_dFC %>%
  inner_join(metvals) 
mat_LE <- mat_LE %>%
  inner_join(metvals) 

#Correlate full values and state values for mean and std for dFC with reconfig mean.
x1list <- strvar_labs
x2list <- met_labs
r_reconvar_dFC <- matrix(NA,nrow=length(x1list),ncol=length(x2list))
rownames(r_reconvar_dFC) <- x1list
colnames(r_reconvar_dFC) <- x2list
p_reconvar_dFC <- r_reconvar_dFC
for (x1 in x1list) {
  for (x2 in x2list) {
    r_reconvar_dFC[x1,x2] <- cor(mat_dFC[,x1],mat_dFC[,x2])
    p_reconvar_dFC[x1,x2] <- cor.test(unlist(mat_dFC[,x1]),unlist(mat_dFC[,x2]))$p.value
  }
}
r_reconvar_dFC_thres <- r_reconvar_dFC
r_reconvar_dFC_thres[p_reconvar_dFC >= 0.05] <- NA

#Round thresholded viewings.
r_reconvar_dFC_thres <- round(r_reconvar_dFC_thres,2)
r_reconvar_LE_thres <- round(r_reconvar_LE_thres,2)
r_gall_dFC_thres <- round(r_gall_dFC_thres,2)
r_gall_LE_thres <- round(r_gall_LE_thres,2)

#Save reconfiguration vs strength and variability.
outpath <- paste0(viewpath,'strvar/')
dir.create(outpath,recursive=T)
outfile <- paste0(outpath,'reconvar_dFC.csv')
outtab <- as.data.frame(r_reconvar_dFC)
write.csv(outtab,outfile,row.names=T)
outfile <- paste0(outpath,'reconvar_dFC_P.csv')
outtab <- as.data.frame(p_reconvar_dFC)
write.csv(outtab,outfile,row.names=T)
outfile <- paste0(outpath,'reconvar_dFC_thres.csv')
outtab <- as.data.frame(r_reconvar_dFC_thres)
write.csv(outtab,outfile,row.names=T)
