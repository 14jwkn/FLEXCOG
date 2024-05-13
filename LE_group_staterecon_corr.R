# For the specified k, find the correlations between reconfiguration variables
# and plot.
# Output:
# predcorr.csv Correlations between reconfiguration variables.
# predcorr_round.csv Correlations between reconfiguration variables, rounded.
# predcorr.jpg Correlations between reconfiguration variables, visualized.

#Set libraries.
library(tidyverse)

#Catches arguments.
args <- commandArgs(trailingOnly=T)
k <- args[1]

#Set base path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')

#Get predictor data.
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

#If there are NA in the matrix.
if (any(is.na(predmat))) {
  
  #Save the subjects with errors.
  remsubs <- predmat[rowSums(is.na(predmat))!=0,'Subject']
  outfile <- paste0(plscpath,'suberr.txt')
  write.table(remsubs,file=outfile,sep='\n',row.names=F,col.names=F)
  
  #Change matrices.
  predmat <- predmat[rowSums(is.na(predmat))==0,]
  keepsubs <- predmat[,'Subject']
  keepidx <- as.character(unlist(othercog[,'Subject'],use.names=F)) %in% as.character(unlist(keepsubs,use.names=F))
  othercog <- othercog[keepidx,]
  keepidx <- as.character(unlist(cogmat[,'Subject'],use.names=F)) %in% as.character(unlist(keepsubs,use.names=F))
  cogmat <- cogmat[keepidx,]
}

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

#Prepare data.
data_X <- data.frame(predmat)
rownames(data_X) <- data_X[,'Subject']
data_X <- within(data_X,rm('Subject'))

#Extract labels.
npred <- ncol(data_X)
predlabs <- colnames(data_X)

#Generate correlation table.
corrmat <- matrix(NA,nrow=npred,ncol=npred)
rownames(corrmat) <- predlabs
colnames(corrmat) <- predlabs
corrmat_r <- corrmat
for (cpred1 in predlabs) {
  for (cpred2 in predlabs) {
    corrmat[cpred1,cpred2] <- cor(data_X[,cpred1],data_X[,cpred2],method='spearman')
  }
}
corrmat_round <- round(corrmat,2)

#Save.
outfile <- paste0(outpath,'predcorr.csv')
write.csv(corrmat,outfile)
outfile <- paste0(outpath,'predcorr_round.csv')
write.csv(corrmat_round,outfile)

#Plot correlation matrix as a colored heatmap.
plotmat <- as.data.frame(corrmat) %>%
  rownames_to_column(var='var1') %>%
  gather(var2,Correlation,-var1) %>%
  mutate(var2=factor(var2,levels=unique(plotmat$var2)),  
         var1=factor(var1,levels=unique(plotmat$var1)))
p1 <- ggplot(plotmat,aes(var1,var2,fill=Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low='blue',mid='black',high='red',midpoint=0) +
  scale_y_discrete(limits=rev(levels(plotmat$var2))) +
  theme(axis.text.x=element_text(angle=90,hjust=1),
        axis.title=element_blank()) 
outfile <- paste0(outpath,'predcorr.jpg')
ggsave(outfile,p1,units='in',width=15,height=11)
