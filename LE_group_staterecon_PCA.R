# For the specified k, conduct PCA on the cognition and network reconfiguration
# variables separately to check the relationship between the first PC and PLSC 
# dimension. Also plot first and second PC loadings to see which network reconfiguration 
# variables are similar and examine variance explained.
# Output:
# newlab,'_PLSC.jpg' First PLSC dimension loadings. 
# newlab,'_PCA.jpg' First PC loadings.
# PCA_PLSC_corr.csv Correlation between first PLSC and PC latent variables.
# PCA_PLSC_corr_r.csv Correlation between first PLSC and PC latent variables, rounded.
# ccomp,'_scatter.jpg' Scatter plot for a bivariate comparison.
# ccomp,'_qq.jpg' QQ plot for a bivariate comparison.
# binormtest_p.csv P-values from multiple normality tests for all bivariate comparisons.
# binormtest_pthres.csv P-values from multiple normality tests for all bivariate comparisons thresholded for significance.
# 2PC_load.jpg First and second PC loadings for network reconfiguration variables.

#Set libraries.
library(tidyverse)
library(ExPosition)
library(PTCA4CATA)
library(data4PCCAR)
library(ggplot2)
library(MVN)

# Helper Functions --------------------------------------------------------

#Plot loadings.
load_plot <- function(cload,norm.color,cwidth,cheight,outfile) {
  
  #Plot.
  collab <- cload
  collab[cload==cload] <- norm.color
  cplt <- PrettyBarPlot2(cload,
                         threshold=0, 
                         color4bar=collab,
                         font.size=3) +
    coord_cartesian(clip='off') +
    scale_y_continuous(expand=expansion(mult=c(0,0.1)))
  ggsave(outfile,cplt,units='in',width=cwidth,height=cheight)
}

#Normality testing function.
norm_checker <- function(part1,part2,complab,normpath) {
  
  #Set values.  
  ncomp <- length(complab)
  ctxtsize <- 11
  
  #Collect bivariate normality test outputs.
  testlabs <- c("mardia_skew","mardia_kurtosis","hz", "royston", "dh","energy")
  ntest <- length(testlabs)
  pr_norm <- matrix(NA,nrow=ncomp,ncol=ntest)
  rownames(pr_norm) <- complab
  colnames(pr_norm) <- testlabs
  for (coidx in 1:ncomp) {
    
    #Extract.
    ccomp <- complab[coidx]
    c1 <- part1[[coidx]]
    c2 <- part2[[coidx]]
    
    #Combine into dataframe.
    c1c2 <- data.frame('lab1'=c1,'lab2'=c2)
    c1c2_labs <- colnames(c1c2)
    
    #Plot scatter plot and save.
    cplt <- ggplot(aes(c1,c2),data=c1c2) +
      geom_point(size=1) +
      geom_smooth(method=lm,se=F) +
      labs(x=c1c2_labs[1],y=c1c2_labs[2]) +
      theme(axis.text.x=element_text(size=ctxtsize),
            axis.text.y=element_text(size=ctxtsize))
    outfile <- paste0(normpath,ccomp,'_scatter.jpg')
    ggsave(outfile,cplt,units='in',width=4,height=3)
    
    #Plot bivariate QQ plot and save.
    outfile <- paste0(normpath,ccomp,'_qq.jpg')
    jpeg(file=outfile)
    mvn(c1c2,multivariatePlot='qq')
    dev.off()
    
    #Do each test.
    for (ctest in testlabs) {
      
      #Whether indexing needs to be specific.
      if (grepl('mardia',ctest,fixed=T)) {
        
        #Skew or kurtosis value.
        if (ctest=='mardia_skew') {
          ctest_raw <- 'mardia'
          ctest_lab <- 'Mardia Skewness'
        } else if (ctest=='mardia_kurtosis') {
          ctest_raw <- 'mardia'
          ctest_lab <- 'Mardia Kurtosis'
        }
        
        #Do test.
        fitted <- mvn(c1c2,mvnTest=ctest_raw)$multivariateNormality
        fac <- fitted[fitted$Test==ctest_lab,'p value']
        cpr <- as.numeric(levels(fac)[fac])
      } else {
        
        #Do test.
        fitted <- mvn(c1c2,mvnTest=ctest)$multivariateNormality
        cpr <- fitted[1,'p value']
      }
      
      #Save it.
      pr_norm[ccomp,ctest] <- cpr
    }
  }
  
  #Threshold p-values and save.
  pr_norm_thres <- pr_norm
  pr_norm_thres[pr_norm < 0.05] <- NA
  pr_norm_thres <- signif(pr_norm_thres,2)
  outfile <- paste0(normpath,'binormtest_p.csv')
  write.csv(data.frame(pr_norm),outfile)
  outfile <- paste0(normpath,'binormtest_pthres.csv')
  write.csv(data.frame(pr_norm_thres),outfile)
}

# Driver ------------------------------------------------------------------

#Catches arguments.
args <- commandArgs(trailingOnly=T)
k <- args[1]

#Set base path, PLSC path, and PCA path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
plscpath <- paste0(basepath,'allreconfig/PLSC/')
pcapath <-  paste0(basepath,'allreconfig/PCA/')
dir.create(pcapath,recursive=T)

#Set parameters.
nk <- as.numeric(k)
theme_set(theme_gray())

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Isolate desired cognitive variables.
ourcog <- colnames(othercog)
ourcog <- c('Subject','PMAT24','PicVoc','ProcSpeed','CardSort','Flanker',
            'ListSort','PicSeq','ReadVoc','WordMem','LineOrient') 
cogmat <- othercog[,ourcog]
ncog <- ncol(cogmat) - 1

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

#Extract labels.
ncog <- ncol(data_Y)
npred <- ncol(data_X)
coglabs <- colnames(data_Y)
predlabs <- colnames(data_X)

#Retrieve first PLSC cognition latent variable and loadings.
infile <- paste0(plscpath,'perm/LVY_scores.csv')
lvy_PLSC <- read.csv(infile,row.names=1)
rowlab <- row.names(lvy_PLSC)
lvy_PLSC <- lvy_PLSC[,1]
names(lvy_PLSC) <- rowlab
infile <- paste0(plscpath,'boot/LVYall_load_full.csv')
loy_PLSC <- read.csv(infile,row.names=1)
rowlab <- row.names(loy_PLSC)
loy_PLSC <- loy_PLSC[,1]
names(loy_PLSC) <- rowlab

#Retrieve first PLSC reconfiguration latent variable and loadings.
infile <- paste0(plscpath,'perm/LVX_scores.csv')
lvx_PLSC <- read.csv(infile,row.names=1)
rowlab <- row.names(lvx_PLSC)
lvx_PLSC <- lvx_PLSC[,1]
names(lvx_PLSC) <- rowlab
infile <- paste0(plscpath,'boot/LVXall_load_full.csv')
lox_PLSC <- read.csv(infile,row.names=1)
rowlab <- row.names(lox_PLSC)
lox_PLSC <- lox_PLSC[,1]
names(lox_PLSC) <- rowlab

#Plot first PLSC loadings.
lolist <- list(loy_PLSC,
               lox_PLSC[grepl('occur',names(lox_PLSC))],
               lox_PLSC[grepl('dwell',names(lox_PLSC))],
               lox_PLSC[grepl('numtrans',names(lox_PLSC))],
               lox_PLSC[grepl('trans_',names(lox_PLSC))],
               lox_PLSC[grepl('jump',names(lox_PLSC))],
               lox_PLSC[grepl('sim',names(lox_PLSC))])
oldlabs <- c('cog','occur','dwell','numtrans','trans','jump','sim')
newlabs <- c('cog','occur','dwell','transnum','transpro','transdist','idio')
widthlist <- c(2,1.5,1.5,0.75,10,10,1.5)
heightlist <- c(3,3,3,3,3,3,3)
nset <- length(newlabs)
for (seidx in 1:nset) {
  
  #Extract.
  oldlab <- oldlabs[seidx]
  newlab <- newlabs[seidx]
  cload <- lolist[[seidx]]
  names(cload) <- gsub(oldlab,newlab,names(cload))
  cwidth <- widthlist[seidx]
  cheight <- heightlist[seidx]
  outfile <- paste0(pcapath,newlab,'_PLSC.jpg')
  
  #Plot.
  norm.color <- '#AE9F0F'
  load_plot(cload,norm.color,cwidth,cheight,outfile)
}

#Conduct PCA for cognition.
set.seed(12345)
data_XY <- data_Y
cogpca <- epPCA(data_XY,scale=T,center=T,graphs=F)
lvy_PCA <- -(cogpca$ExPosition.Data$fi)
rowlab <- row.names(lvy_PCA)
lvy_PCA <- lvy_PCA[,1]
names(lvy_PCA) <- rowlab
loy_PCA <- -(cogpca$ExPosition.Data$fj)
rowlab <- row.names(loy_PCA)
loy_PCA <- loy_PCA[,1]
names(loy_PCA) <- rowlab

#Conduct PCA for reconfiguration.
set.seed(12345)
data_XY <- data_X
reconpca <- epPCA(data_XY,scale=T,center=T,graphs=F)
lvx_PCA <- reconpca$Fixed.Data$ExPosition.Data$fi
rowlab <- row.names(lvx_PCA)
lvx_PCA <- lvx_PCA[,1]
names(lvx_PCA) <- rowlab
lox_PCA <- reconpca$Fixed.Data$ExPosition.Data$fj
rowlab <- row.names(lox_PCA)
lox_PCA <- lox_PCA[,1]
names(lox_PCA) <- rowlab

#Plot first PC loadings.
lolist <- list(loy_PCA,
               lox_PCA[grepl('occur',names(lox_PCA))],
               lox_PCA[grepl('dwell',names(lox_PCA))],
               lox_PCA[grepl('numtrans',names(lox_PCA))],
               lox_PCA[grepl('trans_',names(lox_PCA))],
               lox_PCA[grepl('jump',names(lox_PCA))],
               lox_PCA[grepl('sim',names(lox_PCA))])
oldlabs <- c('cog','occur','dwell','numtrans','trans','jump','sim')
newlabs <- c('cog','occur','dwell','transnum','transpro','transdist','idio')
widthlist <- c(2,1.5,1.5,0.75,10,10,1.5)
heightlist <- c(3,3,3,3,3,3,3)
nset <- length(newlabs)
for (seidx in 1:nset) {
  
  #Extract.
  oldlab <- oldlabs[seidx]
  newlab <- newlabs[seidx]
  cload <- lolist[[seidx]]
  names(cload) <- gsub(oldlab,newlab,names(cload))
  cwidth <- widthlist[seidx]
  cheight <- heightlist[seidx]
  outfile <- paste0(pcapath,newlab,'_PCA.jpg')
  
  #Plot.
  norm.color <- '#C10088'
  load_plot(cload,norm.color,cwidth,cheight,outfile2)
}

#Produce correlation matrix between first PLSC dimension and first PC.
corr_collect <- matrix(NA,2,2)
rownames(corr_collect) <- c('Pearson','Spearman')
colnames(corr_collect) <- c('Cognition','Reconfiguration')
corr_collect['Pearson','Cognition'] <- cor(lvy_PLSC,lvy_PCA,method='pearson')
corr_collect['Spearman','Cognition'] <- cor(lvy_PLSC,lvy_PCA,method='spearman')
corr_collect['Pearson','Reconfiguration'] <- cor(lvx_PLSC,lvx_PCA,method='pearson')
corr_collect['Spearman','Reconfiguration'] <- cor(lvx_PLSC,lvx_PCA,method='spearman')
corr_collect_r <- round(corr_collect,2)
outfile <- paste0(pcapath,'PCA_PLSC_corr.csv')
write.csv(corr_collect,outfile)
outfile <- paste0(pcapath,'PCA_PLSC_corr_r.csv')
write.csv(corr_collect_r,outfile)

#Check the normality for the comparisons.
part1 <- list(lvy_PLSC,lvx_PLSC)
part2 <- list(lvy_PCA,lvx_PCA)
complab <- c('plsy_pcay','plsx_pcax')
normpath <- paste0(pcapath,'normcheck/')
dir.create(normpath,recursive=T)
norm_checker(part1,part2,complab,normpath)

#Plot the PCs for reconfiguration to see if there are blocks of relationships, 
#coloring by dynamics, transition distance, and idiosyncrasy.
plotmat <- reconpca$Fixed.Data$ExPosition.Data$fj
plotlabs <- rownames(plotmat)
oldlab <- 'numtrans'
newlab <- 'transnum'
plotlabs <- gsub(oldlab,newlab,plotlabs)
oldlab <- 'trans_'
newlab <- 'transpro_'
plotlabs <- gsub(oldlab,newlab,plotlabs)
oldlab <- 'jump'
newlab <- 'transdist'
plotlabs <- gsub(oldlab,newlab,plotlabs)
oldlab <- 'sim'
newlab <- 'idio'
plotlabs <- gsub(oldlab,newlab,plotlabs)
plotlabs <- gsub('([0-9]+)t([0-9]+)','\\1-\\2',plotlabs)
rownames(plotmat) <- plotlabs
catlabs <- gsub("_.*", "", plotlabs)
col <- createColorVectorsByDesign(makeNominalData(as.matrix(catlabs)))
fig <- createFactorMap(plotmat,max.overlaps=1000,
                       col.points = col$oc,
                       col.labels = col$oc)
fig$zeMap + addArrows(plotmat, color = col$oc)
outfile <- paste0(pcapath,'2PC_load.jpg')
ggsave(outfile,units='in',width=14,height=14)
dev.off()

#Extract variance explained for first two PCs.
cogexp <- cogpca$Fixed.Data$ExPosition.Data$t
reconexp <- reconpca$Fixed.Data$ExPosition.Data$t
cog12 <- sum(cogexp[1:2])
recon12 <- sum(reconexp[1:2])
print('Variance Explained for PC 1 to 2:')
cat('Cognition:',cog12,'\n')
cat('Reconfiguration:',recon12,'\n')
