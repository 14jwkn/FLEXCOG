# For the specified k, investigate the relationship between the main and distance-based LVs,
# between the main and frequency-based LVs, and the first and second LVs of the main,
# distance-based, and frequency-based PLSCs. Also investigate normality of comparisons.
# Output:
# sepcorr.csv Correlations between main LV1 with distance-based LV1, and main LV2 with distance-based LV2 and frequency-based LV1.
# corr_1v2.csv Correlations between the first two latent variables of main, distance-based, and frequency-based PLSC.
# ccomp,'_scatter.jpg' Scatter plot for a bivariate comparison.
# ccomp,'_qq.jpg' QQ plot for a bivariate comparison.
# binormtest_p.csv P-values from multiple normality tests for all bivariate comparisons.
# binormtest_pthres.csv P-values from multiple normality tests for all bivariate comparisons thresholded for significance.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(MVN)

# Helper Functions --------------------------------------------------------

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

#Set base path and output path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
outpath <-  paste0(basepath,'allreconfig/PLSC_freqdist/')
dir.create(outpath,recursive=T)

#Read in main LVs.
inpath <-  paste0(basepath,'allreconfig/PLSC/perm/')
infile <- paste0(inpath,'LVY_scores.csv')
mainLVY <- read.csv(infile,row.names=1)
infile <- paste0(inpath,'LVX_scores.csv')
mainLVX <- read.csv(infile,row.names=1)

#Read in frequency-based LVs.
inpath <-  paste0(basepath,'allreconfig/PLSC_freq/perm/')
infile <- paste0(inpath,'LVY_scores.csv')
tradLVY <- read.csv(infile,row.names=1)
infile <- paste0(inpath,'LVX_scores.csv')
tradLVX <- read.csv(infile,row.names=1)

#Read in distance-based LVs.
inpath <-  paste0(basepath,'allreconfig/PLSC_dist/perm/')
infile <- paste0(inpath,'LVY_scores.csv')
newLVY <- read.csv(infile,row.names=1)
infile <- paste0(inpath,'LVX_scores.csv')
newLVX <- read.csv(infile,row.names=1)

#Set correlation method.
corrmet <- 'spearman'

#Correlate main LV1 with distance-based LV1, and main LV2 with distance-based LV2
#and frequency-based LV1.
corrmat <- matrix(NA,3,2)
rownames(corrmat) <- c('m1_d1','m2_d2','m2_f1')
colnames(corrmat) <- c('Cognition','Reconfiguration')
corrmat[1,1] <- cor(mainLVY[,1],newLVY[,1],method=corrmet)
corrmat[1,2] <- cor(mainLVX[,1],newLVX[,1],method=corrmet)
corrmat[2,1] <- cor(mainLVY[,2],newLVY[,2],method=corrmet)
corrmat[2,2] <- cor(mainLVX[,2],newLVX[,2],method=corrmet)
corrmat[3,1] <- cor(mainLVY[,2],tradLVY[,1],method=corrmet)
corrmat[3,2] <- cor(mainLVX[,2],tradLVX[,1],method=corrmet)

#Save.
outfile <- paste0(outpath,'sepcorr.csv')
write.csv(corrmat,outfile)

#Find correlations between first two latent variables.
corrmat <- matrix(NA,3,4)
rownames(corrmat) <- c('Full','Frequency-Based','Distance-Based')
colnames(corrmat) <- c('Cognition','Reconfiguration','Cross_X1','Cross_Y1')
corrmat[1,1] <- cor(mainLVY[,1],mainLVY[,2],method=corrmet)
corrmat[2,1] <- cor(tradLVY[,1],tradLVY[,2],method=corrmet)
corrmat[3,1] <- cor(newLVY[,1],newLVY[,2],method=corrmet)
corrmat[1,2] <- cor(mainLVX[,1],mainLVX[,2],method=corrmet)
corrmat[2,2] <- cor(tradLVX[,1],tradLVX[,2],method=corrmet)
corrmat[3,2] <- cor(newLVX[,1],newLVX[,2],method=corrmet)
corrmat[1,3] <- cor(mainLVX[,1],mainLVY[,2],method=corrmet)
corrmat[2,3] <- cor(tradLVX[,1],tradLVY[,2],method=corrmet)
corrmat[3,3] <- cor(newLVX[,1],newLVY[,2],method=corrmet)
corrmat[1,4] <- cor(mainLVX[,2],mainLVY[,1],method=corrmet)
corrmat[2,4] <- cor(tradLVX[,2],tradLVY[,1],method=corrmet)
corrmat[3,4] <- cor(newLVX[,2],newLVY[,1],method=corrmet)

#Save.
outfile <- paste0(outpath,'corr_1v2.csv')
write.csv(corrmat,outfile)

#Check every comparison for normality.
part1 <- list(mainLVY[,1],mainLVX[,1],mainLVY[,2],mainLVX[,2],mainLVY[,2],mainLVX[,2],
              mainLVY[,1],tradLVY[,1],newLVY[,1],
              mainLVX[,1],tradLVX[,1],newLVX[,1],
              mainLVX[,1],tradLVX[,1],newLVX[,1],
              mainLVX[,2],tradLVX[,2],newLVX[,2])
part2 <- list(newLVY[,1],newLVX[,1],newLVY[,2],newLVX[,2],tradLVY[,1],tradLVX[,1],
              mainLVY[,2],tradLVY[,2],newLVY[,2],
              mainLVX[,2],tradLVX[,2],newLVX[,2],
              mainLVY[,2],tradLVY[,2],newLVY[,2],
              mainLVY[,1],tradLVY[,1],newLVY[,1])
complab <- c('m1y_d1y','m1x_d1x','m2y_d2y','m2x_d2x','m2y_f1y','m2x_f1x',
             'm1y_m2y','f1y_f2y','d1y_d2y',
             'm1x_m2x','f1x_f2x','d1x_d2x',
             'm1x_m2y','f1x_f2y','d1x_d2y',
             'm2x_m1y','f2x_f1y','d2x_d1y')
normpath <- paste0(outpath,'normcheck/')
dir.create(normpath,recursive=T)
norm_checker(part1,part2,complab,normpath)
print('Saved.')
