# For the specified k, generate scatter and QQ plots for every bivariate comparison
# between sFC and dFC and do multiple normality tests. 
# Output:
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
        fitted <- mvn(c_XY,mvnTest=ctest_raw)$multivariateNormality
        fac <- fitted[fitted$Test==ctest_lab,'p value']
        cpr <- as.numeric(levels(fac)[fac])
      } else {
        
        #Do test.
        fitted <- mvn(c_XY,mvnTest=ctest)$multivariateNormality
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

#Set parameters.
nk <- as.numeric(k)
nroi <- 360
nconn <- ((nroi*(nroi-1))/2)

#Set output path.
subgroup <- 'full'
outpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                  subgroup,'/',k,'/dFC_sFC/')

#Append to comparison list.
part1 <- list()
part2 <- list()
complab <- c()

#Read in state-wise dFC and sFC.
for (kidx in 1:nk) {
  infile <- paste0(outpath,'dFC_sFC_',kidx,'_1D.csv')
  inmat <- read.csv(infile,row.names=1)
  part1[[kidx]] <- unlist(inmat[1,])
  part2[[kidx]] <- unlist(inmat[2,])
  complab <- c(complab,paste0('dFC_sFC_',as.character(kidx)))
}

#Read in full dFC and sFC.
infile <- paste0(outpath,'dFC_sFC_full_1D.csv')
inmat <- read.csv(infile,row.names=1)
part1[[(nk+1)]] <- inmat[1,]
part2[[(nk+1)]] <- inmat[2,]
complab <- c(complab,'dFC_sFC_full')

#Check every comparison for normality.
normpath <- paste0(outpath,'normcheck/')
dir.create(normpath,recursive=T)
norm_checker(part1,part2,complab,normpath)
print('Saved.')
