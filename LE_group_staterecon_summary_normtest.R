# For the specified k, for pairs of summary scores relating to g or processing 
# speed, conduct bivariate normality tests, plot scatterplots, and plot QQ plots.
# Output:
# binormtest_p.csv Bivariate normality test p-values for variable pairs.
# binormtest_pthres.csv Bivariate normality test p-values for significant latent variable pairs thresholded to significant ones.
# clab,'_scatter.jpg' Scatterplot for variable pairs.
# clab,'_scatterlim.jpg' Scatterplot for variable pairs, with limits for the specific network reconfiguration set for comparisons.
# clab,'_qq.jpg' Bivariate QQ plot for variable pairs.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(MVN)

#Catches arguments.
args = commandArgs(trailingOnly=T)
k <- args[1]
subset <- args[2]

#Set base path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
sumpath <-  paste0(basepath,'allreconfig_summary/')
outpath <- paste0(setpath,'/normtest/',subset,'/')
dir.create(outpath,recursive=T)

#Set parameters.
nk <- as.numeric(k)
if (subset == 'g') {
  targlist <- c('g_dynwithin_pos','g_dynwithin_neg',
                'g_auto','g_jumppos_between','g_jumpneg_between',
                'g_simneg')
} else if (subset == 'pr') {
  targlist <- c('pr_dyn2345','pr_dyn_exit1_pos','pr_dyn_exit1_neg','pr_dyn_exit6_pos',
                'pr_16auto',
                'pr_simpos')
}
ntarg <- length(targlist)

#Read in subjects and get length.
subjects <- scan('r_full_submain.txt',character(),quote='') 
nsub <- length(subjects)

#Read in reconfiguration summary values.
data_XY <- matrix(NA,nrow=nsub,ncol=ntarg)
for (taidx in 1:ntarg) {
  targtype <- targlist[taidx]
  infile <- paste0(sumpath,targtype,'_values.csv')
  cvals <- read_csv(infile,col_names=T) %>% 
    rename_at(1,~'Subject')
  data_XY[,taidx] <- unlist(cvals[,'Value'])
}
colnames(data_XY) <- targlist

# Bivariate Normality Tests -----------------------------------------------

#Set up labels.
testlabs <- c('mardia_skew','mardia_kurtosis','hz','royston','dh','energy')
ntest <- length(testlabs)
sumlabs <- c()
for (cv in targlist) {
  for (cv2 in targlist) {
    clab <- paste0(cv,'_',cv2)
    if ((cv != cv2) & !(clab %in% sumlabs)) {
      sumlabs <- c(sumlabs,clab)
    }
  }
}
nsum <- length(sumlabs)

#Set up matrix to collect p-values.
pr_norm <- matrix(NA,nrow=nsum,ncol=ntest)
rownames(pr_norm) <- sumlabs
colnames(pr_norm) <- testlabs
for (cv in targlist) {
  for (cv2 in targlist) {
    clab <- paste0(cv,'_',cv2)
    if ((cv != cv2) & !(clab %in% sumlabs)) {
      c_XY <- data_XY[,c(cv,cv2)]
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
        pr_norm[clab,ctest] <- cpr
      }
    }
  }
}

#Threshold p-values.
pr_norm_thres <- pr_norm
pr_norm_thres[pr_norm < 0.05] <- NA
pr_norm_thres <- signif(pr_norm_thres,2)

#Save the matrices.
outfile <- paste0(outpath,'binormtest_p.csv')
write.csv(data.frame(pr_norm),outfile)
outfile <- paste0(outpath,'binormtest_pthres.csv')
write.csv(data.frame(pr_norm_thres),outfile)

# Bivariate Diagnostic Plots ----------------------------------------------

#Set theme.
theme_set(theme_bw())

#Get maximum and minimum values for the plots.
maxval <- max(max(data_XY[,targlist]))
minval <- min(min(data_XY[,targlist])) 
vinc <- (maxval - minval)*0.1
maxval <- maxval + vinc
minval <- minval + vinc

#Go through each variable pair.
for (cv in targlist) {
  for (cv2 in targlist) {
    clab <- paste0(cv,'_',cv2)
    if ((cv != cv2) & !(clab %in% sumlabs)) {
      c_XY <- data_XY[,c(cv,cv2)]
      
      #Generate scatterplot and save.
      cwidth <- 10
      cheight <- 10
      outfile <- paste0(outpath,clab,'_scatter.jpg')
      scatterp <- ggplot(data_XY,aes(!!!ensyms(x=cv,y=cv2))) +
        geom_point() +
        geom_smooth(method='lm',se=F,fullrange=T,color='black') +
        xlab(cv) +
        ylab(cv2) +
        theme(plot.title=element_blank(),
              plot.margin=grid::unit(c(1,1,1,1),'mm')) 
      ggsave(outfile,scatterp,
             units='in',width=cwidth,height=cheight)
      
      #Generate scatterplot with limits and save.
      outfile <- paste0(outpath,clab,'_scatterlim.jpg')
      scatterp_lim <- ggplot(data_XY,aes(!!!ensyms(x=cv,y=cv2))) +
        geom_point() +
        geom_smooth(method='lm',se=F,fullrange=T,color='black') +
        scale_x_continuous(limits=c(minval,maxval)) +   
        scale_y_continuous(limits=c(minval,maxval)) +
        xlab(cv) +
        ylab(cv2) +
        theme(plot.title=element_blank(),
              plot.margin=grid::unit(c(1,1,1,1),'mm')) 
      ggsave(outfile,scatterp_lim,
             units='in',width=cwidth,height=cheight)
      
      #Generate bivariate QQ plot and save.
      outfile <- paste0(outpath,clab,'_qq.jpg')
      jpeg(file=outfile)
      mvn(c_XY,multivariatePlot='qq')
      dev.off()
    }
  }
}
