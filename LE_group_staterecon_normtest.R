# For the specified k, for pairs of cognitive and network reconfiguration variables, 
# conduct bivariate normality tests, plot scatterplots, and plot QQ plots.
# Output:
# binormtest_p.csv Bivariate normality test p-values for variable pairs.
# binormtest_pthres.csv Bivariate normality test p-values for significant latent variable pairs thresholded to significant ones.
# ccog,'_',cvar,'_scatter.jpg' Scatterplot for variable pairs.
# ccog,'_',cvar,'_scatterlim.jpg' Scatterplot for variable pairs, with limits for the specific network reconfiguration set for comparisons.
# ccog,'_',cvar,'_qq.jpg' Bivariate QQ plot for variable pairs.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(MVN)

#Catches arguments.
k <- args[1]

#Set base path and out path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
outpath <- paste0(basepath,'allreconfig/normtest/')
dir.create(outpath,recursive=T)

#Set parameters.
nk <- as.numeric(k)

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Isolate desired cognitive variables.
ourcog <- c('gPCA','gEFA','gCFA')
ncog <- length(ourcog)
cogmat <- othercog[,c('Subject',ourcog)]

#Get predictor data from dynamics, state similarity, and state jump.
infile <- paste0(basepath,'avg_substat.csv')
dynmat <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
colnames(dynmat) <- gsub('sptrans_','trans_',colnames(dynmat))
infile <- paste0(basepath,'statesim/',simstat,'sim_',statetype,'_',distmet,'.csv')
simmat <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
infile <- paste0(basepath,'statejump/',jumpstat,'_',statetype,'_',distmet,'.csv')
colnames(simmat) <- gsub('State','sim',colnames(simmat))
jumpmat <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
colnames(jumpmat) <- c('Subject',paste0('jump_',colnames(jumpmat)[2:ncol(jumpmat)]))
predmat <- dynmat %>%
  inner_join(simmat) %>%
  inner_join(jumpmat)
colnames(predmat) <- gsub('-','t',colnames(predmat))

#Remove confounds from predictors.
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

#Remove confounds from cognitive variables.
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

# Bivariate Normality Tests -----------------------------------------------

#Set up variables.
data_XY <- cbind(data_Y,data_X)

#Set up labels.
testlabs <- c('mardia_skew','mardia_kurtosis','hz','royston','dh','energy')
cogtest_labs <- c()
for (ccog in coglabs) {
  for (ctest in testlabs) {
    cogtest_labs <- c(cogtest_labs,paste0(ccog,'_',ctest))
  }
}
ncogtest <- length(cogtest_labs)

#Set up matrix to collect p-values.
pr_norm <- matrix(NA,nrow=npred,ncol=ncogtest)
rownames(pr_norm) <- predlabs
colnames(pr_norm) <- cogtest_labs
for (cpred in predlabs) {
  cat(cpred,'\n')
  for (ccog in coglabs) {
    for (ctest in testlabs) {
      ccogtest <- paste0(ccog,'_',ctest)
      c_XY <- data_XY[,c(ccog,cpred)]
      
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
      pr_norm[cpred,ccogtest] <- cpr
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

#Set up variables.
vtypes <- c('occur','dwell','numtrans','trans','sim','jump')
nvtypes <- length(vtypes)

#Go through each variable.
for (vtidx in 1:nvtypes) {
  ctype <- vtypes[vtidx]
  cat(ctype,'\n')
  
  #Extract all variables which match.
  cvarlist <- predlabs[grepl(ctype,predlabs)]
  if (ctype=='occur') {
    cvarlabs <- gsub('occur_','',cvarlist)
  } else if (ctype=='dwell') {
    cvarlabs <- gsub('dwell_','',cvarlist)
  } else if (ctype=='numtrans') {
    cvarlabs <- c('Transition Number')
  } else if (ctype=='trans') {
    cvarlist <- cvarlist[!grepl('numtrans',cvarlist)]
    cvarlabs <- gsub('trans_','',cvarlist)
    cvarlabs <- gsub('t','-',cvarlabs)
  } else if (ctype=='sim') {
    cvarlabs <- gsub('sim_','',cvarlist)
  } else if (ctype=='jump') {
    cvarlabs <- gsub('jump_','',cvarlist)
    cvarlabs <- gsub('t','-',cvarlabs)
  } 
  nvars <- length(cvarlist)
  
  #Find X limits for the set.
  xlim <- c()
  for (vidx in 1:nvars) {
    cvar <- cvarlist[vidx]
    cx <- data_XY[,cvar]
    xlim <- c(xlim,max(cx))
    xlim <- c(xlim,min(cx))
  }
  xmax <- max(xlim)
  xmin <- min(xlim)
  
  #Initiate plots and go through cognition.
  for (cidx in 1:ncog) {
    ccog <- coglabs[cidx]
    
    #Set up Y limits.
    ymax <- max(data_XY[,ccog])
    ymin <- min(data_XY[,ccog])
    ydiff <- ymax - ymin
    yinc <- ydiff*0.1
    ymax <- ymax + yinc
    ymin <- ymin - yinc
    
    #Go through predictors.
    for (vidx in 1:nvars) {
      cvar <- cvarlist[vidx]
      clab <- cvarlabs[vidx]
      
      #Generate scatterplot and save.
      cwidth <- 10
      cheight <- 10
      outfile <- paste0(outpath,ccog,'_',cvar,'_scatter.jpg')
      scatterp <- ggplot(data_XY,aes(!!!ensyms(x=cvar,y=ccog))) +
        geom_point() +
        geom_smooth(method='lm',se=F,fullrange=T,color='black') +
        xlab(cvar) +
        ylab(ccog) +
        theme(plot.title=element_blank(),
              plot.margin=grid::unit(c(1,1,1,1),'mm')) 
      ggsave(outfile,scatterp,
             units='in',width=cwidth,height=cheight)
      
      #Generate scatterplot with limits and save.
      outfile <- paste0(outpath,ccog,'_',cvar,'_scatterlim.jpg')
      scatterp_lim <- ggplot(data_XY,aes(!!!ensyms(x=cvar,y=ccog))) +
        geom_point() +
        geom_smooth(method='lm',se=F,fullrange=T,color='black') +
        scale_x_continuous(limits=c(xmin,xmax)) +   
        scale_y_continuous(limits=c(ymin,ymax)) +
        xlab(cvar) +
        ylab(ccog) +
        theme(plot.title=element_blank(),
          plot.margin=grid::unit(c(1,1,1,1),'mm')) 
      ggsave(outfile,scatterp_lim,
             units='in',width=cwidth,height=cheight)

      #Generate bivariate QQ plot and save.
      outfile <- paste0(outpath,ccog,'_',cvar,'_qq.jpg')
      jpeg(file=outfile)
      c_XY <- data_XY[,c(ccog,cvar)]
      mvn(c_XY,multivariatePlot='qq')
      dev.off()
    }
  }
}
