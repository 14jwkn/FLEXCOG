# For the specified k, conduct PLSC where X = network
# reconfiguration metrics and Y = 10 cognitive tests of interest. Do PLSC to find
# direction of loadings to flip, conduct permutation testing on singular values,
# investigate bivariate normality and outliers for significant latent variable pairs and correlate, 
# conduct bootstrap resampling on loadings, and investigate reproducibility of
# singular values and loadings.
# Adapts code from:
# https://github.com/McIntosh-Lab/PLS-CCA
# Output:
# permsig.csv Table containing singular value, covariance explained, and permutation p-values.
# permsig.jpg Plot indicating which dimensions are significant.
# 'baseLVY',llab,'.jpg' Bar plot of Y cognitive loadings with default directions.
# yflip.csv Table indicating whether to flip the dimension for interpretation.
# 'LV',as.character(i),'_corr.jpg' Scatterplot for significant latent variable pairs.
# 'LV',as.character(i),'_qq.jpg' Bivariate QQ plot for significant latent variable pairs.
# binormtest_p.csv Bivariate normality test p-values for significant latent variable pairs.
# binormtest_pthres.csv Bivariate normality test p-values for significant latent variable pairs thresholded to significant ones.
# LVX_scores.csv Significant X latent variable scores.
# LVY_scores.csv Significant Y latent variable scores.
# 'LVY',llab,'_load.jpg' Bar plot of Y cognitive loadings where stable loadings are colored for one dimension.
# 'LVX',llab,'_load_',namedyn,'.jpg' Bar plot of X reconfiguration loadings where stable loadings are colored for one dimension.
# LVYall_load.jpg Bar of Y cognitive loadings for all dimensions with stable loadings colored.
# LVXall_load.jpg Bar of X reconfiguration loadings for all dimensions with stable loadings colored.
# LVXall_load.csv X reconfiguration loadings thresholded by bootstrap ratio.
# LVYall_load.csv Y cognitive loadings thresholded by bootstrap ratio.
# LVXall_load_full.csv X reconfiguration loadings.
# LVYall_load_full.csv Y cognitive loadings.
# 'LVY',llab,'_bootrat.jpg' Bar plot of Y cognitive bootstrap ratios where stable loadings are colored for one dimension.
# 'LVX',llab,'_bootrat_',namedyn,'.jpg' Bar plot of X reconfiguration bootstrap ratios where stable loadings are colored for one dimension.
# LVYall_bootrat.jpg Bar of Y cognitive bootstrap ratios for all dimensions with stable loadings colored.
# LVXall_bootrat.jpg Bar of X reconfiguration bootstrap ratios for all dimensions with stable loadings colored.
# LVXall_bootrat.csv X reconfiguration bootstrap ratios thresholded by bootstrap ratio.
# LVYall_bootrat.csv Y cognitive bootstrap ratios thresholded by bootstrap ratio.
# LVXall_bootrat_full.csv X reconfiguration bootstrap ratios.
# LVYall_bootrat_full.csv Y cognitive bootstrap ratios.
# olab,'.csv' Mean and Z-score reproducibility values for singular values and loadings.
# svrep.jpg Box plot of singular value reproducibility across split-halfs colored based on Z-score.
# lvxrep.jpg Box plot of X loading reproducibility across split-halfs colored based on Z-score.
# lvyrep.jpg Box plot of Y loading reproducibility across split-halfs colored based on Z-score.
# repsummary.jpg Box plots of all reproducibility scores for significant dimensions.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(TExPosition)
library(PTCA4CATA)
library(data4PCCAR)
library(MVN)

#Catches arguments.
args <- commandArgs(trailingOnly=T)
k <- args[1]

#Set base path and PLSC path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
plscpath <-  paste0(basepath,'allreconfig/PLSC/')

#Set parameters.
nk <- as.numeric(k)
theme_set(theme_gray())
iternum <- '10000'
permver <- 'byColumns'
bootcv <- '2.5'

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

# PLSC --------------------------------------------------------------------

#Run PLSC.
data1 <- data_X
data2 <- data_Y
pls.res <- tepPLS(data1,data2, 
                  center1=T,scale1=T, 
                  center2=T,scale2=T, 
                  graphs = FALSE)

# Permutation Significance ------------------------------------------------

#Generate outpath.
outpath <- paste0(plscpath,'perm/')
dir.create(outpath,recursive=T)

#PLSC permutation model fit.
set.seed(12345)
perm.plsc <- perm4PLSC(data1, data2, 
                       center1=T,scale1=T, 
                       center2=T,scale2=T, 
                       nIter=as.numeric(iternum),
                       permType=permver,
                       compact=F)

#Save singular value table.
singexp <- pls.res$TExPosition.Data$pdq$t
singval <- pls.res$TExPosition.Data$pdq$Dv
nsv <- min(ncol(data1),ncol(data2))
nsv.perm <- length(singval)
singcor <- matrix(NA,nrow=1,ncol=nsv)
for (i in 1:nsv.perm) {
  clx <- pls.res$TExPosition.Data$lx[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  singcor[,i] <- cor(clx,cly,method='spearman')
}
singp <- perm.plsc$pEigenvalues
singsig <- as.character(perm.plsc$pEigenvalues < 0.05)
svmat <- matrix(NA,nrow=5,ncol=nsv)
svmat[1,1:length(singexp)] <- singexp
svmat[2,1:length(singval)] <- singval
svmat[3,1:length(singcor)] <- singcor
svmat[4,1:length(singp)] <- singp
svmat[5,1:length(singsig)] <- singsig
rownames(svmat) <- c('Explained','SV','Correlation','P-value','Significance')
outfile <- paste0(outpath,'permsig.csv')
write.csv(svmat,outfile,row.names=T)

#Plot.
sigcolors <- rep('nosig',nsv)
sigcolors[perm.plsc$pEigenvalues < 0.05] <- 'sig'
lvplot <- data.frame(singexp)
colnames(lvplot) <- 'Explained'
lvplot['Significant'] <- sigcolors
lvplot['LV'] <- factor(paste0('LV',as.character(1:nsv)),
                       levels=paste0('LV',as.character(1:nsv)))
ggobj <- ggplot(lvplot,aes(x=LV,y=Explained,group=1)) +
  geom_line() +
  geom_point(aes(colour=factor(Significant),size=factor(Significant))) +
  scale_colour_manual(name='P-value',
                      values=c("black","red"),
                      labels=c('≥ 0.05','< 0.05')) +
  theme(axis.title.x=element_blank()) +
  ylab('Variance Explained') 
outfile <- paste0(outpath,'permsig.jpg')
ggsave(outfile,ggobj,units='in',width=5,height=2)

#Set significant LVs.
permsig <- (perm.plsc$pEigenvalues < 0.05)
svlist <- which(permsig)
nsv <- length(svlist)

#Generate outpath.
outpath <- paste0(plscpath,'observations/')
dir.create(outpath,recursive=T)

#Set colors.
norm.color <- '#696969'
ysig.color <- '#4E79A7'
xsig.color <- '#4E79A7'
pos.color <- '#FF0000'
neg.color <- '#4E79A7'

#Plot Y loadings.
for (lidx in svlist) {
  
  #Get significance label.
  if (permsig[lidx]) {
    addsig = '*'
  } else {
    addsig = ''
  }
  
  #Plot.
  llab <- as.character(lidx)
  ly <- pls.res$TExPosition.Data$fj[,lidx]
  if (permsig[lidx]) {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(ysig.color,ncog),
                           font.size=3) 
  } else {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(norm.color,ncog),
                           font.size=3)
  }
  
  #Save.
  outfile <- paste0(outpath,'baseLVY',llab,'.jpg')
  ggsave(outfile,plty,units='in',width=5,height=3)
}

#Generate flippers.
outfile <- paste0(outpath,'yflip.csv')
if (file.exists(outfile)) {
  yflip <- as.vector(t(read.csv(outfile,header=F)))
} else {
  
  #Set up flippers based on specific loadings.
  yflip <- rep(F,each=nsv)

  #k = 2.
  if (nk==2) {
    if (pls.res$TExPosition.Data$fj['ProcSpeed',1] < 0) {
      yflip[1] <- T
    } 
    if (pls.res$TExPosition.Data$fj['PicSeq',2] < 0) {
      yflip[2] <- T
    } 
  
  #k = 3-4.
  } else if (nk==3|nk==4) {
    if (mean(pls.res$TExPosition.Data$fj[,1]) < 0) {
      yflip[1] <- T
    } 
    if (pls.res$TExPosition.Data$fj['ProcSpeed',2] < 0) {
      yflip[2] <- T
    } 

  #k = 5-11.
  } else if (nk==5|nk==6|nk==7|nk==8|nk==9|nk==10|nk==11) {
    if (mean(pls.res$TExPosition.Data$fj[,1]) < 0) {
      yflip[1] <- T
    } 
    if (pls.res$TExPosition.Data$fj['ProcSpeed',2] < 0) {
      yflip[2] <- T
    } 
    if (pls.res$TExPosition.Data$fj['PicSeq',3] < 0) {
      yflip[3] <- T
    }
  
  #k = 12.
  } else if (nk == '12') {
    if (mean(pls.res$TExPosition.Data$fj[,1]) < 0) {
      yflip[1] <- T
    } 
    if (pls.res$TExPosition.Data$fj['PicSeq',2] < 0) {
      yflip[2] <- T
    } 
    if (pls.res$TExPosition.Data$fj['ProcSpeed',3] < 0) {
      yflip[3] <- T
    }
  
  #Default.
  } else {
    if (mean(pls.res$TExPosition.Data$fj[,1]) < 0) {
      yflip[1] <- T
    } 
  }
  
  #Generate file.
  outflip <- t(yflip) 
  outflip[outflip == TRUE] = 'T'
  outflip[outflip == FALSE] = 'F'
  write.table(outflip,outfile,sep=',',row.names=F,col.names=F)
}

#Generate outpath.
outpath <- paste0(plscpath,'perm/')
dir.create(outpath,recursive=T)

#Produce max values based on flipped versions.
clc_xlim <- c()
clc_ylim <- c()
for (i in 1:nsv) {
  
  #Extract correlation.
  clx <- pls.res$TExPosition.Data$lx[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  if (yflip[i]) {
    clx <- -clx
    cly <- -cly
  }
  clc_xlim <- c(clc_xlim,max(clx),min(clx))
  clc_ylim <- c(clc_ylim,max(cly),min(cly))
}
clc_xmax <- max(clc_xlim)
clc_xmin <- min(clc_xlim)
clc_ymax <- max(clc_ylim)
clc_ymin <- min(clc_ylim)

#Plot scatterplots and do bivariate normality tests based on flipped versions.
testlabs <- c('mardia_skew','mardia_kurtosis','hz','royston','dh','energy')
ntest <- length(testlabs)
pr_norm <- matrix(NA,nrow=nsv,ncol=ntest)
rownames(pr_norm) <- paste0('LV',as.character(1:nsv))
colnames(pr_norm) <- testlabs
clc_txtsize <- 11
for (i in 1:nsv) {
  
  #Extract correlation.
  clx <- pls.res$TExPosition.Data$lx[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  if (yflip[i]) {
    clx <- -clx
    cly <- -cly
  }
  ccor <- cor(clx,cly,method='spearman')
  clc <- data.frame('LVX'=clx,'LVY'=cly)
  clc_labs <- paste0(colnames(clc),as.character(i))
  
  #Plot.
  pltc <- ggplot(aes(LVX,LVY),data=clc) +
    geom_point(size=1) +
    geom_smooth(method=lm,se=F) +
    scale_x_continuous(limits=c(clc_xmin,clc_xmax)) +
    scale_y_continuous(limits=c(clc_ymin,clc_ymax)) +
    labs(x=clc_labs[1],y=clc_labs[2]) +
    theme(axis.text.x=element_text(size=clc_txtsize),
          axis.text.y=element_text(size=clc_txtsize))
  
  #Save.
  outfile <- paste0(outpath,'LV',as.character(i),'_corr.jpg')
  ggsave(outfile,pltc,units='in',width=4,height=3)
  
  #Set current bivariate relationship.
  c_XY <- clc
  
  #Generate bivariate QQ plot and save.
  outfile <- paste0(outpath,'LV',as.character(i),'_qq.jpg')
  jpeg(file=outfile)
  mvn(c_XY,multivariatePlot='qq')
  dev.off()
  # grid.echo()
  # qqp <- grid.grab()
  
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
    pr_norm[paste0('LV',as.character(i)),ctest] <- cpr
  }
}

#Threshold p-values and save.
pr_norm_thres <- pr_norm
pr_norm_thres[pr_norm < 0.05] <- NA
pr_norm_thres <- signif(pr_norm_thres,2)
outfile <- paste0(outpath,'binormtest_p.csv')
write.csv(data.frame(pr_norm),outfile)
outfile <- paste0(outpath,'binormtest_pthres.csv')
write.csv(data.frame(pr_norm_thres),outfile)

#Save scores for flipped versions.
lx_scores <- matrix(NA,nrow=nrow(data_X),ncol=nsv)
rownames(lx_scores) <- rownames(data_X)
ly_scores <- lx_scores
for (i in 1:nsv) {
  
  #Extract scores.
  clx <- pls.res$TExPosition.Data$lx[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  if (yflip[i]) {
    clx <- -clx
    cly <- -cly
  }
  lx_scores[,i] <- clx
  ly_scores[,i] <- cly
}
colnames(lx_scores) <- paste0('LVX',as.character(1:nsv))
colnames(ly_scores) <- paste0('LVY',as.character(1:nsv))
outfile <- paste0(outpath,'LVX_scores.csv')
write.csv(lx_scores,outfile)
outfile <- paste0(outpath,'LVY_scores.csv')
write.csv(ly_scores,outfile)

# Bootstrap Stability -----------------------------------------------------

#Generate outpath.
outpath <- paste0(plscpath,'boot/')
dir.create(outpath,recursive=T)

#PLSC bootstrap feature stability.
bootres.plsc <- Boot4PLSC(data1,data2,
                          eig=T,
                          nIter=as.numeric(iternum),
                          critical.value=as.numeric(bootcv),
                          nf2keep=nsv)


#Generate color for significant bootstrap ratios.
bootcolor.X <- ifelse(bootres.plsc[["bootRatiosSignificant.i"]],xsig.color,norm.color)
bootcolor.Y <- ifelse(bootres.plsc[["bootRatiosSignificant.j"]],ysig.color,norm.color)

#Plot X and Y loadings.
scaler <- 1
plot_list <- list()
pidx <- 1
lxtab <- matrix(NA,nrow=nsv,ncol=npred)
colnames(lxtab) <- predlabs
lytab <- matrix(NA,nrow=nsv,ncol=ncog)
colnames(lytab) <- coglabs
lxtab_full <- lxtab
lytab_full <- lytab
for (lidx in 1:nsv) {
  
  #Get significance label.
  if (permsig[lidx]) {
    addsig = '*'
  } else {
    addsig = ''
  }
  
  #Get and flip target Y.
  ly <- pls.res$TExPosition.Data$fj[,lidx]
  if (yflip[lidx]) {
    ly <- -ly
  }
  
  #Scale.
  ly <- ly*scaler
  
  #Plot.
  llab <- as.character(lidx)
  if (permsig[lidx]) {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=bootcolor.Y[,lidx],
                           font.size=3) 
  } else {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(norm.color,ncog),
                           font.size=3)
  }
  
  #Save.
  plot_list[[pidx]] <- plty
  pidx <- pidx + 1
  
  #Plot second version.
  collab <- ly
  collab[ly > 0] <- pos.color
  collab[ly < 0] <- neg.color
  collab[!(bootres.plsc[["bootRatiosSignificant.j"]][,lidx])] <- norm.color
  if (permsig[lidx]) {
    plty2 <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=collab,
                           font.size=3) 
  } else {
    plty2 <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(norm.color,ncog),
                           font.size=3)
  }
  outfile <- paste0(outpath,'LVY',llab,'_load.jpg')
  ggsave(outfile,plty2,units='in',width=2,height=3)
  
  #Get and flip target X.
  lx <- pls.res$TExPosition.Data$fi[,lidx] 
  if (yflip[lidx]) {
    lx <- -lx
  }
  
  #Scale.
  lx <- lx*scaler
  
  #Plot.
  if (permsig[lidx]) {
    pltx <- PrettyBarPlot2(lx,
                           threshold=0, 
                           color4bar=bootcolor.X[,lidx],
                           font.size=3)
  } else {
    pltx <- PrettyBarPlot2(lx,
                           threshold=0, 
                           color4bar=rep(norm.color,npred),
                           font.size=3)
  }  
  
  #Save.
  plot_list[[pidx]] <- pltx
  pidx <- pidx + 1
  
  #Save the loading tables.
  lx_thres <- lx
  lxtab_full[lidx,] <- lx_thres
  lx_thres[!bootres.plsc[["bootRatiosSignificant.i"]][,lidx]] <- NA
  lxtab[lidx,] <- lx_thres
  ly_thres <- ly
  lytab_full[lidx,] <- ly_thres
  ly_thres[!bootres.plsc[["bootRatiosSignificant.j"]][,lidx]] <- NA
  lytab[lidx,] <- ly_thres
  
  #Plot second version.
  dynlabs <- c('occur','dwell','numtrans','trans_','sim_','jump_')
  ndyn <- length(dynlabs)
  for (didx in 1:ndyn) {
    
    #Extract.
    cdyn <- dynlabs[didx]
    dlx <- data.frame(t(lx)) %>%
      select(matches(cdyn))
    if (grepl('_',cdyn)) {
      namedyn <- gsub('_','',cdyn)
    } else {
      namedyn <- cdyn
    }
    dsig <- bootres.plsc[["bootRatiosSignificant.i"]][colnames(dlx),lidx]
    
    #Set up labels for plotting.
    if (namedyn == 'occur') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'dwell') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'numtrans') {
      inlabs <- paste0('Transition Number')
      nwidth <- 1
      nheight <- 3
    } else if (namedyn == 'trans') {
      inlabs <- gsub('trans_','',colnames(dlx))
      inlabs <- gsub('t','-',colnames(dlx))
      nwidth <- 10
      nheight <- 3
    } else if (namedyn == 'sim') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'jump') {
      inlabs <- gsub('jump_','',colnames(dlx))
      inlabs <- gsub('t','-',colnames(dlx))
      nwidth <- 10
      nheight <- 3
    }
    
    #Plot.
    collab <- dlx
    collab[dlx > 0] <- pos.color
    collab[dlx < 0] <- neg.color
    collab[!(dsig)] <- norm.color
    dlx <- unlist(as.vector(dlx))
    names(dlx) <- inlabs
    if (permsig[lidx]) {
      pltdx <- PrettyBarPlot2(dlx,
                              threshold=0, 
                              color4bar=collab,
                              font.size=3) 
    } else {
      pltdx <- PrettyBarPlot2(dlx,
                              threshold=0, 
                              color4bar=rep(norm.color,length(dlx)),
                              font.size=3)
    }
    outfile <- paste0(outpath,'LVX',llab,'_load_',namedyn,'.jpg')
    ggsave(outfile,pltdx,units='in',width=nwidth,height=nheight)
  }
}

#Format plots.
p1 <- plot_list[c(1,3,5)]
p2 <- plot_list[c(2,4,6)]

#Combine plots.
outfile <- paste0(outpath,'LVYall_load.jpg')
ggsave(outfile,arrangeGrob(grobs=p1,ncol=1),
       units='in',width=2,height=8)
outfile <- paste0(outpath,'LVXall_load.jpg')
ggsave(outfile,arrangeGrob(grobs=p2,ncol=1),
       units='in',width=30,height=8)

#Save tables.
outfile <- paste0(outpath,'LVXall_load.csv')
write.csv(lxtab,outfile)
outfile <- paste0(outpath,'LVYall_load.csv')
write.csv(lytab,outfile)
outfile <- paste0(outpath,'LVXall_load_full.csv')
write.csv(t(lxtab_full),outfile)
outfile <- paste0(outpath,'LVYall_load_full.csv')
write.csv(t(lytab_full),outfile)

#Plot X and Y bootstrap ratios.
scaler <- 1
plot_list <- list()
pidx <- 1
lxtab <- matrix(NA,nrow=nsv,ncol=npred)
colnames(lxtab) <- predlabs
lytab <- matrix(NA,nrow=nsv,ncol=ncog)
colnames(lytab) <- coglabs
lxtab_full <- lxtab
lytab_full <- lytab
for (lidx in 1:nsv) {
  
  #Get significance label.
  if (permsig[lidx]) {
    addsig = '*'
  } else {
    addsig = ''
  }
  
  #Get and flip target Y.
  ly <- bootres.plsc$bootRatios.j[,lidx]
  if (yflip[lidx]) {
    ly <- -ly
  }
  
  #Scale.
  ly <- ly*scaler
  
  #Plot.
  llab <- as.character(lidx)
  if (permsig[lidx]) {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=bootcolor.Y[,lidx],
                           font.size=3) 
  } else {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(norm.color,ncog),
                           font.size=3)
  }
  
  #Save.
  plot_list[[pidx]] <- plty
  pidx <- pidx + 1
  
  #Plot second version.
  collab <- ly
  collab[ly > 0] <- pos.color
  collab[ly < 0] <- neg.color
  collab[!(bootres.plsc[["bootRatiosSignificant.j"]][,lidx])] <- norm.color
  if (permsig[lidx]) {
    plty2 <- PrettyBarPlot2(ly,
                            threshold=0, 
                            color4bar=collab,
                            font.size=3) 
  } else {
    plty2 <- PrettyBarPlot2(ly,
                            threshold=0, 
                            color4bar=rep(norm.color,ncog),
                            font.size=3)
  }
  outfile <- paste0(outpath,'LVY',llab,'_bootrat.jpg')
  ggsave(outfile,plty2,units='in',width=2,height=3)
  
  #Get and flip target X.
  lx <- bootres.plsc$bootRatios.i[,lidx]
  if (yflip[lidx]) {
    lx <- -lx
  }
  
  #Scale.
  lx <- lx*scaler
  
  #Plot.
  if (permsig[lidx]) {
    pltx <- PrettyBarPlot2(lx,
                           threshold=0, 
                           color4bar=bootcolor.X[,lidx],
                           font.size=3)
  } else {
    pltx <- PrettyBarPlot2(lx,
                           threshold=0, 
                           color4bar=rep(norm.color,npred),
                           font.size=3)
  }  
  
  #Save.
  plot_list[[pidx]] <- pltx
  pidx <- pidx + 1
  
  #Save the loading tables.
  lx_thres <- lx
  lxtab_full[lidx,] <- lx_thres
  lx_thres[!bootres.plsc[["bootRatiosSignificant.i"]][,lidx]] <- NA
  lxtab[lidx,] <- lx_thres
  ly_thres <- ly
  lytab_full[lidx,] <- ly_thres
  ly_thres[!bootres.plsc[["bootRatiosSignificant.j"]][,lidx]] <- NA
  lytab[lidx,] <- ly_thres
  
  #Plot second version.
  dynlabs <- c('occur','dwell','numtrans','trans_','sim_','jump_')
  ndyn <- length(dynlabs)
  for (didx in 1:ndyn) {
    
    #Extract.
    cdyn <- dynlabs[didx]
    dlx <- data.frame(t(lx)) %>%
      select(matches(cdyn))
    if (grepl('_',cdyn)) {
      namedyn <- gsub('_','',cdyn)
    } else {
      namedyn <- cdyn
    }
    dsig <- bootres.plsc[["bootRatiosSignificant.i"]][colnames(dlx),lidx]
    
    #Set up labels for plotting.
    if (namedyn == 'occur') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'dwell') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'numtrans') {
      inlabs <- paste0('Transition Number')
      nwidth <- 1
      nheight <- 3
    } else if (namedyn == 'trans') {
      inlabs <- gsub('trans_','',colnames(dlx))
      inlabs <- gsub('t','-',colnames(dlx))
      nwidth <- 10
      nheight <- 3
    } else if (namedyn == 'sim') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'jump') {
      inlabs <- gsub('jump_','',colnames(dlx))
      inlabs <- gsub('t','-',colnames(dlx))
      nwidth <- 10
      nheight <- 3
    }
    
    #Plot.
    collab <- dlx
    collab[dlx > 0] <- pos.color
    collab[dlx < 0] <- neg.color
    collab[!(dsig)] <- norm.color
    dlx <- unlist(as.vector(dlx))
    names(dlx) <- inlabs
    if (permsig[lidx]) {
      pltdx <- PrettyBarPlot2(dlx,
                              threshold=0, 
                              color4bar=collab,
                              font.size=3) 
    } else {
      pltdx <- PrettyBarPlot2(dlx,
                              threshold=0, 
                              color4bar=rep(norm.color,length(dlx)),
                              font.size=3)
    }
    outfile <- paste0(outpath,'LVX',llab,'_bootrat_',namedyn,'.jpg')
    ggsave(outfile,pltdx,units='in',width=nwidth,height=nheight)
  }
}

#Format plots.
p1 <- plot_list[c(1,3,5)]
p2 <- plot_list[c(2,4,6)]

#Combine plots.
outfile <- paste0(outpath,'LVYall_bootrat.jpg')
ggsave(outfile,arrangeGrob(grobs=p1,ncol=1),
       units='in',width=2,height=8)
outfile <- paste0(outpath,'LVXall_bootrat.jpg')
ggsave(outfile,arrangeGrob(grobs=p2,ncol=1),
       units='in',width=30,height=8)

#Save tables.
outfile <- paste0(outpath,'LVXall_bootrat.csv')
write.csv(lxtab,outfile)
outfile <- paste0(outpath,'LVYall_bootrat.csv')
write.csv(lytab,outfile)
outfile <- paste0(outpath,'LVXall_bootrat_full.csv')
write.csv(t(lxtab_full),outfile)
outfile <- paste0(outpath,'LVYall_bootrat_full.csv')
write.csv(t(lytab_full),outfile)

# Reproducibility ---------------------------------------------------------

#Create replication outpath.
outpath <- paste0(plscpath,'replicate/')
dir.create(outpath,recursive=T)

#Create a matrix for the SV, LVX, and LVY replicability values.
niter <- as.numeric(iternum)
maxsv <- min(ncol(data_X),ncol(data_Y))
svrep <- matrix(NA,nrow=niter,ncol=maxsv)
lvxrep <- matrix(NA,nrow=niter,ncol=maxsv)
lvyrep <- matrix(NA,nrow=niter,ncol=maxsv)

#For each of the splits.
nsample <- nrow(data_Y)
nsplit <- floor(nsample/2)
set.seed(12345)
for (spidx in 1:niter) {
  
  #Extract splits.
  idx <- sample(1:nsample,replace=F)
  idx_1 <- idx[1:nsplit]
  idx_2 <- idx[(ceiling(nsample/2)+1):nsample]
  x1 <- data_X[idx_1,]
  y1 <- data_Y[idx_1,]
  x2 <- data_X[idx_2,]
  y2 <- data_Y[idx_2,]
  
  #Do PLS.
  pls1 <- tepPLS(x1,y1, 
                 center1=T,scale1=T, 
                 center2=T,scale2=T, 
                 graphs=F)
  pls2 <- tepPLS(x2,y2, 
                 center1=T,scale1=T, 
                 center2=T,scale2=T, 
                 graphs=F)
  
  #Calculate second correlation matrix.
  rxy2 <- cor(x2,y2) 
  
  #Get first and second u and v.
  u1 <- pls1$TExPosition.Data$pdq$p 
  v1 <- pls1$TExPosition.Data$pdq$q 
  u2 <- pls2$TExPosition.Data$pdq$p 
  v2 <- pls2$TExPosition.Data$pdq$q
  
  #Calculate SV reproducibility value using the second set Rxy and first set u and v.
  svrep_val <- t(u1)%*%rxy2%*%v1 
  svrep_val <- diag(svrep_val)
  svrep[spidx,] <- svrep_val
  
  #Calculate LVX and LVY reproducibility value using the first and second set u and v.
  lvxrep_val <- t(u1)%*%u2
  lvxrep_val <- abs(diag(lvxrep_val))
  lvxrep[spidx,] <- lvxrep_val
  lvyrep_val <- t(v1)%*%v2
  lvyrep_val <- abs(diag(lvyrep_val))
  lvyrep[spidx,] <- lvyrep_val
}

#Find the mean and std columnwise and z-score.
svrep_mean <- colMeans(svrep)
svrep_std <- apply(svrep,2,sd)
svrep_z <- svrep_mean / svrep_std
lvxrep_mean <- colMeans(lvxrep)
lvxrep_std <- apply(lvxrep,2,sd)
lvxrep_z <- lvxrep_mean / lvxrep_std
lvyrep_mean <- colMeans(lvyrep)
lvyrep_std <- apply(lvyrep,2,sd)
lvyrep_z <- lvyrep_mean / lvyrep_std

#Save each file.
outlabs <- c('svrep','svrep_z','lvxrep','lvxrep_z','lvyrep','lvyrep_z')
nout <- length(outlabs)
outvars <- list(svrep,svrep_z,lvxrep,lvxrep_z,lvyrep,lvyrep_z)
for (oidx in 1:nout) {
  
  #Extract.
  olab <- outlabs[[oidx]]
  ovar <- outvars[[oidx]]
  if (is.matrix(ovar)) {
    colnames(ovar) <- paste0('LV',1:maxsv)
  } else {
    ovar <- t(data.frame(ovar))
    colnames(ovar) <- paste0('LV',1:maxsv)
  }
  
  #Save.
  outfile <- paste0(outpath,olab,'.csv')
  write.csv(ovar,outfile,row.names=F)
}

#Create replication plot outpath.
plotpath <- paste0(plscpath,'replicate/plots/')
dir.create(plotpath,recursive=T)

#Label LVs.
lvlabs <- paste0('LV',1:maxsv)

#Plot SV replicability.
if (length(unique(svrep_z > 1.95)) == 2) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
} else if (unique(svrep_z > 1.95) == F) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('≤ 1.95')))
} else if (unique(svrep_z > 1.95) == T) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('> 1.95')))
}
svrep_color['LV'] <- lvlabs
plotmat <- data.frame(svrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(svrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
svrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
  geom_boxplot(outlier.size=0.1) +
  ggtitle('SV Reproducibility') +
  ylab('Reproducibility') +
  theme(axis.title.x = element_blank())
outfile <- paste0(plotpath,'svrep.jpg')
ggsave(outfile,svrep_plot,units='in')

#Plot LVX replicability.
if (length(unique(lvxrep_z > 1.95)) == 2) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
} else if (unique(svrep_z > 1.95) == F) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('≤ 1.95')))
} else if (unique(svrep_z > 1.95) == T) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('> 1.95')))
}
lvxrep_color['LV'] <- lvlabs
plotmat <- data.frame(lvxrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(lvxrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
lvxrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
  geom_boxplot(outlier.size=0.1) +
  ggtitle('LVX Reproducibility') +
  ylab('Reproducibility') +
  theme(axis.title.x = element_blank())
outfile <- paste0(plotpath,'lvxrep.jpg')
ggsave(outfile,lvxrep_plot,units='in')

#Plot LVY replicability.
if (length(unique(lvyrep_z > 1.95)) == 2) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
} else if (unique(svrep_z > 1.95) == F) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('≤ 1.95')))
} else if (unique(svrep_z > 1.95) == T) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('> 1.95')))
}
lvyrep_color['LV'] <- lvlabs
plotmat <- data.frame(lvyrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(lvyrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
lvyrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
  geom_boxplot(outlier.size=0.1) +
  ggtitle('LVY Reproducibility') +
  ylab('Reproducibility') +
  theme(axis.title.x = element_blank())
outfile <- paste0(plotpath,'lvyrep.jpg')
ggsave(outfile,lvyrep_plot,units='in')

#Isolate the significant LVs only.
siglv <- svlist
lvlabs <- paste0('LV',as.character(siglv))
svrep <- svrep[,siglv]
svrep_z <- svrep_z[siglv]
lvxrep <- lvxrep[,siglv]
lvxrep_z <- lvxrep_z[siglv]
lvyrep <- lvyrep[,siglv]
lvyrep_z <- lvyrep_z[siglv]
nsv <- length(siglv)

#Set summary text size.
sumtxt <- 10

#SV.
if (length(unique(svrep_z > 1.95)) == 2) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
  ncolors <- 'ft'
} else if (unique(svrep_z > 1.95) == F) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('≤ 1.95')))
  ncolors <- 'fone'
} else if (unique(svrep_z > 1.95) == T) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('> 1.95')))
  ncolors <- 'tone'
}
svrep_color['LV'] <- lvlabs
plotmat <- data.frame(svrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(svrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
if (ncolors == 'ft') {
  svrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    ggtitle('Singular Values') +
    ylab('Reproducibility') +
    theme(axis.title.x = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'fone') {
  svrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#F8766D')) +
    ggtitle('Singular Values') +
    ylab('Reproducibility') +
    theme(axis.title.x = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'tone') {
  svrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#00BFC4')) +
    ggtitle('Singular Values') +
    ylab('Reproducibility') +
    theme(axis.title.x = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
}

#LVX.
if (length(unique(lvxrep_z > 1.95)) == 2) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
  ncolors <- 'ft'
} else if (unique(lvxrep_z > 1.95) == F) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('≤ 1.95')))
  ncolors <- 'fone'
} else if (unique(lvxrep_z > 1.95) == T) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('> 1.95')))
  ncolors <- 'tone'
}
lvxrep_color['LV'] <- lvlabs
plotmat <- data.frame(lvxrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(lvxrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
if (ncolors == 'ft') {
  lvxrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    ggtitle('X Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'fone') {
  lvxrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#F8766D')) +
    ggtitle('X Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'tone') {
  lvxrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#00BFC4')) +
    ggtitle('X Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
}

#LVY.
if (length(unique(lvyrep_z > 1.95)) == 2) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
  ncolors <- 'ft'
} else if (unique(lvyrep_z > 1.95) == F) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('≤ 1.95')))
  ncolors <- 'fone'
} else if (unique(lvyrep_z > 1.95) == T) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('> 1.95')))
  ncolors <- 'tone'
}
lvyrep_color['LV'] <- lvlabs
plotmat <- data.frame(lvyrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(lvyrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
if (ncolors == 'ft') {
  lvyrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    ggtitle('Y Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'fone') {
  lvyrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#F8766D')) +
    ggtitle('Y Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'tone') {
  lvyrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#00BFC4')) +
    ggtitle('Y Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
}  
(svrep_plot + lvyrep_plot + lvxrep_plot)
ggsave(paste0(plotpath,'repsummary.jpg'),units='in',width=5.25,height=3)
dev.off()
