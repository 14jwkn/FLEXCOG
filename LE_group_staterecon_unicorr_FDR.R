# For the specified k and cognitive variable, show box plots for the distributions
# of reconfiguration metrics across subjects. Additionally, explore coloring by
# significant correlations for box plots and bar plots.
# Output:
# unipath,'distplots/,namedyn,'.jpg' Box plots for each reconfiguration metric.
# unipath,'distplots/,curr_lab,'.jpg' Colored box plots for each metric.
# unipath,'barplots/',namedyn,'.jpg' Colored bar plots for each metric.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(scico)
library(PTCA4CATA)

#Catches arguments.
k <- args[1]
ccog <- args[2]

#Set base path and univariate path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
unipath <- paste0(basepath,'allreconfig/unicorr/',ccog,'/')

#Set parameters.
nk <- as.numeric(k)
pos.color <- '#FF0000'
neg.color <- '#4E79A7'
norm.color <- '#696969'
theme_set(theme_minimal())

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
meanmat <- dynmat %>%
  inner_join(simmat) %>%
  inner_join(jumpmat)
colnames(meanmat) <- gsub('-','t',colnames(meanmat))

#Get thresholded correlations.
infile <- paste0(unipath,'corr_fullfdr_thres.csv')
thcorr <- read.csv(infile) %>%
  replace(is.na(.),0)
thcorr[,1] <- gsub('sptrans','trans',thcorr[,1])
maxcorr <- max(thcorr[,2]) + 0.05
mincorr <- min(thcorr[,2]) - 0.05

#Get actual correlations and significance.
infile <- paste0(unipath,'corr.csv')
realcorr <- read.csv(infile) 
infile <- paste0(unipath,'pcorr_fullfdr.csv')
sigcorr <- read.csv(infile)
sigcorr[,2] <- ifelse((sigcorr[,2] < 0.05),T,F)

#Go through each set of values.
dynlabs <- c('occur','dwell','numtrans','trans_','sim_','jump_')
ndyn <- length(dynlabs)
for (didx in 1:ndyn) {
  
  #Extract current set of values, and their max and min.
  cdyn <- dynlabs[didx]
  cmean <- meanmat %>%
    select(matches(cdyn))
  maxval <- max(cmean)
  minval <- min(cmean)
  if (grepl('_',cdyn)) {
    namedyn <- gsub('_','',cdyn)
  } else {
    namedyn <- cdyn
  }
  
  #Set up labels for plotting.
  if (namedyn == 'occur') {
    inlabs <- as.character(1:nk)
    barlabs <- as.character(1:nk)
    nwidth <- 3
    nheight <- 3
    bwidth <- 1.5
    bheight <- 3
    xbreaks <- seq(minval,maxval,by=0.05)
    xlabels <- round(xbreaks,2)
  } else if (namedyn == 'dwell') {
    inlabs <- as.character(1:nk)
    barlabs <- as.character(1:nk)
    nwidth <- 3
    nheight <- 3
    bwidth <- 1.5
    bheight <- 3
    xbreaks <- seq(minval,maxval,by=1)
    xlabels <- round(xbreaks,0)
  } else if (namedyn == 'numtrans') {
    inlabs <- paste0(' ')
    barlabs <- inlabs
    nwidth <- 2
    nheight <- 5
    bwidth <- 1
    bheight <- 3
    xbreaks <- seq(minval,maxval,by=50)
    xlabels <- round(xbreaks,0)
  } else if (namedyn == 'trans') {
    inlabs <- gsub('trans_','',colnames(cmean))
    inlabs <- gsub('t','-',inlabs)
    barlabs <- inlabs
    nwidth <- 12
    nheight <- 5
    bwidth <- 10
    bheight <- 3
    xbreaks <- seq(minval,maxval,by=0.05)
    xlabels <- round(xbreaks,2)
  } else if (namedyn == 'sim') {
    inlabs <- as.character(1:nk)
    barlabs <- paste0('State ',1:nk)
    nwidth <- 3
    nheight <- 3
    bwidth <- 1.5
    bheight <- 2.5
    xbreaks <- seq(minval,maxval,by=1)
    xlabels <- round(xbreaks,0)
  } else if (namedyn == 'jump') {
    inlabs <- gsub('jump_','',colnames(cmean))
    inlabs <- gsub('t','-',inlabs)
    barlabs <- inlabs
    nwidth <- 12
    nheight <- 5
    bwidth <- 10
    bheight <- 3
    xbreaks <- seq(minval,maxval,by=1)
    xlabels <- round(xbreaks,0)
  }
  npred <- length(inlabs)
  
  #Format plotting matrix.
  if ((cdyn == 'trans_')|(cdyn == 'jump_')) {
    plotmat <- cmean %>%
      gather(label,stat) %>%
      arrange(label)
    plotmat[,'label'] <- factor(plotmat$label)
    num_one <- nrow(plotmat)/nk
    plotmat[,'state'] <- factor(rep(paste0(1:nk),each=num_one))
    plotmat[,'label'] <- gsub('trans_','',plotmat$label)
    plotmat[,'label'] <- gsub('jump_','',plotmat$label)
    plotmat[,'label'] <- gsub('t','-',plotmat$label)
  } else {
    plotmat <- cmean %>%
      gather(label,stat) %>%
      arrange(label)
  }
  
  #Extract current set.
  ccorr <- thcorr[,c('Label',ccog)]
  ccorr <- data.frame(t(ccorr))
  colnames(ccorr) <- ccorr[1,]
  ccorr <- ccorr %>%
    select(matches(cdyn))
  ccorr <- ccorr[-c(1),]
  coglong <- as.data.frame(lapply(ccorr,rep,nrow(cmean))) %>%
    gather(label,corr) %>%
    arrange(label) 
  
  #Format plotting matrix.
  if ((cdyn == 'trans_')|(cdyn == 'jump_')) {
    plotcorr <- cbind(plotmat,coglong[,'corr'])
    colnames(plotcorr) <- c('label','stat','state','corr')
    plotcorr$corr <- as.numeric(plotcorr$corr)
  } else {
    plotcorr <- cbind(plotmat,coglong[,'corr'])
    colnames(plotcorr) <- c('label','stat','corr')
    plotcorr$corr <- as.numeric(plotcorr$corr)
  }
  
  #Plot box plots.
  if ((cdyn == 'trans_')|(cdyn == 'jump_')) {
    boxmat <- ggplot(plotcorr,aes(x=label,y=stat)) + 
      geom_boxplot(position='dodge',width=0.9) +
      facet_grid(~state,scales='free',space='free_x') +
      theme(axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            strip.text.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank()) +
      scale_y_continuous(breaks=xbreaks,
                         labels=xlabels) +
      scale_x_discrete(expand=c(0.2,0)) +
      ggtitle(namedyn) 
  } else {
    boxmat <- ggplot(plotcorr,aes(x=label,y=stat)) + 
      geom_boxplot(position='dodge',width=0.9) +
      theme(axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank()) +
      scale_y_continuous(breaks=xbreaks,
                         labels=xlabels) +
      scale_x_discrete(labels=inlabs) +
      ggtitle(namedyn) 
  }
  
  #Plot box plots with coloring.
  if ((cdyn == 'trans_')|(cdyn == 'jump_')) {
    curr_lab <- paste0(cdyn,ccog)
    boxcorr <- ggplot(plotcorr,aes(x=label,y=stat,fill=corr)) + 
      geom_boxplot(position='dodge',width=0.9) +
      facet_grid(~state,scales='free',space='free_x') +
      theme(axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            strip.text.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank()) +
      scale_fill_scico(palette='berlin',
                       midpoint=0,
                       limits=c(mincorr,maxcorr)) +
      scale_y_continuous(breaks=xbreaks,
                         labels=xlabels) +
      scale_x_discrete(expand=c(0.2,0)) +
      ggtitle(curr_lab) + 
      labs(fill='corr')
  } else {
    curr_lab <- paste0(namedyn,'_',ccog)
    boxcorr <- ggplot(plotcorr,aes(x=label,y=stat,fill=corr)) + 
      geom_boxplot(position='dodge',width=0.9) +
      theme(axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank()) +
      scale_fill_scico(palette='berlin',
                       midpoint=0,
                       limits=c(mincorr,maxcorr)) +
      scale_y_continuous(breaks=xbreaks,
                         labels=xlabels) +
      scale_x_discrete(labels=inlabs) +
      ggtitle(curr_lab) + 
      labs(fill='Correlation')
  }
  
  #Set outpath and save.
  plotpath <- paste0(unipath,'distplots/')
  dir.create(plotpath,recursive=T)
  ggsave(paste0(plotpath,namedyn,'.jpg'),boxmat,units='in',width=nwidth,height=nheight,bg='white')
  ggsave(paste0(plotpath,curr_lab,'.jpg'),boxcorr,units='in',width=nwidth,height=nheight,bg='white')
  
  #Extract bars.
  cbar <- realcorr[,c('Label',ccog)]
  cbar <- data.frame(t(cbar))
  colnames(cbar) <- cbar[1,]
  cbar <- cbar %>%
    select(matches(cdyn))
  cbar <- cbar[-c(1),]
  
  #Extract significance.
  csig <- sigcorr[,c('Label',ccog)]
  csig <- data.frame(t(csig))
  colnames(csig) <- csig[1,]
  csig <- csig %>%
    select(matches(cdyn))
  csig <- csig[-c(1),]

  #Create color labels.
  collab <- cbar
  collab[cbar > 0] <- pos.color
  collab[cbar < 0] <- neg.color
  collab[!as.logical(csig)] <- norm.color
  
  #Plot.
  pbar <- as.numeric(cbar)
  names(pbar) <- barlabs
  barp <- PrettyBarPlot2(pbar,
                         threshold=0, 
                         color4bar=collab,
                         font.size=3) 
  
  #Set outpath and save.
  plotpath <- paste0(unipath,'barplots/')
  dir.create(plotpath,recursive=T)
  ggsave(paste0(plotpath,namedyn,'.jpg'),barp,units='in',width=bwidth,height=bheight)
}
