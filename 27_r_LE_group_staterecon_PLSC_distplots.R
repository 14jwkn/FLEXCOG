# For the subject group and the specified k, color the distribution of reconfiguration
# values across subjects in box plots by the PLSC loading with different colors 
# for stable values. 
# Output:
# curr_lab,'.jpg' Colored box plot for a given reconfiguration metric.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(scico)

#Catches arguments.
subgroup <- args[1]
k <- args[2]

#Set base path, PLSC path, and plot path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/')
plscpath <-  paste0(basepath,'allreconfig/PLSC/')

#Set parameters.
nk <- as.numeric(k)
theme_set(theme_minimal())
bootver <- 'load'

#Retrieve the loadings or bootstrap ratios for each LV.
infile <- paste0(plscpath,'boot/LVXall_',bootver,'.csv')
xboot <- read.csv(infile) %>%
  replace(is.na(.),0)
maxboot <- max(xboot[,2:ncol(xboot)]) + 0.5
minboot <- min(xboot[,2:ncol(xboot)]) - 0.5
svlist <- xboot[,1]
nsv <- length(svlist)
svlabs <- paste0('LV',svlist)
xboot[,1] <- svlabs
colnames(xboot) <- c('LV',colnames(xboot)[2:ncol(xboot)])

#Retrieve the subject mean values for state dynamics, similarity, and jump.
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
    nwidth <- 5
    nheight <- 5
    xbreaks <- seq(minval,maxval,by=0.05)
    xlabels <- round(xbreaks,2)
  } else if (namedyn == 'dwell') {
    inlabs <- as.character(1:nk)
    nwidth <- 5
    nheight <- 5
    xbreaks <- seq(minval,maxval,by=1)
    xlabels <- round(xbreaks,0)
  } else if (namedyn == 'numtrans') {
    inlabs <- paste0(' ')
    nwidth <- 2
    nheight <- 5
    xbreaks <- seq(minval,maxval,by=50)
    xlabels <- round(xbreaks,0)
  } else if (namedyn == 'trans') {
    inlabs <- gsub('trans_','',colnames(cmean))
    nwidth <- 12
    nheight <- 5
    xbreaks <- seq(minval,maxval,by=0.05)
    xlabels <- round(xbreaks,2)
  } else if (namedyn == 'sim') {
    inlabs <- as.character(1:nk)
    nwidth <- 5
    nheight <- 5
    xbreaks <- seq(minval,maxval,by=1)
    xlabels <- round(xbreaks,0)
  } else if (namedyn == 'jump') {
    inlabs <- gsub('jump_','',colnames(cmean))
    nwidth <- 12
    nheight <- 5
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
  } else {
    plotmat <- cmean %>%
      gather(label,stat) %>%
      arrange(label)
  }
  
  #Go through each LV.
  for (lidx in svlist) {
    
    #Extract current set.
    clv <- svlabs[lidx]
    cboot <- xboot %>%
      select(matches(paste0(cdyn,'|LV')))
    cboot <- data.frame(t(cboot))
    colnames(cboot) <- cboot[1,]
    cboot <- cboot[-c(1),]
    cboot <- cboot[,clv]
    coglong <- as.data.frame(lapply(cboot,rep,nrow(cmean))) 
    colnames(coglong) <- 1:npred
    coglong <- coglong %>%
      gather(label,boot) 
    coglong[,'label'] <- as.numeric(coglong[,'label'])
    coglong <- coglong %>%
      arrange(label)
    coglong[,'label'] <- as.character(coglong[,'label'])
    
    #Format plotting matrix.
    if ((cdyn == 'trans_')|(cdyn == 'jump_')) {
      plotboot <- cbind(plotmat,coglong[,'boot'])
      colnames(plotboot) <- c('label','stat','state','boot')
      plotboot$boot <- as.numeric(plotboot$boot)
    } else {
      plotboot <- cbind(plotmat,coglong[,'boot'])
      colnames(plotboot) <- c('label','stat','boot')
      plotboot$boot <- as.numeric(plotboot$boot)
    }
    
    #Plot box plots.
    curr_lab <- paste0(bootver,'_',namedyn,'_',clv)
    if ((cdyn == 'trans_')|(cdyn == 'jump_')) {
      boxboot <- ggplot(plotboot,aes(x=label,y=stat,fill=boot)) + 
        geom_boxplot(position='dodge',width=0.9) +
        facet_grid(~state,scales='free',space='free_x') +
        theme(axis.title.y=element_blank(),
              axis.title.x=element_blank(),
              strip.text.x=element_blank(),
              panel.grid.minor.y=element_blank(),
              panel.grid.minor.x=element_blank()) +
        scale_fill_scico(palette='berlin',
                         midpoint=0,
                         limits=c(minboot,maxboot)) +
        scale_y_continuous(breaks=xbreaks,
                           labels=xlabels) +
        scale_x_discrete(expand=c(0.2,0)) +
        ggtitle(curr_lab) + 
        labs(fill='Bootstrap')
    } else {
      boxboot <- ggplot(plotboot,aes(x=label,y=stat,fill=boot)) + 
        geom_boxplot(position='dodge',width=0.9) +
        theme(axis.title.y=element_blank(),
              axis.title.x=element_blank(),
              panel.grid.minor.y=element_blank(),
              panel.grid.minor.x=element_blank()) +
        scale_fill_scico(palette='berlin',
                         midpoint=0,
                         limits=c(minboot,maxboot)) +
        scale_y_continuous(breaks=xbreaks,
                           labels=xlabels) +
        scale_x_discrete(labels=inlabs) +
        ggtitle(curr_lab) + 
        labs(fill='Bootstrap')
    }
    
    #Set outpath and save.
    plotpath <- paste0(plscpath,'distplots/',clv,'/')
    dir.create(plotpath,recursive=T)
    ggsave(paste0(plotpath,curr_lab,'.jpg'),boxboot,units='in',width=nwidth,height=nheight,bg='white')
  }
}
