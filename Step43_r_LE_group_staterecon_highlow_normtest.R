# For the subject group and the specified k, plot the Q-Q and scatter plots for the
# sets of summary network reconfiguration metrics and state dFC(t) strength and variability.
# Output:
# 'highlow_',ctype,'_qqplot.jpg' Q-Q plots for either the summary reconfiguration metrics, dFC strength, or dFC variability.
# 'highlow_',ctarg,'_scatter.jpg' Scatter plots for the specific relationship between dFC strength or variability with summary reconfiguration metrics with uniform axis limits.
# 'highlow_',ctarg,'_scatter_nolim.jpg' Scatter plots for the specific relationship between dFC strength or variability with summary reconfiguration metrics with automatic axis limits.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]
k <- args[2]

#Set base path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/')
outpath <- paste0(basepath,'allreconfig_subsets/normtest/')

#Set parameters.
nk <- as.numeric(k)
targlist <- c('g_jumppos','g_auto','g_simneg',
              'pr_126far','pr_auto','pr_simpos')
ntarg <- length(targlist)
targlabs <- c('g Between','g Within','g Idiosyncrasy',
              'PS Between','PS Within','PS Idiosyncrasy')

#Read in strvar.
infile <- paste0(basepath,'dFC_strvar.csv')
mat_dFC <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject') 
strvar_labs <- c(colnames(mat_dFC)[2:ncol(mat_dFC)])

#Get confounds.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Define confounds.
con_list <- c('Age','Gender')
con_pred <- paste0(con_list,collapse='+')
con_mat <- othercog[,con_list]

#Join confounds and values.
val_pred <- strvar_labs
resmat <- cbind(mat_dFC[,val_pred],con_mat)
colnames(resmat) <- c(val_pred,con_list)

#For each value.
for (val_lab in val_pred) {
  
  #Produce residuals.
  inform <- as.formula(paste0(val_lab,'~',con_pred))
  respred <- resid(lm(inform,data=resmat))
  
  #Replace.
  mat_dFC[,val_lab] <- respred
}

#Read in reconfiguration summary values which have already been confound regressed.
nsub <- nrow(mat_dFC)
metvals <- matrix(NA,nrow=nsub,ncol=ntarg)
for (taidx in 1:ntarg) {
  targtype <- targlist[taidx]
  infile <- paste0(hlpath,targtype,'/',targtype,'_values.csv')
  cvals <- read_csv(infile,col_names=T) %>% 
    rename_at(1,~'Subject')
  metvals[,taidx] <- unlist(cvals[,'Value'])
}
colnames(metvals) <- targlist
metvals <- cbind(mat_dFC[,'Subject'],metvals)

#Combine data into one table.
data_XY <- metvals %>%
  inner_join(mat_dFC) 
mainlabs <- colnames(data_XY)
vlabs <- mainlabs[!grepl('gPCA',mainlabs)]
vlabs <- vlabs[!grepl('Subject',vlabs)]

# QQ Plots ----------------------------------------------------------------

#Set theme.
theme_set(theme_bw())

#Go through each variable.
vartypes <- c('mean','std','summary')
nvtypes <- length(vartypes)
for (vtidx in 1:nvtypes) {
  ctype <- vartypes[vtidx]
  
  #Extract all variables which match.
  if (ctype=='summary') {
    cvarlist <- targlist
    cvarlabs <- targlabs
  } else {
    cvarlist <- vlabs[grepl(ctype,vlabs)]
    if (ctype=='mean') {
      cvarlabs <- gsub('mean_','',cvarlist)
    } else if (ctype=='std') {
      cvarlabs <- gsub('std_','',cvarlist)
    } 
  }
  nvars <- length(cvarlist)
  
  #Initiate plots.
  plist <- list()
  pidx <- 1
  for (vidx in 1:nvars) {
    cvar <- cvarlist[vidx]
    clab <- cvarlabs[vidx]
    
    #Generate plot.
    p1 <- ggplot(data_XY,aes_string(sample=cvar)) +
      stat_qq() +
      stat_qq_line() + 
      ggtitle(clab) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust=0.5),
            plot.margin=grid::unit(c(1,1,1,1),'mm')) 
    
    #Append plot.
    plist[[pidx]] <- p1
    pidx <- pidx + 1
  }
  
  #Save plots.
  if (ctype=='summary') {
    ncols <- 6
    cwidth <- 4*ncols
    cheight <- 4
  } else {
    ncols <- nk
    cwidth <- 4*ncols
    cheight <- 4
  } 
  outfile <- paste0(outpath,'highlow_',ctype,'_qqplot.jpg')
  ggsave(outfile,arrangeGrob(grobs=plist,ncol=ncols),
         units='in',width=cwidth,height=cheight)
}

# Scatter Plots -----------------------------------------------------------

#Go through each variable.
pad <- 0
for (taidx in 1:ntarg) {
  ctarg <- targlist[taidx]
  
  #Extract Y limits.
  ymax <- max(data_XY[,ctarg])
  ymin <- min(data_XY[,ctarg])
  ydiff <- ymax - ymin
  yinc <- ydiff*pad
  ymax <- ymax + yinc
  ymin <- ymin - yinc
  
  #Extract all variables which match selected mean and std.
  if ((x1=='g_jumppos')|(x1=='pr_126far')) {
    cvarlist <- vlabs[grepl('mean_',vlabs)]
    cvarlabs <- gsub('mean_','STR ',cvarlist)
  } else if ((x1=='g_auto')|(x1=='pr_auto')) {
    cvarlist <- vlabs[grepl('std_',vlabs)]
    cvarlabs <- gsub('std_','VAR ',cvarlist)
  } else if ((x1=='g_simneg')|(x1=='pr_simpos')) {
    cvarlist <- vlabs[grepl('mean_',vlabs)|grepl('std_',vlabs)]
    cvarlabs <- gsub('mean_','STR ',cvarlist)
    cvarlabs <- gsub('std_','VAR ',cvarlabs)
  }
  nvars <- length(cvarlist)
  
  #Find X limits for each set.
  xlim_mean <- c()
  xlim_std <- c()
  for (vidx in 1:nvars) {
    cvar <- cvarlist[vidx]
    cx <- data_XY[,cvar]
    if (grepl('mean_',cvar)) {
      xlim_mean <- c(xlim_mean,max(cx))
      xlim_mean <- c(xlim_mean,min(cx))
    } else if (grepl('std_',cvar)) {
      xlim_std <- c(xlim_std,max(cx))
      xlim_std <- c(xlim_std,min(cx))
    }
  }
  xmax_mean <- max(xlim_mean)
  xmin_mean <- min(xlim_mean)
  xmax_std <- max(xlim_std)
  xmin_std <- min(xlim_std)
  
  #Initiate plots with and without limits.
  plist <- list()
  plist_nolim <- list()
  pidx <- 1
  for (vidx in 1:nvars) {
    cvar <- cvarlist[vidx]
    clab <- cvarlabs[vidx]
    
    #Generate plot.
    if (grepl('mean',cvar)) {
      xmax <- xmax_mean
      xmin <- xmin_mean
    } else if (grepl('std',cvar)) {
      xmax <- xmax_std
      xmin <- xmin_std
    }
    xdiff <- xmax - xmin
    xinc <- xdiff*pad
    xmax <- xmax + xinc
    xmin <- xmin - xinc
    p1 <- ggplot(data_XY,aes_string(x=cvar,y=ctarg)) +
      geom_point() +
      geom_smooth(method='lm',se=F,fullrange=T,color='black') +
      scale_x_continuous(limits=c(xmin,xmax)) +   
      scale_y_continuous(limits=c(ymin,ymax)) +
      ggtitle(clab) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust=0.5),
            plot.margin=grid::unit(c(1,1,1,1),'mm')) 
    p2 <- ggplot(data_XY,aes_string(x=cvar,y=ctarg)) +
      geom_point() +
      geom_smooth(method='lm',se=F,fullrange=T,color='black') +
      ggtitle(clab) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust=0.5),
            plot.margin=grid::unit(c(1,1,1,1),'mm')) 
    
    #Append plot.
    plist[[pidx]] <- p1
    plist_nolim[[pidx]] <- p2
    pidx <- pidx + 1
  }
  
  #Save plots.
  defnum <- 4
  ncols <- nk
  cwidth <- defnum*ncols
  if ((x1=='g_jumppos')|(x1=='pr_126far')) {
    nrows <- 1
  } else if ((x1=='g_auto')|(x1=='pr_auto')) {
    nrows <- 1
  } else if ((x1=='g_simneg')|(x1=='pr_simpos')) {
    nrows <- 2
  }
  cheight <- defnum*nrows
  outfile <- paste0(outpath,'highlow_',ctarg,'_scatter.jpg')
  ggsave(outfile,arrangeGrob(grobs=plist,ncol=ncols),
         units='in',width=cwidth,height=cheight)
  outfile <- paste0(outpath,'highlow_',ctarg,'_scatter_nolim.jpg')
  ggsave(outfile,arrangeGrob(grobs=plist_nolim,ncol=ncols),
         units='in',width=cwidth,height=cheight)
}
