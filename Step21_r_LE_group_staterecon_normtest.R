# For the subject group and the specified k, plot the Q-Q and scatter plots for
# the dynamics metrics and the cognitive variable.
# Output:
# ctype,'_qqplot.jpg' Q-Q plots for the specific group of dynamics or cognitive metrics.
# ccog,'_',ctype,'_scatter.jpg' Scatter plots for the specific group of cognitive-dynamic relationships for g with uniform axis limits.
# ccog,'_',ctype,'_scatter_nolim.jpg' Scatter plots for the specific group of cognitive-dynamic relationships for g with automatic axis limits.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)

#Catches arguments.
subgroup <- args[1]
k <- args[2]

#Set base path and out path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
outpath <- paste0(basepath,'allreconfig/normtest/')
dir.create(outpath,recursive=T)

#Set parameters.
nk <- as.numeric(k)
ccog <- 'gPCA'

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Isolate desired cognitive variables.
ourcog <- colnames(othercog)
matchlist <- paste0('Subject|gPCA|CardSort|Flanker|ListSort|PicSeq|PV|PS|ReadEng|P24_CR|',
                      'IWRD_TOT|VSPLOT_TC')
ourcog <- str_subset(ourcog,matchlist)
cogmat <- othercog[,ourcog]
ncog <- ncol(cogmat) - 1

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

#Transform Y column names.
if (septype == 'sep') {
  data_Y <- data_Y %>%
    rename('PMAT24'='P24_CR',
           'PicVoc'='PV',
           'ProcSpeed'='PS',
           'WordMem'='IWRD_TOT',
           'ReadVoc'='ReadEng',
           'LineOrient'='VSPLOT_TC')
}

#Extract labels.
ncog <- ncol(data_Y)
npred <- ncol(data_X)
coglabs <- colnames(data_Y)
predlabs <- colnames(data_X)

# QQ Plots ----------------------------------------------------------------

#Set theme.
theme_set(theme_bw())

#Go through each variable.
data_XY <- cbind(data_Y,data_X)
vlabs <- c(coglabs,predlabs)
vartypes <- c('cog','occur','dwell','numtrans','trans','sim','jump')
nvtypes <- length(vartypes)
for (vtidx in 1:nvtypes) {
  ctype <- vartypes[vtidx]
  
  #Extract all variables which match.
  if (ctype=='cog') {
    cvarlist <- c('gPCA','CardSort','Flanker','ListSort','PicSeq','PicVoc','ProcSpeed',
                  'ReadEng','PMAT24','WordMem','LineOrient')
    cvarlabs <- cvarlist
  } else {
    cvarlist <- vlabs[grepl(ctype,vlabs)]
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
  if (ctype=='cog') {
    ncols <- 3
    cwidth <- 4*ncols
    cheight <- 4
  } else if ((ctype=='occur')|(ctype=='dwell')|(ctype=='sim')) {
    ncols <- nk
    cwidth <- 4*ncols
    cheight <- 4
  } else if ((ctype=='trans')|(ctype=='jump')) {
    ncols <- nk
    cwidth <- 4*ncols
    cheight <- 4*ncols
  } else if (ctype=='numtrans') {
    ncols <- 1
    cwidth <- 4
    cheight <- 4
  }
  outfile <- paste0(outpath,ctype,'_qqplot.jpg')
  ggsave(outfile,arrangeGrob(grobs=plist,ncol=ncols),
         units='in',width=cwidth,height=cheight)
}

# Scatter Plots -----------------------------------------------------------

#Set theme.
theme_set(theme_bw())

#Extract.
data_XY <- cbind(data_Y,data_X)
vlabs <- predlabs
vartypes <- c('occur','dwell','numtrans','trans','sim','jump')
nvtypes <- length(vartypes)
ccog <- 'gPCA'

#Set up limits.
ymax <- max(data_XY[,ccog])
ymin <- min(data_XY[,ccog])
ydiff <- ymax - ymin
yinc <- ydiff*0.1
ymax <- ymax + yinc
ymin <- ymin - yinc

#Go through each variable.
for (vtidx in 1:nvtypes) {
  ctype <- vartypes[vtidx]
  
  #Extract all variables which match.
  cvarlist <- vlabs[grepl(ctype,vlabs)]
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
  
  #Initiate plots with and without limits.
  plist <- list()
  plist_nolim <- list()
  pidx <- 1
  for (vidx in 1:nvars) {
    cvar <- cvarlist[vidx]
    clab <- cvarlabs[vidx]
    
    #Generate plot.
    p1 <- ggplot(data_XY,aes_string(x=cvar,y=ccog)) +
      geom_point() +
      geom_smooth(method='lm',se=F,fullrange=T,color='black') +
      scale_x_continuous(limits=c(xmin,xmax)) +   
      scale_y_continuous(limits=c(ymin,ymax)) +
      ggtitle(clab) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust=0.5),
            plot.margin=grid::unit(c(1,1,1,1),'mm')) 
    p2 <- ggplot(data_XY,aes_string(x=cvar,y=ccog)) +
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
  if (ctype=='cog') {
    ncols <- 3
    cwidth <- defnum*ncols
    cheight <- defnum
  } else if ((ctype=='occur')|(ctype=='dwell')|(ctype=='sim')) {
    ncols <- nk
    cwidth <- defnum*ncols
    cheight <- defnum
  } else if ((ctype=='trans')|(ctype=='jump')) {
    ncols <- nk
    cwidth <- defnum*ncols
    cheight <- defnum*ncols
  } else if (ctype=='numtrans') {
    ncols <- 1
    cwidth <- defnum
    cheight <- defnum
  }
  outfile <- paste0(outpath,ccog,'_',ctype,'_scatter.jpg')
  ggsave(outfile,arrangeGrob(grobs=plist,ncol=ncols),
         units='in',width=cwidth,height=cheight)
  outfile <- paste0(outpath,ccog,'_',ctype,'_scatter_nolim.jpg')
  ggsave(outfile,arrangeGrob(grobs=plist_nolim,ncol=ncols),
         units='in',width=cwidth,height=cheight)
}
