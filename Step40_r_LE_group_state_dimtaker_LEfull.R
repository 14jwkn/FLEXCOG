# For the subject group and the specified k, plot the 75% and 25% Q1-3 whiskers
# for the PCs for each state. For three random subjects, plot LE(t) as dots where
# colors are states along three PC dimensions.
# Output:
# q3med_PC12.jpg PC1-PC2 dimension Q3 whisker plots for each state.
# q3med_PC23.jpg PC2-PC3 dimension Q3 whisker plots for each state.
# q3med_PC13.jpg PC1-PC3 dimension Q3 whisker plots for each state.
# 'sub',widx,'_LEvis.jpg' PC1-3 dimension dot plots for subjects colored by state.

#Set libraries.
library(tidyverse)
library(reticulate)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(gg3D)

#Set virtual environment.
use_virtualenv('../envs/FlexEnv')

#Load python packages.
h5py <- import('h5py')
np <- import ('numpy')
pd <- import('pandas')

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]
k <- args[2]
print(paste0(subgroup,k))

#Set base path and out path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
outpath <-  paste0(basepath,'lowdim_LEfull/')

#Read subjects.
if (subgroup == 'full') {
  infile <- 'r_full_submain.txt'
} else if (subgroup == 'half') {
  infile <- 'r_half_submain.txt'
}
subjects <- scan(infile,character(),quote='')

#Set parameters.
nwin <- 4792
nk <- as.numeric(k)
nsubs <- length(subjects)
ndim <- 3
theme_set(theme_classic())
pwidth <- 7.5
pheight <- 6.5
errwidth <- 0.04
errheight <- 0.04
linsize <- 2
dotsize <- 8
rowfull <- nwin*nsubs + nk

#Extract centroids and data.
infile <- paste0(outpath,'LEfull_PCreduce.h5')
instore <- pd$HDFStore(infile,'r')
inkey <- '/dFClow'
allmat <- instore$select(inkey)
inkey <- '/dFClow_lab'
clulab <- instore$select(inkey)
instore$close()
allmat <- cbind(clulab,allmat)
colnames(allmat) <- c('Subject','Cluster','PC1','PC2','PC3')
centmat <- allmat[1:nk,]
datmat <- allmat[(nk+1):rowfull,]

#Resort centroids.
classfile <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/',
                    subgroup,'/allk/',subgroup,'_sortclass.csv')
classmat <- read.csv(classfile,row.names=1)
colnames(classmat) <- 1:ncol(classmat)
classmat <- classmat[k,]
classmat <- as.vector(t(na.omit(t(classmat))))
centmat <- centmat[classmat,]
centmat[,'Cluster'] <- 1:nk

# Q3 Median ---------------------------------------------------------------

#Separate out Q3 max and min values for each cluster and dimension.
q3labs <- c(paste0('q3max_',1:ndim),paste0('q3med_',1:ndim),paste0('q3min_',1:ndim))
q3lim = matrix(NA,nrow=(length(q3labs)),ncol=nk)
rownames(q3lim) = q3labs
colnames(q3lim) = 1:nk
for (kidx in 1:nk) {
  
  #Extract.
  cmat <- datmat[(datmat[,'Cluster'] == kidx),]
  for (didx in 1:ndim) {
    
    #Extract.
    cdim <- paste0('PC',didx)
    outlab <- paste0('q3max_',didx)
    q3lim[outlab,kidx] <- quantile(cmat[,cdim],prob=0.75)
    outlab <- paste0('q3med_',didx)
    q3lim[outlab,kidx] <- quantile(cmat[,cdim],prob=0.50)
    outlab <- paste0('q3min_',didx)
    q3lim[outlab,kidx] <- quantile(cmat[,cdim],prob=0.25)
  }
}

#Set up.
plotmat <- cbind(centmat,t(q3lim))
plotmat[,'Cluster'] <- factor(plotmat[,'Cluster'])

#Do PC1 and PC2.
xdim <- '1'
ydim <- '2'
pc12 <- ggplot(plotmat,aes_string(x=paste0('q3med_',xdim),y=paste0('q3med_',ydim),color='Cluster')) + 
  geom_point(size=dotsize) +
  geom_errorbar(aes_string(ymin=paste0('q3min_',ydim),ymax=paste0('q3max_',ydim)),
                width=errwidth,size=linsize) +
  geom_errorbarh(aes_string(xmin=paste0('q3min_',xdim),xmax=paste0('q3max_',xdim)),
                 height=errheight,size=linsize) +
  xlab(paste0('PC',xdim)) +
  ylab(paste0('PC',ydim)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
print(pc12)

#Do PC2 and PC3.
xdim <- '3'
ydim <- '2'
pc23 <- ggplot(plotmat,aes_string(x=paste0('q3med_',xdim),y=paste0('q3med_',ydim),color='Cluster')) + 
  geom_point(size=dotsize) +
  geom_errorbar(aes_string(ymin=paste0('q3min_',ydim),ymax=paste0('q3max_',ydim)),
                width=errwidth,size=linsize) +
  geom_errorbarh(aes_string(xmin=paste0('q3min_',xdim),xmax=paste0('q3max_',xdim)),
                 height=errheight,size=linsize) +
  xlab(paste0('PC',xdim)) +
  ylab(paste0('PC',ydim)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
print(pc23)

#Do PC1 and PC3.
xdim <- '3'
ydim <- '1'
pc13 <- ggplot(plotmat,aes_string(x=paste0('q3med_',xdim),y=paste0('q3med_',ydim),color='Cluster')) + 
  geom_point(size=dotsize) +
  geom_errorbar(aes_string(ymin=paste0('q3min_',ydim),ymax=paste0('q3max_',ydim)),
                width=errwidth,size=linsize) +
  geom_errorbarh(aes_string(xmin=paste0('q3min_',xdim),xmax=paste0('q3max_',xdim)),
                 height=errheight,size=linsize) +
  xlab(paste0('PC',xdim)) +
  ylab(paste0('PC',ydim)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) 
print(pc13)

#Save these versions.
outfile <- paste0(outpath,'q3med_PC12.jpg')
ggsave(outfile,pc12,units='in',width=pwidth,height=pheight)
outfile <- paste0(outpath,'q3med_PC23.jpg')
ggsave(outfile,pc23,units='in',width=pwidth,height=pheight)
outfile <- paste0(outpath,'q3med_PC13.jpg')
ggsave(outfile,pc13,units='in',width=pwidth,height=pheight)

# Subject Plots -----------------------------------------------------------

#Select random subjects.
nwant <- 3
set.seed(12345)
wantsubs <- sample(subjects,nwant)

#Collect max and min values for each dimension.
totlabs <- c(paste0('max_',1:ndim),paste0('min_',1:ndim))
totlim <- matrix(NA,nrow=length(totlabs),ncol=1)
rownames(totlim) = totlabs
colnames(totlim) = 1
for (widx in 1:nwant) {
  
  #Extract.
  csub <- wantsubs[widx]
  subidx <- match(csub,subjects)
  print(paste0('Doing: ',csub))
  
  #Produce its indices.
  datstart = (nwin*(subidx-1)) + 1
  datend = (nwin*subidx)
  
  #Extract its values.
  cmat = datmat[datstart:datend,]
  cmat[,'Cluster'] <- factor(cmat[,'Cluster'])
  
  #For each dimension.
  for (didx in 1:ndim) {
    
    #Extract.
    cdim <- paste0('PC',didx)
    outlab <- paste0('max_',didx)
    totlim[outlab,1] <- max(c(totlim[outlab,1],max(cmat[,cdim])),na.rm=T)
    outlab <- paste0('min_',didx)
    totlim[outlab,1] <- min(c(totlim[outlab,1],min(cmat[,cdim])),na.rm=T)
  }
}

#Adjust.
for (didx in 1:ndim) {
  
  #Do separately.
  if (didx == 1) {
    
    #Extract.
    cdim <- paste0('PC',didx)
    outlab <- paste0('max_',didx)
    totlim[outlab,1] <- totlim[outlab,1] - 0.6
    outlab <- paste0('min_',didx)
    totlim[outlab,1] <- totlim[outlab,1] + 0.3
    
  } else {
    
    #Extract.
    cdim <- paste0('PC',didx)
    outlab <- paste0('max_',didx)
    totlim[outlab,1] <- totlim[outlab,1] - 0.3
    outlab <- paste0('min_',didx)
    totlim[outlab,1] <- totlim[outlab,1] + 0.3
  }
}

#Plot each subject's data.
for (widx in 1:nwant) {
  
  #Extract.
  csub <- wantsubs[widx]
  subidx <- match(csub,subjects)
  
  #Produce its indices.
  datstart = (nwin*(subidx-1)) + 1
  datend = (nwin*subidx)
  
  #Extract its values.
  cmat = datmat[datstart:datend,]
  cmat[,'Cluster'] <- factor(cmat[,'Cluster'])
  
  #Plot.
  theta <- 360
  phi <- 0
  xdim <- '1'
  ydim <- '3'
  zdim <- '2'
  ggobj <- ggplot(cmat,aes_string(x=paste0('PC',xdim),
                                  y=paste0('PC',ydim),
                                  z=paste0('PC',zdim),
                                  xmin=totlim[paste0('min_',xdim),1],
                                  xmax=totlim[paste0('max_',xdim),1],
                                  ymin=totlim[paste0('min_',ydim),1],
                                  ymax=totlim[paste0('max_',ydim),1],
                                  zmin=totlim[paste0('min_',zdim),1],
                                  zmax=totlim[paste0('max_',zdim),1],
                                  colour='Cluster')) +
    axes_3D(theta=theta,phi=phi) +
    stat_3D(theta=theta,phi=phi,size=0.5) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line=element_blank())
  
  # Save.
  outfile <- paste0(outpath,'sub',widx,'_LEvis.jpg')
  ggsave(outfile,ggobj,units='in',width=7,height=7)
}
