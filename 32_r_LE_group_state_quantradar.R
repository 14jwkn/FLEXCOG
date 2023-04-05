# For the subject group and the specified k, separately average the positive and
# negative LE values for each region for each network to produce a score.
# Output:
# clab,'.jpg' Contains the network radar plots for average positive and negative LE values.

#Set libraries.
library(tidyverse)
library(reticulate)
library(ggplot2)
library(gridExtra)
library(ggradar)

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

#Select parameters.
maxk <- 6
klist <- as.character(c(4:6))
knum <- length(klist)
ncolor <- 'red'
pcolor <- 'blue'
bcolor <- c(ncolor,pcolor)

#Read in net labels.
netlabels <- scan('colenetlabels.txt',character(),quote='') 

#Produce namelist.
namelist <- c('VIS1','VIS2','SMN',
              'CON','DAN','LAN',
              'FPN','AUD','DMN',
              'PMM','VMM',
              'ORA')
nnet <- length(namelist)

#Produce networks where 1 is within and 0 is outside and the number of regions
#classified to each network.
netmat <- matrix(0,nrow=360,ncol=nnet)
colnames(netmat) <- namelist
for (nidx in 1:nnet) {
  
  #Extract.
  netid <- as.character(nidx)
  
  #Make it 1 where the region is classified into the network.
  netmat[netlabels == netid,nidx] = 1
}
netsum <- t(data.frame(colSums(netmat)))

#Read in the sorting.
classfile = paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/',
                   subgroup,'/allk/',subgroup,'_sortclass.csv')
classmat = read.csv(classfile,row.names=1)
rownames(classmat) <- klist
colnames(classmat) <- as.character(c(1:maxk))

#Read in the best iterations.
bestfile = paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/',
                  'group/best_iter/',subgroup,'/best_iter.csv')
bestiter = read.csv(bestfile,row.names=1)

#Extract best iteration.
iteration <- as.character(bestiter[subgroup,paste0(clean,'_',order,'_',roiname,'_',k)])

#Read in the LE.
inpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                 subgroup,'/',k,'/')
infile <- paste0(inpath,'subcent_',iteration,'.h5')
inkey <- paste0('/',subgroup)
instore <- h5py$File(infile,'r')
clustmat <- np$array(instore[inkey])
instore$close()

#Get max and min LE values.
kmax.LE <- max(clustmat)
kmin.LE <- min(clustmat)

#Extract the original states, universal states, and universal index.
orig <- classmat[k,]
orig <- orig[!is.na(orig)]
uval <- colnames(classmat[k,!is.na(classmat[k,])])
nk <- as.numeric(k)
sort.tab <- matrix(NA,nrow=nk,ncol=3)
colnames(sort.tab) <- c('Orig','Uval','Uidx')
sort.tab[,'Orig'] <- orig
sort.tab[,'Uval'] <- uval
sort.tab[,'Uidx'] <- c(1:nk)

#Generate max and min positive and absolute negative values across states.
pos.avg.lim <- c()
neg.avg.lim <- c()
for (uidx in c(1:nk)) {
  
  #Extract original state and universal state.
  origstate <- orig[uidx]
  ustate <- uval[uidx]
  
  #Read in the original state, threshold, and turn negatives to positives.
  kstate <- clustmat[,origstate]
  pos.kstate <- kstate
  pos.kstate[kstate < 0] <- 0
  neg.kstate <- kstate
  neg.kstate[kstate > 0] <- 0
  neg.kstate <- -neg.kstate
  
  #Find values.
  pos.kstate.rep <- matrix(pos.kstate,nrow=360,ncol=nnet)
  pos.kstate.values <- pos.kstate.rep * netmat
  neg.kstate.rep <- matrix(neg.kstate,nrow=360,ncol=nnet)
  neg.kstate.values <- neg.kstate.rep * netmat
  
  #Find averages.
  pos.kstate.avg <- t(data.frame(colSums(pos.kstate.values)/colSums(!!pos.kstate.values)))
  pos.kstate.avg[is.na(pos.kstate.avg)] <- 0
  neg.kstate.avg <- t(data.frame(colSums(neg.kstate.values)/colSums(!!neg.kstate.values)))
  neg.kstate.avg[is.na(neg.kstate.avg)] <- 0
  
  #Find limits.
  pos.avg.lim <- c(pos.avg.lim,max(pos.kstate.avg),min(pos.kstate.avg))
  neg.avg.lim <- c(neg.avg.lim,max(neg.kstate.avg),min(neg.kstate.avg))
}

#Max-min values.
pavgmax <- max(pos.avg.lim)
pavgmmin <- min(pos.avg.lim)
navgmax <- max(neg.avg.lim)
navgmin <- min(neg.avg.lim)

#Produce network radar plots.
bavgplt <- list()
pidx <- 1
for (uidx in c(1:nk)) {
  
  #Extract original state and universal state.
  origstate <- orig[uidx]
  ustate <- uval[uidx]
  
  #Read in the original state, threshold, and turn negatives to positives.
  kstate <- clustmat[,origstate]
  pos.kstate <- kstate
  pos.kstate[kstate < 0] <- 0
  neg.kstate <- kstate
  neg.kstate[kstate > 0] <- 0
  neg.kstate <- -neg.kstate
  
  #Find values.
  pos.kstate.rep <- matrix(pos.kstate,nrow=360,ncol=nnet)
  pos.kstate.values <- pos.kstate.rep * netmat
  neg.kstate.rep <- matrix(neg.kstate,nrow=360,ncol=nnet)
  neg.kstate.values <- neg.kstate.rep * netmat
  
  #Find averages.
  pos.kstate.avg <- t(data.frame(colSums(pos.kstate.values)/colSums(!!pos.kstate.values)))
  pos.kstate.avg[is.na(pos.kstate.avg)] <- 0
  neg.kstate.avg <- t(data.frame(colSums(neg.kstate.values)/colSums(!!neg.kstate.values)))
  neg.kstate.avg[is.na(neg.kstate.avg)] <- 0
  
  #Add a column for groups.
  adder = data.frame('Group'=0)
  badder = data.frame('Group'=c('P','N'))
  both.avg <- rbind(pos.kstate.avg,neg.kstate.avg)
  both.avg <- cbind(badder,both.avg)
  
  #Create radar for averages.
  maxrad <- max(pavgmax,navgmax)
  midrad <- maxrad/2
  minrad <- 0
  ggobj <- ggradar(both.avg,
                   values.radar = c('','',''),
                   grid.min=minrad,grid.mid=midrad,grid.max=maxrad,
                   axis.label.size=3,
                   group.line.width=1,
                   group.point.size=0,
                   group.colours=bcolor,
                   legend.position='none')
  bavgplt[[pidx]] <- ggobj
  
  #Increment.
  pidx <- pidx + 1
}

#Plot and save.
outpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                 subgroup,'/',k,'/netradar/')
dir.create(outpath,recursive=T)
clab <- 'bothavg'
cvar <- bavgplt
outfile <- paste0(outpath,clab,'.jpg')
ggsave(outfile,arrangeGrob(grobs=cvar,ncol=1),units='in',width=3,height=10)
