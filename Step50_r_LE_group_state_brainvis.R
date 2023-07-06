# For the subject group, plot the ROIs with positive LE values for each state
# across the ks of interest on the brain.
# Output:
# 'max',as.character(maxk),'_',onehemi,'_',onelat,'_ksortvis_brain.jpg' Brain plots of positive LE values for each state.
# 'max',as.character(maxk),'_ksortvis_brain_clbr.jpg' Corresponding color bar.

#Set libraries.
library(tidyverse)
library(reticulate)
library(ggplot2)
library(gridExtra)
library(ggseg)
library(ggsegGlasser)
library(scico)

#Set virtual environment.
use_virtualenv('../envs/FlexEnv')

#Load python packages.
h5py <- import('h5py')
np <- import ('numpy')
pd <- import('pandas')

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]

#Set k to investigate.
maxk <- 12
klist <- as.character(c(2:maxk))
knum <- length(klist)

#Read in the sorting.
classfile = paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/',
                   subgroup,'/allk/',subgroup,'_sortclass.csv')
classmat = read.csv(classfile,row.names=1)
rownames(classmat) <- klist
colnames(classmat) <- as.character(c(1:12))

#Read in the best iterations.
bestfile = paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/',
                 'group/best_iter/',subgroup,'/best_iter.csv')
bestiter = read.csv(bestfile,row.names=1)

#Find cross k max and min LE values.
lim_LE <- c()
for (kidx in 1:knum) {
  
  #Extract.
  k = klist[kidx]
  
  #Extract best iteration.
  iteration <- as.character(bestiter[subgroup,k])
  
  #Read in the LE.
  inpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
  infile <- paste0(inpath,'subcent_',iteration,'.h5')
  inkey <- paste0('/',subgroup)
  instore <- h5py$File(infile,'r')
  clustmat <- np$array(instore[inkey])
  instore$close()
  
  #Add.
  lim_LE <- c(lim_LE,min(clustmat))
  lim_LE <- c(lim_LE,max(clustmat))
}
allk_min.LE <- min(lim_LE)
allk_max.LE <- max(lim_LE)

#Make a brain image.
base_atlas <- as.data.frame(na.omit(cbind(glasser$data$region,glasser$data$hemi)))
colnames(base_atlas) <- c('region', 'hemi')

#Read in our atlas organization.
atlasfile <- '../inputs/data/atlas/atlasorder.csv'
atlasord <- read_csv(atlasfile)[,'labelname'] %>%
  mutate_all(sub,pattern='_ROI',replacement='') %>%
  separate(col='labelname',
           sep='_',
           into=c('hemi','region')) %>%
  mutate(hemi=case_when(
    hemi=='L'~'left',
    hemi=='R'~'right')) %>%
  select('region','hemi')

#Generate layout
m = matrix(NA,nrow=knum,ncol=maxk)
pstart <- 0
for (kidx in 1:knum) {

  #Extract.
  k <- as.numeric(klist[kidx])

  #Append.
  pstart <- pstart + 1
  pend <- pstart + k - 1
  m[kidx,1:k] <- c(pstart:pend)

  #Increment.
  pstart <- pend
}

#For each hemisphere and side.
for (onehemi in c('left','right')) {
for (onelat in c('medial','lateral')) {

#For each k.
plot_list <- list()
pidx <- 1
for (kidx in 1:knum) {
  
  #Extract.
  k = klist[kidx]
  
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
  
  #For each universal index.
  for (uidx in c(1:nk)) {
    
    #Extract original state and universal state.
    origstate <- orig[uidx]
    ustate <- uval[uidx]
    
    #Read in the original state and threshold.
    kstate <- clustmat[,origstate]
    kstate[kstate < 0] <- 0
    
    #Visualize.
    vmat <- matrix(kstate)
    colnames(vmat) <- 'Loading'
    
    #Bind atlas and brain image labels.
    curr_v <- cbind(atlasord,vmat)
    atlas_v <- inner_join(base_atlas,curr_v)
    
    #Plot.
    ggobj <- atlas_v %>%
      ggplot() +
      geom_brain(atlas=glasser,
                 hemi=onehemi,
                 side=onelat,
                 aes(fill=Loading)) +
      scale_fill_scico(palette='berlin',
                       midpoint=0,
                       limits=c(allk_min.LE,allk_max.LE)) +
      theme_classic() +
      annotate(geom='text',x=-12,y=200,label=ustate,size=2.1,fontface='bold') +
      theme(axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'none',
            plot.margin=grid::unit(c(0,0,0,0),'mm')) 
    plot_list[[pidx]] <- ggobj
    pidx <- pidx + 1
  }
}

#Combine.
outpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/',
                  subgroup,'/allk/')
outfile <- paste0(outpath,'max',as.character(maxk),'_',onehemi,'_',onelat,'_ksortvis_brain.jpg')
ggsave(outfile,arrangeGrob(grobs=plot_list,layout_matrix=m),
       units='in',width=8,height=6)
}
}

#Save the color bar.
ggobj <- atlas_v %>%
  ggplot() +
  geom_brain(atlas=glasser,
             hemi=onehemi,
             side=onelat,
             aes(fill=Loading)) +
  scale_fill_scico(palette='berlin',
                   midpoint=0,
                   limits=c(0,allk_max.LE),
                   guide=guide_colorbar(frame.colour='black',ticks.colour='black',
                                        frame.linewidth=2,ticks=F)) +
  theme_classic() +
  annotate(geom='text',x=-12,y=200,label=ustate,size=2.1,fontface='bold') +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0),'mm')) 
outfile <- paste0(outpath,'max',as.character(maxk),'_ksortvis_brain_clbr.jpg')
ggsave(outfile,ggobj,units='in',width=5,height=5)
