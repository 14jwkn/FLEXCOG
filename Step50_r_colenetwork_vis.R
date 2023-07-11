# For each network, delineate it in brain ROIs by color.
# Output:
# onehemi,'_',onelat,'_colenetwork.jpg' Brain plots for networks for one hemisphere and side.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(ggseg)
library(ggsegGlasser)
library(scico)

#Set outpath.
outpath <- '../outputs/r_network/'
dir.create(outpath,recursive=T)

#Read in net labels.
netlabels <- scan('colenetlabels.txt',character(),quote='') 

#Produce namelist.
namelist <- c('Primary Visual (VIS1)','Secondary Visual (VIS2)','Somatomotor (SMN)',
            'Cingulo-Opercular (CON)','Dorsal Attention (DAN)','Language (LAN)',
            'Frontoparietal (FPN)','Auditory (AUD)','Default Mode (DMN)',
            'Posterior Multimodal (PMM)','Ventral Multimodal (VMM)',
            'Orbitoaffective (OA)')
nnet <- length(namelist)

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

#For each hemisphere and side.
for (onehemi in c('left','right')) {
for (onelat in c('medial','lateral')) {

#For each network.
plot_list <- list()
pidx <- 1
for (nidx in 1:nnet) {
  
  #Extract.
  netid <- as.character(nidx)
  netname <- namelist[nidx]
  
  #Make it 1 where the region is classified into the network.
  vmat <- matrix(0,nrow=360,ncol=1)
  vmat[netlabels == netid,] = 1
  colnames(vmat) <- 'Class'
  
  #Bind atlas and brain image labels.
  curr_v <- cbind(atlasord,vmat)
  atlas_v <- inner_join(base_atlas,curr_v)
  
  #Plot.
  ggobj <- atlas_v %>%
    ggplot() +
    geom_brain(atlas=glasser, 
               hemi=onehemi,
               side=onelat,
               aes(fill=Class)) +
    scale_fill_scico(palette='berlin',
                     midpoint=0) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0,vjust=0,size=9),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none',
          plot.margin=grid::unit(c(0,0,0,0),'mm')) +
    labs(title=netname)
  plot_list[[pidx]] <- ggobj
  pidx <- pidx + 1
}

#Combine.
outfile <- paste0(outpath,onehemi,'_',onelat,'_colenetwork.jpg')
ggsave(outfile,arrangeGrob(grobs=plot_list,ncol=4),
       units='in',width=7,height=4)
}
}
