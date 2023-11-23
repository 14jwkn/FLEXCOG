# Plot the cluster validity index for each k to produce the elbow plot.
# Output:
# curr_lab,'.jpg' Cluster validity index elbow plot.

#Load packages.
library(tidyverse)
library(reticulate)
library(ggplot2)
library(gridExtra)

#Set virtual environment.
use_virtualenv('../envs/FlexEnv')

#Load python packages.
h5py <- import('h5py')
np <- import ('numpy')
pd <- import('pandas')

#Catch arguments.
args = commandArgs(trailingOnly=T)
subfile <- args[1]
print(paste0('Doing: ',subfile))

#Set initial values.
theme_set(theme_classic())
nclust <- c(2:12)
maxk <- max(nclust)
subjects <- scan(subfile,character(),quote='')

#Define inpaths.
subgroup = 'full'
inkey <- paste0('/',subgroup)
inpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                 subgroup,'/')
groupfile <- paste0(inpath,'metric_group.h5')

#Read in metric values.
store <- h5py$File(groupfile,'r')
group <- np$array(store[inkey])
store$close()

#Define outpath.
outpath <- paste0(inpath,'elbow_plots/')
dir.create(outpath,recursive=T)

#Plot group versions.
met <- group[,2]
plotmat <- tibble('k'=nclust,'CVI_Metric'=met)
subplot <- plotmat %>%
  ggplot(aes(x=k,y=CVI_Metric)) +
  geom_line() +
  scale_x_continuous(breaks=c(1:maxk)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
subplot
curr_lab <- paste0('cvi_out')
outfile <- paste0(outpath,curr_lab,'.jpg')
ggsave(outfile,units='in',width=4,height=4)
print('Saved output.')
