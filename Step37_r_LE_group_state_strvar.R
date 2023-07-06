# For the subject group and the specified k, find the Pearson's correlation between
# edge-averaged dFC(t) strength and variability across all states, and within
# each state and plot corresponding scatter plots.
# Output:
# strvar.csv Strength/variability correlations across all states, and within each state.
# ysel,'_',xsel,'.jpg' Strength/variability scatter plots across all states, and within each state.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(scales)

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]
k <- args[2]

#Set base path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')

#Set parameters.
nk <- as.numeric(k)

#Read in data.
infile <- paste0(basepath,'dFC_strvar.csv')
mat_dFC <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject') 
strvar_labs <- c(colnames(mat_dFC)[2:ncol(mat_dFC)])

#Create collector.
collector <- matrix(NA,1,1+nk)
rownames(collector) <- c('dFC_R')
colnames(collector) <- c('full',1:nk)

#Correlate full values.
x1 <- 'full_mean'
x2 <- 'full_std'
collector['dFC_R','full'] <- cor(mat_dFC[,x1],mat_dFC[,x2])

#Go through each state.
for (kidx in 1:nk) {
  
  #Extract.
  kval <- as.character(kidx)
  x1 <- paste0('mean_',kval)
  x2 <- paste0('std_',kval)
  
  #Correlate.
  collector['dFC_R',kval] <- cor(mat_dFC[,x1],mat_dFC[,x2])
}

#Save strength vs variability.
outpath <- paste0(basepath,'strvar/')
dir.create(outpath,recursive=T)
outfile <- paste0(outpath,'strvar.csv')
outtab <- as.data.frame(collector)
write.csv(outtab,outfile,row.names=T)

#Get global max and mins.
allstr <- c('full_mean',paste0('mean_',1:nk))
allstd <- c('full_std',paste0('std_',1:nk))
limstr <- mat_dFC[,allstr]
limstd <- mat_dFC[,allstd]
maxstr <- max(limstr)
minstr <- min(limstr)
maxstd <- max(limstd)
minstd <- min(limstd)

#Set up limits.
xmax <- maxstd
xmin <- minstd
ymax <- maxstr
ymin <- minstr
xstep <- 0.01
ystep <- 0.01

#Set up scatter plot.
xsel <- 'full_std'
ysel <- 'full_mean'
plotmat <- cbind(mat_dFC[,xsel],mat_dFC[,ysel])
colnames(plotmat) <- c(xsel,ysel)
plotmat <- data.frame(plotmat)

#Set theme.
theme_set(theme_bw())
txtsize <- 16
acc <- 0.01

#Plot.
p1 <- ggplot(plotmat,aes_string(x=xsel,y=ysel)) +
  geom_point(size=1) +
  geom_smooth(method=lm,se=F) +
  scale_x_continuous(labels=label_number(accuracy=acc)) +
  scale_y_continuous(labels=label_number(accuracy=acc)) +
  theme(axis.text.x=element_text(size=txtsize),
        axis.text.y=element_text(size=txtsize))

#Save strength vs variability.
outpath <- paste0(basepath,'strvar/')
dir.create(outpath,recursive=T) 
outfile <- paste0(outpath,ysel,'_',xsel,'.jpg')
ggsave(outfile,p1,units='in',width=4,height=3,bg='white')

#Go through each state.
for (kidx in 1:nk) {
  
  #Extract.
  kval <- as.character(kidx)
  xsel <- paste0('std_',kval)
  ysel <- paste0('mean_',kval)
  
  #Set up scatter plot.
  plotmat <- cbind(mat_dFC[,xsel],mat_dFC[,ysel])
  colnames(plotmat) <- c(xsel,ysel)
  plotmat <- data.frame(plotmat)
  
  #Set theme.
  theme_set(theme_bw())
  
  #Plot.
  p1 <- ggplot(plotmat,aes_string(x=xsel,y=ysel)) +
    geom_point(size=1) +
    geom_smooth(method=lm,se=F) +
    scale_x_continuous(labels=label_number(accuracy=acc)) +
    scale_y_continuous(labels=label_number(accuracy=acc)) +
    theme(axis.text.x=element_text(size=txtsize),
          axis.text.y=element_text(size=txtsize))
  
  #Save.
  outfile <- paste0(outpath,ysel,'_',xsel,'.jpg')
  ggsave(outfile,p1,units='in',width=4,height=3,bg='white')
}
