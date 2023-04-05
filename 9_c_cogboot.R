# For the subject group, do bootstrap resampling for the g PC1 loadings.
# Output:
# gPCA_load.csv Loadings for cognitive tests on PC1.
# gPCA_load_bootstab.csv Bootstrap ratios for cognitive tests on PC1.
# gPCA_load_bootstab.jpg Plots of loadings colored based on bootstrap ratios for PC1.


#Set libraries.
library(tidyverse)
library(InPosition)
library(PTCA4CATA)

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read.csv(infile,header=T,row.names=1) 
othercog <- othercog[subjects,]

#Isolate desired cognitive variables.
ourcog <- colnames(othercog)
matchlist <- c('CardSort','Flanker','ListSort','PicSeq','PV','PS','ReadEng','P24_CR',
               'IWRD_TOT','VSPLOT_TC')
ourcog <- ourcog[ourcog %in% matchlist]
cogmat <- othercog[,ourcog] %>%
  rename('PMAT24'='P24_CR',
         'PicVoc'='PV',
         'ProcSpeed'='PS',
         'WordMem'='IWRD_TOT',
         'LineOrient'='VSPLOT_TC')

#Do PCA.
set.seed(12345)
pcamat <- epPCA.inference.battery(cogmat,scale=T,center=T,
                                  k=10,
                                  test.iters=10000,critical.value=2.5)
pcascore <- -(pcamat$Fixed.Data$ExPosition.Data$fi[,1])
pcaload <- -(pcamat$Fixed.Data$ExPosition.Data$fj[,1])
pcaboot <- -(pcamat$Inference.Data$fj.boots$tests$boot.ratios[,1])
pcaboot.sig <- pcamat$Inference.Data$fj.boots$tests$sig.boot.ratios[,1]

#Save the matrices.
pcaload.pack <- data.frame(pcaload)
pcaboot.pack <- data.frame(pcaboot)
outpath <- '../outputs/c_cognition/'+subgroup+'/'
write.csv(pcaload.pack,paste0(outpath,'gPCA_load.csv'),row.names=T)
write.csv(pcaboot.pack,paste0(outpath,'gPCA_load_bootstab.csv'),row.names=T)

#Set colors.
norm.color <- '#696969'
pos.color <- '#FF0000'
neg.color <- '#4E79A7'

#Plot.
collab <- pcaload
collab[pcaload > 0] <- pos.color
collab[pcaload < 0] <- neg.color
collab[!(pcaboot.sig)] <- norm.color
plty <- PrettyBarPlot2(pcaload,
                      threshold=0, 
                      color4bar=collab,
                      font.size=3) 
outfile <- paste0(outpath,'gPCA_load_bootstab.jpg')
ggsave(outfile,plty,units='in',width=2,height=3)
