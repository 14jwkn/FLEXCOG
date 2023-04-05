# For the subject group and the specified k, read in the strength/variability
# data each state and rank subjects according to the average values.
# Output:
# 'high_',cset,'.txt' Participants ranked from high to low in average scores.
# 'low_',cset,'.txt' Participants ranked from low to high in average scores.

#Set libraries.
library(tidyverse)
library(ggplot2)

#Catches arguments.
subgroup <- args[1]
k <- args[2]

#Set base path and set path,
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
setpath <-  paste0(basepath,'allreconfig_subsets/strvar/')
dir.create(setpath,recursive=T)

#Set parameters.
nk <- as.numeric(k)
limval <- 1

#Get strength/variability data.
infile <- paste0(basepath,'dFC_strvar.csv')
mat_dFC <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject') 

# Rank Subjects ----------------------------------------------------------

#Go through each subset.
statetype_list <- c('dFC')
nstatetype <- length(statetype_list)
strvar_list <- c('mean','std')
nstrvar <- length(strvar_list)
for (kidx in 1:nk) {
  ck <- as.character(kidx)
  for (stidx in 1:nstatetype) {
    cstatetype <- statetype_list[stidx]
    cmat <- mat_dFC
    for (svidx in 1:nstrvar) {
      cstrvar <- strvar_list[svidx]
      clab <- paste0(cstrvar,'_',ck)
      
      #Extract and rank.
      crank <- cmat[,c('Subject',clab)]
      colnames(crank) <- c('Subject','Value')
      crank <- crank %>%
        arrange(desc(Value))
      highsub <- as.character(as.matrix(crank[,'Subject']))
      crank <- crank %>%
        arrange(Value)
      lowsub <- as.character(as.matrix(crank[,'Subject']))
      
      #Save text files.
      cset <- paste0(cstatetype,'_',clab)
      outpath <- setpath
      outfile <- file(paste0(outpath,'high_',cset,'.txt'))
      writeLines(highsub,outfile)
      close(outfile)
      outfile <- file(paste0(outpath,'low_',cset,'.txt'))
      writeLines(lowsub,outfile)
      close(outfile)
    }
  }
}
