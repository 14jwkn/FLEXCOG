# For the subject group, gather g, cognitive tests used to find g, and confounds
# into one file.
# Output:
# pred_all.csv g, cognitive tests used to find g, and confounds (i.e. Age, Gender).

#Load packages.
library(tidyverse)

#Catches arguments.
args = commandArgs(trailingOnly=T)
subgroup <- args[1]
print(paste0('Doing: ',subgroup))

#Set outpath.
outpath <- paste0('../outputs/c_cognition/',subgroup,'/')
dir.create(outpath,recursive=T)

#Read in subjects.
if (subgroup == 'full') {
  subjects <- scan('r_full_submain.txt',character(),quote='') 
} else if (subgroup == 'half') {
  subjects <- scan('r_half_submain.txtt',character(),quote='') 
} 

#Read in g.
infile <- paste0('../outputs/c_cognition/'+subgroup+'/cogscores.csv')
gmat <- read.csv(infile,header=T,row.names=1)
gmat <- gmat[subjects,]
gmat <- rownames_to_column(gmat,'Subject') %>%
  rename('g'=gPCA) %>%
  select(Subject,g) %>%
  mutate(Subject=as.character(Subject))

#Read in other cognitive variables and predictors.
inpath <- '../inputs/data/hcp/'
infile <- paste0(inpath,'HCP1200_Data_Dictionary.csv')
hcpopen <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject))
infile <- paste0(inpath,'RESTRICTED_HCP1200_Data_Dictionary.csv')
hcpclose <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject))
infile <- paste0(inpath,'unrestricted_hcp_freesurfer.csv')
hcpfs <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject))
hcpfull <- hcpopen %>%
  full_join(hcpclose) %>%
  full_join(hcpfs)
hcpselect <- hcpfull %>%
  select(Subject,
         PMAT24_A_CR, 
         PicVocab_Unadj,
         ProcSpeed_Unadj,
         CardSort_Unadj,
         Flanker_Unadj,
         ListSort_Unadj,
         PicSeq_Unadj,
         ReadEng_Unadj,
         IWRD_TOT,
         VSPLOT_TC,
         Gender,
         Age_in_Yrs) %>%
  rename('P24_CR'=PMAT24_A_CR, 
         'PV'=PicVocab_Unadj,
         'PS'=ProcSpeed_Unadj,
         'CardSort'=CardSort_Unadj,
         'Flanker'=Flanker_Unadj,
         'ListSort'=ListSort_Unadj,
         'PicSeq'=PicSeq_Unadj,
         'ReadEng'=ReadEng_Unadj,
         'Age'=Age_in_Yrs) %>%
  mutate_if(is.character,as.factor) %>%
  mutate(Gender=recode_factor(Gender,
                              'F'='F',
                              'M'='M'))

#Merge matrices.
predall <- gmat %>%
  inner_join(hcpselect) 

#Save the matrix.
write.csv(predall,paste0(outpath,'pred_all.csv'),row.names=F)

