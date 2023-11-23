# Gather cognitive tests and confounds into one file.
# Output:
# pred_all.csv cognitive tests and confounds (i.e., Age, Gender).

#Load packages.
library(tidyverse)

#Set outpath.
subgroup <- 'full'
outpath <- paste0('../outputs/c_cognition/',subgroup,'/')
dir.create(outpath,recursive=T)

#Read in subjects.
subjects <- scan('r_full_submain.txt',character(),quote='') 

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
  rename('PMAT24'=PMAT24_A_CR, 
         'PicVoc'=PicVocab_Unadj,
         'ProcSpeed'=ProcSpeed_Unadj,
         'CardSort'=CardSort_Unadj,
         'Flanker'=Flanker_Unadj,
         'ListSort'=ListSort_Unadj,
         'PicSeq'=PicSeq_Unadj,
         'ReadVoc'=ReadEng_Unadj,
         'WordMem'=IWRD_TOT,
         'LineOrient'=VSPLOT_TC,
         'Age'=Age_in_Yrs) %>%
  mutate_if(is.character,as.factor) %>%
  mutate(Gender=recode_factor(Gender,
                              'F'='F',
                              'M'='M'))

#Save the matrix.
write.csv(hcpselect,paste0(outpath,'pred_all.csv'),row.names=F)
