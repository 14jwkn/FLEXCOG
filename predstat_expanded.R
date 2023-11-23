# Gather cognitive tests and confounds into one file, with the addition of
# gPCA, gEFA, and gCFA.
# Output:
# pred_all.csv Cognitive tests and confounds (i.e., Age, Gender), with the addition of gPCA, gEFA, and gCFA.

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

#Read in gPCA.
infile <- paste0('../outputs/c_cognition/full/PCA_scores.csv')
gPCA <- read.csv(infile,header=T,row.names=1)
gPCA <- gPCA[subjects,1]
gPCA <- rownames_to_column(gPCA,'Subject') 
colnames(gPCA) <- c('Subject','gPCA')
gPCA <- gPCA %>%
  mutate(Subject=as.character(Subject))
hcpselect <- hcpselect %>%
  left_join(gPCA)

#Read in gEFA.
infile <- '../outputs/c_cognition/ourfa_subs/dubois_efa_scores.csv'
du_EFA <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject)) %>%
  select(Subject,gEFA)
hcpselect <- hcpselect %>%
  left_join(du_EFA)

#Read in gCFA.
infile <- '../outputs/c_cognition/ourfa_subs/dubois_cfa_scores.csv'
du_CFA <- read_csv(infile) %>%
  mutate(Subject=as.character(Subject)) %>%
  select(Subject,gCFA)
hcpselect <- hcpselect %>%
  left_join(du_CFA)

#Save the matrix.
write.csv(hcpselect,paste0(outpath,'pred_all.csv'),row.names=F)
