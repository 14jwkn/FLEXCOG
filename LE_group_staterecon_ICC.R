# For the specified k, find ICC(1,1) between the session averages across runs
# for each network reconfiguration metric.
# Output:
# cout,'.jpg' ICC matrix plots for transition distance, transition probability, k-wise metrics, transition number, and a sample matrix for formatting. Also, ICC interpretation legend.

#Set libraries.
library(tidyverse)
library(psych)

#Catches arguments.
args <- commandArgs(trailingOnly=T)
k <- args[1]

#Set base path and output path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
outpath <-  paste0(basepath,'allreconfig/ICC/')
dir.create(outpath,recursive=T)

#Set parameters.
nk <- as.numeric(k)

#Set runs and sessions.
runs <- c('REST1_LR','REST1_RL','REST2_LR','REST2_RL')
runlabs <- c('r1_lr','r1_rl','r2_lr','r2_rl')
ses <- c('ses1','ses2')
ses_runs <- list(c('REST1_LR','REST1_RL'),c('REST2_LR','REST2_RL'))
ses_runlabs <- list(c('r1_lr','r1_rl'),c('r2_lr','r2_rl'))
nrun <- length(runs)
nses <- length(ses)

#Get session averages for dynamics, state similarity, and state jump.
reconlist <- list()
for (sidx in 1:nses) {
  cses <- ses[sidx]
  cses_out <- paste0('s',as.character(sidx))
  cses_runs <- ses_runs[[sidx]]
  cses_runlabs <- ses_runlabs[[sidx]]
  nses_runs <- length(cses_runs)
  sescollect <- matrix(0,nrow=nsub,ncol=nfeat)
  colnames(sescollect) <- featlabs
  for (ridx in 1:nses_runs) {
    crun <- cses_runs[ridx]
    crun_lab <- cses_runlabs[ridx]
    
    #Get predictor data.
    infile <- paste0(basepath,crun_lab,'_substat.csv')
    dynmat <- read_csv(infile,col_names=T) %>% 
      rename_at(1,~'Subject')
    colnames(dynmat) <- gsub('sptrans_','transpro_',colnames(dynmat))
    
    infile <- paste0(basepath,'statesim/',crun,'_meansim.csv')
    simmat <- read_csv(infile,col_names=T) %>% 
      rename_at(1,~'Subject')
    colnames(simmat) <- gsub('State','idio',colnames(simmat))
    
    infile <- paste0(basepath,'statejump/',crun,'_transmean.csv')
    jumpmat <- read_csv(infile,col_names=T) %>% 
      rename_at(1,~'Subject')
    colnames(jumpmat) <- c('Subject',paste0('transdist_',colnames(jumpmat)[2:ncol(jumpmat)]))
    predmat <- dynmat %>%
      inner_join(simmat) %>%
      inner_join(jumpmat)
    colnames(predmat) <- gsub('-','t',colnames(predmat))
    nsub <- nrow(predmat)
    featlabs <- colnames(predmat)[-1]
    nfeat <- length(featlabs)
    
    #Add run to average.
    sescollect <- sescollect + predmat[,-1]
  }
  
  #Divide by number of runs to average, add back subject column.
  sescollect <- sescollect/nses_runs
  sescollect <- cbind(predmat[,'Subject'],sescollect)
  
  #Append.
  reconlist[[cses_out]] <- sescollect
}

#Calculate the ICC. We want ICC(1,1).
icc_collect <- matrix(NA,nfeat,1)
rownames(icc_collect) <- featlabs
colnames(icc_collect) <- 's1vs2'
cmat1 <- reconlist[['s1']]
cmat2 <- reconlist[['s2']]

#For each feature.
for (fidx in 1:nfeat) {
  cfeat <- featlabs[fidx]
  ccom <- cbind(cmat1[,cfeat],cmat2[,cfeat])
  
  #Calculate ICC. 
  cobj <- ICC(ccom)[['results']]
  icc_collect[cfeat,'s1vs2'] <- cobj[cobj[,'type']=='ICC1','ICC']
}

#Convert to interpretations.
icc_collect_i <- data.frame(icc_collect)
ccomp <- icc_collect
icc_collect_i[ccomp<0] <- 'Negative'
icc_collect_i[(0<=ccomp)&(ccomp<0.2)] <- 'Low'
icc_collect_i[(0.2<=ccomp)&(ccomp<0.4)] <- 'Fair'
icc_collect_i[(0.4<=ccomp)&(ccomp<0.6)] <- 'Moderate'
icc_collect_i[(0.6<=ccomp)&(ccomp<0.8)] <- 'Substantial'
icc_collect_i[(0.8<=ccomp)&(ccomp<=1)] <- 'Perfect'

# Output ------------------------------------------------------------------

#Generate tables for transition probability and distance.
pro_collect <- data.frame(matrix(NA,nrow=nk,ncol=nk))
rownames(pro_collect) <- paste0('s',as.character(1:nk))
colnames(pro_collect) <- paste0('s',as.character(1:nk))
pro_collect_i <- pro_collect
dist_collect <- pro_collect
dist_collect_i <- pro_collect
for (kidx1 in 1:nk) {
  for (kidx2 in 1:nk) {
    
    #Probability.
    cpro <- paste0('transpro_',as.character(kidx1),'t',as.character(kidx2))
    pro_collect[kidx1,kidx2] <- icc_collect[cpro,'s1vs2']
    pro_collect_i[kidx1,kidx2] <- icc_collect_i[cpro,'s1vs2']
    
    #Distance.
    cdist <- paste0('transdist_',as.character(kidx1),'t',as.character(kidx2))
    dist_collect[kidx1,kidx2] <- icc_collect[cdist,'s1vs2']
    dist_collect_i[kidx1,kidx2] <- icc_collect_i[cdist,'s1vs2']
  }
}

#k metrics.
k_met <-  data.frame(matrix(NA,nrow=3,ncol=nk))
rownames(k_met) <- c('Occurrence','Dwell Time','Idiosyncrasy')
colnames(k_met) <- paste0('s',as.character(1:nk))
k_met_i <- k_met

#Fill in.
k_met['Occurrence',] <- icc_collect[grepl('occur',rownames(icc_collect)),'s1vs2']
k_met['Dwell Time',] <- icc_collect[grepl('dwell',rownames(icc_collect)),'s1vs2']
k_met['Idiosyncrasy',] <- icc_collect[grepl('idio',rownames(icc_collect)),'s1vs2']
k_met_i['Occurrence',] <- icc_collect_i[grepl('occur',rownames(icc_collect)),'s1vs2']
k_met_i['Dwell Time',] <- icc_collect_i[grepl('dwell',rownames(icc_collect)),'s1vs2']
k_met_i['Idiosyncrasy',] <- icc_collect_i[grepl('idio',rownames(icc_collect)),'s1vs2']

#Transition number.
nt <- data.frame('Transition Number'=icc_collect['numtrans','s1vs2'])
nt_i <- data.frame('Transition Number'=icc_collect_i['numtrans','s1vs2'])

# Generate Plot -----------------------------------------------------------

#Set values.
theme_set(theme_classic())
ilevels <- c('Negative','Low','Fair','Moderate','Substantial','Perfect')
icolors <- c('#440154FF','#414487FF','#2A788EFF','#22A884FF','#7AD151FF','#FDE725FF')
names(icolors) <- ilevels
klabs <- paste0('S',1:nk)
dynlabs <- c('OC','DT','IS')
cmat_list <- list(dist_collect,pro_collect,k_met,nt,dist_collect)
cmat_i_list <- list(dist_collect_i,pro_collect_i,k_met_i,nt_i,dist_collect_i)
row_list <- list(klabs,klabs,dynlabs,c('TN'),klabs)
col_list <- list(klabs,klabs,klabs,c('TN'),klabs)
out_row_list <- list(klabs,klabs,dynlabs,c('TN'),rep('TN',times=nk))
out_col_list <- list(klabs,klabs,klabs,c('TN'),rep('TN',times=nk))
dim_list <- list(c(3,2,3),c(3,2,3),c(3,1,3),c(1,(2/3),4),c(3,2,3))
out_list <- list('transdist','transpro','misc','numtrans','sample')
nplots <- length(cmat_list)
for (pltidx in 1:nplots) {
  
  #Collect.
  cmat <- cmat_list[[pltidx]]
  cmat_i <- cmat_i_list[[pltidx]]
  crow_labs <- row_list[[pltidx]]
  ccol_labs <- col_list[[pltidx]]
  cdim <- dim_list[[pltidx]]
  cout <- out_list[[pltidx]]
  out_row_labs <- out_row_list[[pltidx]]
  out_col_labs <- out_col_list[[pltidx]]
  
  #Plot.
  rown <- nrow(cmat)
  coln <- ncol(cmat)
  outmat <- data.frame(matrix(NA,rown*coln,4))
  colnames(outmat) <- c('Row','Column','Value','Interpretation')
  curridx <- 1
  for (rowidx in 1:rown) {
    for (colidx in 1:coln) {
      outmat[curridx,1] <- crow_labs[rowidx]
      outmat[curridx,2] <- ccol_labs[colidx]
      outmat[curridx,3] <- format(round(cmat[rowidx,colidx],2),nsmall=2)
      outmat[curridx,4] <- cmat_i[rowidx,colidx]
      curridx <- curridx + 1
    }
  }
  outmat[,1] <- factor(outmat[,1],levels=crow_labs)
  outmat[,2] <- factor(outmat[,2],levels=ccol_labs)
  outmat[,4] <- factor(outmat[,4],levels=ilevels)
  p1 <- ggplot(outmat,aes(Column,Row,fill=Interpretation)) +
    geom_tile(color='black',size=0.5) +
    scale_fill_manual(values=icolors) +
    geom_text(aes(label=Value),size=cdim[3],fontface='bold') +
    scale_x_discrete(expand=c(0,0),labels=out_col_labs) +
    scale_y_discrete(limits=rev,expand=c(0,0),labels=rev(out_row_labs)) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          legend.position='none') 
  
  #Save.
  outfile <- paste0(outpath,cout,'.jpg')
  ggsave(outfile,p1,units='in',width=cdim[1],height=cdim[2])
}

#Save legend.
ilevels2 <- c('Negative (ICC < 0)',
              'Low (0 < ICC < 0.2)',
              'Fair (0.2 < ICC < 0.4)',
              'Moderate (0.4 < ICC < 0.6)',
              'Substantial (0.6 < ICC < 0.8)',
              'Perfect (0.8 < ICC < 1)')
icolors2 <- c('#440154FF','#414487FF','#2A788EFF','#22A884FF','#7AD151FF','#FDE725FF')
names(icolors2) <- ilevels2
cout <- 'interpretation_legend'
dummat <- data.frame('Column'=as.character(c(1:length(icolors2))),
                     'Row'=as.character(c(1:length(icolors2))),
                     'Interpretation'=factor(ilevels2,levels=rev(ilevels2)))
p1 <- ggplot(dummat,aes(Column,Row,fill=Interpretation)) +
  geom_tile(color='black') +
  scale_fill_manual(values=icolors2) +
  theme(legend.key=element_rect(color='black',size=1),
        legend.text=element_text(face='bold'))
outfile <- paste0(outpath,cout,'.jpg')
ggsave(outfile,p1,units='in',width=3,height=3)
