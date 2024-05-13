# For the specified k, conduct PLSC where X = frequency metrics and Y = 10 cognitive 
# tests of interest. Do PLSC to find direction of loadings to flip, conduct 
# permutation testing on singular values, investigate bivariate normality and 
# outliers for significant latent variable pairs and correlate, 
# conduct bootstrap resampling on loadings, investigate reproducibility of
# singular values and loadings, investigate out-of-sample results from cross-validation,
# and variance explained by in-sample and out-of-sample latent variables.
# Adapts code from:
# R package: data4PCCAR
# https://github.com/McIntosh-Lab/PLS-CCA
# https://github.com/ThomasYeoLab/Standalone_Ooi2022_MMP
# https://github.com/juchiyu/Multivariate_functions/tree/main
# Output:
# main_var_explained.csv Table containing variance explained for the main PLSC latent variables.
# permsig.csv Table containing singular value, covariance explained, and permutation p-values.
# permsig.jpg Plot indicating which dimensions are significant.
# 'baseLVY',llab,'.jpg' Bar plot of Y cognitive loadings with default directions.
# yflip.csv Table indicating whether to flip the dimension for interpretation.
# 'LV',as.character(i),'_corr.jpg' Scatterplot for significant latent variable pairs.
# 'LV',as.character(i),'_qq.jpg' Bivariate QQ plot for significant latent variable pairs.
# binormtest_p.csv Bivariate normality test p-values for significant latent variable pairs.
# binormtest_pthres.csv Bivariate normality test p-values for significant latent variable pairs thresholded to significant ones.
# LVX_scores.csv Significant X latent variable scores.
# LVY_scores.csv Significant Y latent variable scores.
# 'LVY',llab,'_load.jpg' Bar plot of Y cognitive loadings where stable loadings are colored for one dimension.
# 'LVX',llab,'_load_',namedyn,'.jpg' Bar plot of X reconfiguration loadings where stable loadings are colored for one dimension.
# LVYall_load.jpg Bar of Y cognitive loadings for all dimensions with stable loadings colored.
# LVXall_load.jpg Bar of X reconfiguration loadings for all dimensions with stable loadings colored.
# LVXall_load.csv X reconfiguration loadings thresholded by bootstrap ratio.
# LVYall_load.csv Y cognitive loadings thresholded by bootstrap ratio.
# LVXall_load_full.csv X reconfiguration loadings.
# LVYall_load_full.csv Y cognitive loadings.
# 'LVY',llab,'_bootrat.jpg' Bar plot of Y cognitive bootstrap ratios where stable loadings are colored for one dimension.
# 'LVX',llab,'_bootrat_',namedyn,'.jpg' Bar plot of X reconfiguration bootstrap ratios where stable loadings are colored for one dimension.
# LVYall_bootrat.jpg Bar of Y cognitive bootstrap ratios for all dimensions with stable loadings colored.
# LVXall_bootrat.jpg Bar of X reconfiguration bootstrap ratios for all dimensions with stable loadings colored.
# LVXall_bootrat.csv X reconfiguration bootstrap ratios thresholded by bootstrap ratio.
# LVYall_bootrat.csv Y cognitive bootstrap ratios thresholded by bootstrap ratio.
# LVXall_bootrat_full.csv X reconfiguration bootstrap ratios.
# LVYall_bootrat_full.csv Y cognitive bootstrap ratios.
# olab,'.csv' Mean and Z-score reproducibility values for singular values and loadings.
# replicate_var_explained.csv Table containing variance explained for the test PLSC latent variables from reproducibility analysis splits.
# svrep.jpg Box plot of singular value reproducibility across split-halfs colored based on Z-score.
# lvxrep.jpg Box plot of X loading reproducibility across split-halfs colored based on Z-score.
# lvyrep.jpg Box plot of Y loading reproducibility across split-halfs colored based on Z-score.
# repsummary.jpg Box plots of all reproducibility scores for significant dimensions.
# cv_var_explained.csv Table containing variance explained for the test PLSC latent variables from cross-validation analysis splits.
# cvout.csv Table containing correlations between test latent variables, and between original and test latent variables for all dimensions from cross-validation.
# 'LV',as.character(i),'_LXHLYH.jpg' Scatter plot between test latent variables for one dimension from cross-validation.
# 'LV',as.character(i),'_LXLXH.jpg' Scatter plot between original and test reconfiguration latent variables for one dimension from cross-validation.
# 'LV',as.character(i),'_LYLYH.jpg' Scatter plot between original and test cognition latent variables for one dimension from cross-validation.
# 'LV',as.character(i),'_',c_XY_lab,'_qq.jpg' QQ plot for the bivariate relationship between test latent variables, and between original and test latent variables for one dimension from cross-validation. 
# 'binormtest_',outlab,'_p.csv' Table containing p-values for normality tests of the bivariate relationship between test latent variables, and between original and test latent variables for all significant dimensions from cross-validation.
# 'binormtest_',outlab,'_pthres.csv' Table containing p-values for normality tests of the bivariate relationship between test latent variables, and between original and test latent variables for all significant dimensions from cross-validation, thresholded for significance.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(TExPosition)
library(PTCA4CATA)
library(data4PCCAR)
library(MVN)

# Helper Functions --------------------------------------------------------

#Residualizer function.
resid_maker <- function(cvar,cont_list,cat_list,cmat,resflag,coef_list=c()) {
  
  #If built-in R function without intercept (demeaning).
  if (resflag == 'demean') {
    
    #Do test with built-in R function.
    conform <- paste0(c(cont_list,cat_list),collapse='+')
    inform <- as.formula(paste0(cvar,'~',conform))
    valresid <- data.frame(resid(lm(inform,data=cmat))) 
    colnames(valresid) <- cvar
    
    #Return the true - predicted value without labels.
    return(valresid)
    
    #If raw value with intercept kept (no demeaning).
  } else if (resflag == 'raw') {
    
    #Isolate target variable.
    valtrue <- cmat %>%
      select(all_of(cvar))
    selmat <- valtrue
    
    #Isolate categorical variables, produce dummy variables, remove reference column.
    if (length(cat_list) != 0) {
      catmat <- cmat %>%
        select(all_of(cat_list)) %>%
        mutate_if(is.character,as.factor) 
      selcat <- as.data.frame(model.matrix(~.,data=catmat))
      selcat[,'(Intercept)'] <- NULL
      selmat <- cbind(selmat,selcat)
    }
    
    #Isolate continuous variables.
    if (length(cont_list) != 0) {
      selcont <- cmat %>%
        select(all_of(cont_list))
      selmat <- cbind(selmat,selcont)
    }
    
    #Fit.
    resid_out <- colnames(selmat)[colnames(selmat)!=cvar]
    resid_form <- paste0(resid_out,collapse='+')
    inform <- as.formula(paste0(cvar,'~',resid_form))
    cfit <- lm(inform,data=selmat)
    
    #If coefficients are not provided.
    if (length(coef_list)==0) {
      
      #Use confound coefficients to build predicted value from confounds.
      valpred <- 0
      for (cconf in resid_out) {
        valpred <- valpred + cfit$coefficients[cconf]*selmat[,cconf]
      }
      outcoef <- cfit$coefficients #Return coefficients.
    
    #If coefficients are provided.
    } else {
      
      #Use confound coefficients to build predicted value from confounds.
      valpred <- 0
      for (cconf in resid_out) {
        valpred <- valpred + coef_list[cconf]*selmat[,cconf]
      }
      outcoef <- coef_list #Return coefficients.
    }
    
    #Produce true - predicted.
    valresid <- data.frame(valtrue-valpred)
    
    #Return the true - predicted value without labels and confound coefficients.
    return(list(valresid,outcoef))
  }
}

#Permutation with family structure.
permF <- function(X,DATAF) {
  
  #Convert family IDs to family indices, get number of unique families, and get
  #the size of each family.
  famidx <- as.integer(factor(DATAF))
  unique_nfam <- max(famidx)
  famsize <- matrix(NA,unique_nfam,1)
  for (fidx in 1:unique_nfam) {
    famsize[fidx,1] <- sum(famidx==fidx)
  }
  
  #Collect permutation groups with unrelated individuals.
  nsub <- length(X)
  perm_group <- matrix(NA,nsub,1)
  
  #For each family.
  for (ufidx in 1:unique_nfam) {
    
    #For the family, get the count of subjects. From 1 to the max number of
    #subjects per family, sample a random permutation group index for each subject 
    #in the family.
    sorting_order <- sample(1:max(famsize),size=famsize[ufidx,1],replace=F)
    
    #Record the permutation group indices for the subjects of the family.
    perm_group[famidx==ufidx,] <- sorting_order
  }
  
  #Set up indices and permuted indices.
  true_idx <- matrix(1:nsub)
  perm_idx <- matrix(NA,nsub)
  
  #For each permutation group, permute subject indices within group.
  for (groupidx in 1:max(famsize)) {
    tmp_idx <- true_idx[perm_group==groupidx]
    tmp_idx <- sample(tmp_idx,replace=F)
    perm_idx[perm_group==groupidx] <- tmp_idx
  }
  
  #Apply the permutation to the values and return.
  permX <- X[perm_idx]
  return(permX)
}

#Bootstrapping with family structure.
bootF <- function(sublist,famlist) {
  
  #Get length.
  nsub <- length(sublist)
  
  #Create family structure.
  unique_famlist <- unique(famlist) #Get unique family IDs.
  nfam <- length(unique_famlist) #Get number of unique families.
  sub_perfam <- list() #Split individuals into families.
  for (fidx in 1:nfam) {
    subidx_list <- sublist[famlist==unique_famlist[fidx]] #Find individuals of the current family.
    sub_perfam[[fidx]] <- subidx_list #Append.
  }
  
  #Generate collector for individual IDs.
  nsub <- length(famlist)
  subidx_collect <- c()
  
  #Sample family IDs with replacement.
  famidx_rand <- sample(1:nfam,replace=T) 
  
  #Retrieve each family.
  for (fidx in 1:nfam) {
    
    #Retrieve subjects for the family ID.
    fidx_rand <- famidx_rand[fidx]
    famsubs <- sub_perfam[[fidx_rand]] 
    
    #Put them in.
    subidx_collect <- c(subidx_collect,famsubs)
  }
  
  #Return the individual IDs.
  return(subidx_collect)
}

#Cross-validation with family structure.
cvF <- function(sublist,famlist,nfold) {
  
  #Get length.
  nsub <- length(sublist)
  
  #Create family structure.
  unique_famlist <- unique(famlist) #Get unique family IDs.
  nfam <- length(unique_famlist) #Get number of unique families.
  sub_perfam <- list() #Split individuals into families.
  for (fidx in 1:nfam) {
    subidx_list <- sublist[famlist==unique_famlist[fidx]] #Find individuals of the current family.
    sub_perfam[[fidx]] <- subidx_list #Append.
  }
  
  #Initialize list of empty folds.
  fold_list <- list() 
  for (cfold in 1:nfold) {
    fold_list[[cfold]] <- vector()
  }
  
  #Split families into folds.
  subfold_size <- ceiling(nsub/nfold) #Find number of subjects that should go in each test fold approximately.
  famidx_rand <- sample(1:nfam,replace=F) #Randomize list of family IDs.
  cfold <- 0 #Current fold initialize.
  for (fidx in 1:nfam) {
    
    #Set fold to put the family in.
    cfold <- (cfold%%nfold) + 1 
    
    #If the number of individuals in the fold is greater than the desired size, increment fold.
    while (length(fold_list[[cfold]]) >= subfold_size) { 
      cfold <- (cfold%%nfold) + 1
    }
    
    #Put a random family into the fold.
    fidx_rand <- famidx_rand[fidx]
    fold_list[[cfold]] <- c(fold_list[[cfold]],sub_perfam[[fidx_rand]])
  }
  
  #Generate a table containing the individuals and which test fold they are in and return.
  testout <- matrix(NA,nsub)
  rownames(testout) <- sublist
  for (foidx in 1:nfold) {
    cfold <- fold_list[[foidx]]
    testout[cfold,] <- foidx
  }
  return(testout)
}

# Run Functions -----------------------------------------------------------

#PLSC permutation with family structure.
perm4PLSC_F <- function(DATA1,DATA2,DATAF,center1=TRUE,center2=TRUE,scale1="ss1", 
                        scale2="ss1",nIter=1000,permType="byMat",compact=FALSE) {
  if (permType != "byColumns") 
    permType <- "byMat"
  DATA1 <- as.matrix(DATA1)
  DATA2 <- as.matrix(DATA2)
  X = ExPosition::expo.scale(DATA1, center1, scale1)
  Y = ExPosition::expo.scale(DATA2, center2, scale2)
  if (NCOL(X) > NCOL(Y)) { #Make Y matrix the one with the most rows.
    tmpoX <- X
    X <- Y
    Y <- tmpoX
  }
  nN <- NROW(X)
  nI <- NCOL(X)
  nJ <- NCOL(Y)
  if (!(nN == NROW(Y))) {
    stop("DATA1 and DATA2 non-conformable")
  }
  maxRank <- min(nI, nJ) #Maximum number of features.
  Sfixed = compS(X, Y, center1 = FALSE, center2 = FALSE, scale1 = FALSE, 
                 scale2 = FALSE) #Compute SCP matrix.
  fixedEigenvalues <- rep(0, maxRank)
  fixedEV <- sv2(Sfixed) #Compute squared singular values of the matrix.
  if (length(fixedEV) > maxRank) { #Depending on the maximum number of features, set the squared singular values we use later.
    fixedEigenvalues <- fixedEV[1:maxRank]
  }
  if (length(fixedEV) == maxRank) {
    fixedEigenvalues <- fixedEV
  }
  if (length(fixedEV) < maxRank) {
    fixedEigenvalues[1:length(fixedEV)] <- fixedEV
  }
  fixedInertia <- sum(fixedEigenvalues) #Sum of singular values.
  permInertia <- rep(NA, nIter) #Collect permutations of the sum of singular values.
  permEigenvalues <- matrix(NA, nrow = nIter, ncol = maxRank) #Collect permutations of the singular values.
  .truc <- function(X, Y, longueur = min(c(dim(X), NCOL(Y))), 
                    permType = permType) {
    valP <- rep(0, longueur)
    if (permType == "byMat") { #Permute by matrix.
      Xrand <- X[permF(1:nN,DATAF), ] #Shuffle X only.
      Yrand <- Y
    }
    if (permType == "byColumns") { #Permute by column.
      Xrand <- apply(X, 2, permF, DATAF) #Shuffle X.
      Yrand <- apply(Y, 2, permF, DATAF) #Shuffle Y.
    }
    Srand <- compS(Xrand, Yrand, FALSE, FALSE, FALSE, FALSE)
    resvp <- sv2(Srand)
    valP[1:length(resvp)] <- resvp
    return(valP)
  }
  laLongueur <- maxRank + 1
  permEigenvalues <- replicate(nIter, .truc(X, Y, laLongueur, 
                                            permType)) #Repeatedly run the function.
  permEigenvalues <- t(permEigenvalues[1:maxRank, ])
  permInertia = rowSums(permEigenvalues)
  pOmnibus = sum(permInertia > fixedInertia)/nIter
  if (pOmnibus == 0) 
    pOmnibus <- 1/nIter
  pEigenvalues <- rowSums(t(permEigenvalues) > (fixedEigenvalues))/nIter
  pEigenvalues[pEigenvalues == 0] <- 1/nIter
  return.list <- structure(list(fixedInertia = fixedInertia, 
                                fixedEigenvalues = fixedEigenvalues, pOmnibus = pOmnibus, 
                                pEigenvalues = pEigenvalues), class = "perm4PLSC")
  if (!compact) {
    return.list$permInertia = permInertia
    return.list$permEigenvalues = permEigenvalues
  }
  return(return.list)
}

#PLSC bootstrapping with family structure.
Boot4PLSC_F <- function (DATA1,DATA2,sublist,famlist,center1=TRUE,center2=TRUE,scale1 = "ss1", 
                         scale2="ss1",Fi=NULL,Fj=NULL,nf2keep=3,nIter=1000, 
                         critical.value=2,eig=FALSE,alphaLevel=0.05) {
  .boot.ratio.test <- function(boot.cube, critical.value = 2) {
    boot.cube.mean <- apply(boot.cube, c(1, 2), mean)
    boot.cube.mean_repeat <- array(boot.cube.mean, dim = c(dim(boot.cube)))
    boot.cube.dev <- (boot.cube - boot.cube.mean_repeat)^2
    s.boot <- (apply(boot.cube.dev, c(1, 2), mean))^(1/2)
    boot.ratios <- boot.cube.mean/s.boot
    significant.boot.ratios <- (abs(boot.ratios) > critical.value)
    rownames(boot.ratios) <- rownames(boot.cube)
    rownames(significant.boot.ratios) <- rownames(boot.cube)
    return(list(sig.boot.ratios = significant.boot.ratios, 
                boot.ratios = boot.ratios))
  }
  nN = NROW(DATA1)
  if (nN != NROW(DATA2)) {
    stop("input matrices not conformable")
  }
  X <- apply(DATA1, 2, scale0, center = center1, scale = scale1)
  Y <- apply(DATA2, 2, scale0, center = center2, scale = scale2)
  nI = NCOL(X)
  nJ = NCOL(Y)
  maxRank <- min(nI, nJ)
  if (maxRank < nf2keep) 
    nf2keep = maxRank
  if (is.null(Fi) | is.null(Fj)) {
    S <- t(X) %*% Y
    svd.S <- svd(S, nu = nf2keep, nv = nf2keep)
    if (nf2keep > length(svd.S$d)) 
      nf2keep = length(svd.S$d)
    Lx <- X %*% svd.S$u
    Ly <- Y %*% svd.S$v
    Fi <- svd.S$u * matrix(svd.S$d, nI, nf2keep, byrow = TRUE)
    Fj <- svd.S$v * matrix(svd.S$d, nJ, nf2keep, byrow = TRUE)
  }
  else {
    nL = min(NCOL(Fi), NCOL(Fj))
    if (nL < nf2keep) 
      nf2keep = nL
    Fi = Fi[, 1:nf2keep]
    Fj = Fj[, 1:nf2keep]
    delta.inv <- 1/sqrt(colSums(Fi^2))
    Lx <- X %*% (Fi * matrix(delta.inv, nI, nf2keep, byrow = TRUE))
    Ly <- Y %*% (Fj * matrix(delta.inv, nJ, nf2keep, byrow = TRUE))
  }
  fj.boot <- array(NA, dim = c(nJ, nf2keep, nIter))
  dimnames(fj.boot)[1] <- list(colnames(Y))
  dimnames(fj.boot)[2] <- list(paste0("Dimension ", 1:nf2keep))
  dimnames(fj.boot)[3] <- list(paste0("Iteration ", 1:nIter))
  fi.boot <- array(NA, dim = c(nI, nf2keep, nIter))
  dimnames(fi.boot)[1] <- list(colnames(X))
  dimnames(fi.boot)[2] <- list(paste0("Dimension ", 1:nf2keep))
  dimnames(fi.boot)[3] <- list(paste0("Iteration ", 1:nIter))
  if (eig) {
    eigenValues <- matrix(0, nrow = nIter, ncol = maxRank)
    colnames(eigenValues) <- paste0("Dimension ", 1:maxRank)
    rownames(eigenValues) <- paste0("Iteration ", 1:nIter)
    fixedEigenvalues <- sv2(compS(X, center1 = center1, 
                                  scale1 = scale1, Y, center2 = center2, scale2 = scale2))
    names(fixedEigenvalues) <- paste0("Dimension ", 1:length(fixedEigenvalues))
  }
  for (ell in 1:nIter) {
    boot.index <- bootF(sublist,famlist)
    fi.boot[, , ell] <- t(X[boot.index, ]) %*% Ly[boot.index, 
    ]
    fj.boot[, , ell] <- t(Y[boot.index, ]) %*% Lx[boot.index, 
    ]
    if (eig) {
      eigenS <- sv2(compS(X[boot.index, ], center1 = center1, 
                          scale1 = scale1, Y[boot.index, ], center2 = center2, 
                          scale2 = scale2))
      index <- min(maxRank, length(eigenS))
      eigenValues[ell, 1:index] <- eigenS
    }
  }
  BR.j <- .boot.ratio.test(fj.boot, critical.value)
  BR.i <- .boot.ratio.test(fi.boot, critical.value)
  return.list <- structure(list(bootstrapBrick.i = fi.boot, 
                                bootRatios.i = BR.i$boot.ratios, bootRatiosSignificant.i = BR.i$sig.boot.ratios, 
                                bootstrapBrick.j = fj.boot, bootRatios.j = BR.j$boot.ratios, 
                                bootRatiosSignificant.j = BR.j$sig.boot.ratios), class = "bootBrick.ij4plsc")
  if (eig) {
    eigenValues <- eigenValues[, colSums(eigenValues) > 
                                 0]
    return.list$eigenValues = eigenValues
    sortedEigenValues <- apply(eigenValues, 2, sort)
    index = round(nIter * (alphaLevel/2))
    if (index == 0) 
      index <- 1
    eigenCI = sortedEigenValues[c(index, nIter - (index - 
                                                    1)), ]
    minCI <- as.character(alphaLevel/2)
    substr(minCI, 1, 2) <- "_"
    minCI <- paste0("LB", minCI)
    maxCI <- as.character(1 - (alphaLevel/2))
    substr(maxCI, 1, 2) <- "_"
    maxCI <- paste0("UB", maxCI)
    rownames(eigenCI) <- c(minCI, maxCI)
    return.list$fixedEigenvalues <- fixedEigenvalues
    return.list$eigenCI <- eigenCI
    class(return.list) <- "bootBrick.ij.eig4plsc"
  }
  return(return.list)
}

# Driver ------------------------------------------------------------------

#Catches arguments.
args <- commandArgs(trailingOnly=T)
k <- args[1]

#Set base path and PLSC path.
subgroup <- 'full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/')
plscpath <-  paste0(basepath,'allreconfig/PLSC_freq/')

#Set parameters.
nk <- as.numeric(k)
theme_set(theme_gray())
iternum <- '10000'
permver <- 'byColumns'
bootcv <- '2.5'

#Get cognitive and other predictor data.
infile <- paste0('../outputs/c_cognition/',subgroup,'/pred_all.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Isolate desired cognitive variables.
ourcog <- colnames(othercog)
ourcog <- c('Subject','PMAT24','PicVoc','ProcSpeed','CardSort','Flanker',
            'ListSort','PicSeq','ReadVoc','WordMem','LineOrient')
cogmat <- othercog[,ourcog]
ncog <- ncol(cogmat) - 1

#Get predictor data.
infile <- paste0(basepath,'avg_substat.csv')
dynmat <- read_csv(infile,col_names=T) %>% 
  rename_at(1,~'Subject')
colnames(dynmat) <- gsub('sptrans_','trans_',colnames(dynmat))
predmat <- dynmat 
colnames(predmat) <- gsub('-','t',colnames(predmat))

#Set labels.
dynlabs <- c('occur','dwell','numtrans','trans_')
ndyn <- length(dynlabs)

#If there are NA in the matrix.
if (any(is.na(predmat))) {
  
  #Save the subjects with errors.
  remsubs <- predmat[rowSums(is.na(predmat))!=0,'Subject']
  outfile <- paste0(plscpath,'suberr.txt')
  write.table(remsubs,file=outfile,sep='\n',row.names=F,col.names=F)
  
  #Change matrices.
  predmat <- predmat[rowSums(is.na(predmat))==0,]
  keepsubs <- predmat[,'Subject']
  keepidx <- as.character(unlist(othercog[,'Subject'],use.names=F)) %in% as.character(unlist(keepsubs,use.names=F))
  othercog <- othercog[keepidx,]
  keepidx <- as.character(unlist(cogmat[,'Subject'],use.names=F)) %in% as.character(unlist(keepsubs,use.names=F))
  cogmat <- cogmat[keepidx,]
}

#Prepare data.
data_X <- data.frame(predmat)
rownames(data_X) <- data_X[,'Subject']
data_X <- within(data_X,rm('Subject'))
data_Y <- data.frame(cogmat)
rownames(data_Y) <- data_Y[,'Subject']
data_Y <- within(data_Y,rm('Subject'))

#Extract labels.
ncog <- ncol(data_Y)
npred <- ncol(data_X)
coglabs <- colnames(data_Y)
predlabs <- colnames(data_X)

#Residualize confounds.
allmat <- cbind(data_Y,data_X,othercog)
xcont_list <- c('Age')
xcat_list <- c('Gender')
ycont_list <- c('Age')
ycat_list <- c('Gender')
coef_Y <- list()
coef_X <- list()
for (yidx in 1:ncog) {
  cY <- colnames(data_Y)[yidx]
  outlist <- resid_maker(cY,ycont_list,ycat_list,allmat,'raw')
  cY_resid <- outlist[[1]]
  cY_coef <- outlist[[2]]
  data_Y[,yidx] <- cY_resid
  coef_Y[[cY]] <- cY_coef
}
for (xidx in 1:npred) {
  cX <- colnames(data_X)[xidx]
  outlist <- resid_maker(cX,xcont_list,xcat_list,allmat,'raw')
  cX_resid <- outlist[[1]]
  cX_coef <- outlist[[2]]
  data_X[,xidx] <- cX_resid
  coef_X[[cX]] <- cX_coef
}
# Family ------------------------------------------------------------------

#Get family structure.
sublist <- rownames(data_Y)
infile <- '../inputs/data/hcp/RESTRICTED_HCP1200_Data_Dictionary.csv'
demomat <- read.csv(infile,row.names=1)
famlist <- demomat[sublist,'Family_ID']
dataF <- famlist

# PLSC --------------------------------------------------------------------

#Generate outpath.
outpath <- plscpath
dir.create(outpath,recursive=T)

#Run PLSC.
data1 <- data_X
data2 <- data_Y
pls.res <- tepPLS(data1,data2, 
                  center1=T,scale1=T, 
                  center2=T,scale2=T, 
                  graphs = FALSE)

#Calculate variance explained using the latent variables and save.
lvx_all <- pls.res$TExPosition.Data$lx
lvy_all <- pls.res$TExPosition.Data$ly
lvx_var <- sum(cor(lvx_all,data_X)^2)/(ncol(data_X)*ncol(lvx_all))
lvy_var <- sum(cor(lvy_all,data_Y)^2)/(ncol(data_Y)*ncol(lvy_all))
outmat <- data.frame(c(lvx_var,lvy_var),row.names=c('X','Y'))
colnames(outmat) <- c('Variance Explained')
outfile <- paste0(outpath,'main_var_explained.csv')
write.csv(outmat,outfile)

# Permutation Significance ------------------------------------------------

#Generate outpath.
outpath <- paste0(plscpath,'perm/')
dir.create(outpath,recursive=T)

#PLSC permutation model fit.
set.seed(12345)
perm.plsc <- perm4PLSC_F(data1,data2,dataF,
                         center1=T,scale1=T, 
                         center2=T,scale2=T, 
                         nIter=as.numeric(iternum),
                         permType=permver,
                         compact=F)

#Save singular value table.
singexp <- pls.res$TExPosition.Data$pdq$t
singval <- pls.res$TExPosition.Data$pdq$Dv
nsv <- min(ncol(data1),ncol(data2))
nsv.perm <- length(singval)
singcor <- matrix(NA,nrow=1,ncol=nsv)
for (i in 1:nsv.perm) {
  clx <- pls.res$TExPosition.Data$lx[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  singcor[,i] <- cor(clx,cly,method='spearman')
}
singp <- perm.plsc$pEigenvalues
singsig <- as.character(perm.plsc$pEigenvalues < 0.05)
svmat <- matrix(NA,nrow=5,ncol=nsv)
svmat[1,1:length(singexp)] <- singexp
svmat[2,1:length(singval)] <- singval
svmat[3,1:length(singcor)] <- singcor
svmat[4,1:length(singp)] <- singp
svmat[5,1:length(singsig)] <- singsig
rownames(svmat) <- c('Explained','SV','Correlation','P-value','Significance')
outfile <- paste0(outpath,'permsig.csv')
write.csv(svmat,outfile,row.names=T)

#Plot.
sigcolors <- rep('nosig',nsv)
sigcolors[perm.plsc$pEigenvalues < 0.05] <- 'sig'
lvplot <- data.frame(singexp)
colnames(lvplot) <- 'Explained'
lvplot['Significant'] <- sigcolors
lvplot['LV'] <- factor(paste0('LV',as.character(1:nsv)),
                       levels=paste0('LV',as.character(1:nsv)))
ggobj <- ggplot(lvplot,aes(x=LV,y=Explained,group=1)) +
  geom_line() +
  geom_point(aes(colour=factor(Significant),size=factor(Significant))) +
  scale_colour_manual(name='P-value',
                      values=c("black","red"),
                      labels=c('≥ 0.05','< 0.05')) +
  theme(axis.title.x=element_blank()) +
  ylab('Variance Explained') 
outfile <- paste0(outpath,'permsig.jpg')
ggsave(outfile,ggobj,units='in',width=5,height=2)

#Set significant LVs.
permsig <- (perm.plsc$pEigenvalues < 0.05)
svlist <- which(permsig)
nsv <- length(svlist)

#Generate outpath.
outpath <- paste0(plscpath,'observations/')
dir.create(outpath,recursive=T)

#Set colors.
norm.color <- '#696969'
ysig.color <- '#4E79A7'
xsig.color <- '#4E79A7'
pos.color <- '#FF0000'
neg.color <- '#4E79A7'

#Plot Y loadings.
for (lidx in svlist) {
  
  #Get significance label.
  if (permsig[lidx]) {
    addsig = '*'
  } else {
    addsig = ''
  }
  
  #Plot.
  llab <- as.character(lidx)
  ly <- pls.res$TExPosition.Data$fj[,lidx]
  if (permsig[lidx]) {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(ysig.color,ncog),
                           font.size=3) 
  } else {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(norm.color,ncog),
                           font.size=3)
  }
  
  #Save.
  outfile <- paste0(outpath,'baseLVY',llab,'.jpg')
  ggsave(outfile,plty,units='in',width=5,height=3)
}

#Generate flippers.
outfile <- paste0(outpath,'yflip.csv')
if (file.exists(outfile)) {
  yflip <- as.vector(t(read.csv(outfile,header=F)))
} else {
  
  #Set up flippers based on specific loadings.
  yflip <- rep(F,each=nsv)

  #k = 2.
  if (nk==2) {
    if (pls.res$TExPosition.Data$fj['ProcSpeed',1] < 0) {
      yflip[1] <- T
    } 
    if (pls.res$TExPosition.Data$fj['PicSeq',2] < 0) {
      yflip[2] <- T
    } 
  
  #k = 3-4.
  } else if (nk==3|nk==4) {
    if (mean(pls.res$TExPosition.Data$fj[,1]) < 0) {
      yflip[1] <- T
    } 
    if (pls.res$TExPosition.Data$fj['ProcSpeed',2] < 0) {
      yflip[2] <- T
    } 

  #k = 5-11.
  } else if (nk==5|nk==6|nk==7|nk==8|nk==9|nk==10|nk==11) {
    if (mean(pls.res$TExPosition.Data$fj[,1]) < 0) {
      yflip[1] <- T
    } 
    if (pls.res$TExPosition.Data$fj['ProcSpeed',2] < 0) {
      yflip[2] <- T
    } 
    if (pls.res$TExPosition.Data$fj['PicSeq',3] < 0) {
      yflip[3] <- T
    }
  
  #k = 12.
  } else if (nk == '12') {
    if (mean(pls.res$TExPosition.Data$fj[,1]) < 0) {
      yflip[1] <- T
    } 
    if (pls.res$TExPosition.Data$fj['PicSeq',2] < 0) {
      yflip[2] <- T
    } 
    if (pls.res$TExPosition.Data$fj['ProcSpeed',3] < 0) {
      yflip[3] <- T
    }
  
  #Default.
  } else {
    if (mean(pls.res$TExPosition.Data$fj[,1]) < 0) {
      yflip[1] <- T
    } 
  }
  
  #Generate file.
  outflip <- t(yflip) 
  outflip[outflip == TRUE] = 'T'
  outflip[outflip == FALSE] = 'F'
  write.table(outflip,outfile,sep=',',row.names=F,col.names=F)
}

#Generate outpath.
outpath <- paste0(plscpath,'perm/')
dir.create(outpath,recursive=T)

#Produce max values based on flipped versions.
clc_xlim <- c()
clc_ylim <- c()
for (i in 1:nsv) {
  
  #Extract correlation.
  clx <- pls.res$TExPosition.Data$lx[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  if (yflip[i]) {
    clx <- -clx
    cly <- -cly
  }
  clc_xlim <- c(clc_xlim,max(clx),min(clx))
  clc_ylim <- c(clc_ylim,max(cly),min(cly))
}
clc_xmax <- max(clc_xlim)
clc_xmin <- min(clc_xlim)
clc_ymax <- max(clc_ylim)
clc_ymin <- min(clc_ylim)

#Plot scatterplots and do bivariate normality tests based on flipped versions.
testlabs <- c('mardia_skew','mardia_kurtosis','hz','royston','dh','energy')
ntest <- length(testlabs)
pr_norm <- matrix(NA,nrow=nsv,ncol=ntest)
rownames(pr_norm) <- paste0('LV',as.character(1:nsv))
colnames(pr_norm) <- testlabs
clc_txtsize <- 11
for (i in 1:nsv) {
  
  #Extract correlation.
  clx <- pls.res$TExPosition.Data$lx[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  if (yflip[i]) {
    clx <- -clx
    cly <- -cly
  }
  ccor <- cor(clx,cly,method='spearman')
  clc <- data.frame('LVX'=clx,'LVY'=cly)
  clc_labs <- paste0(colnames(clc),as.character(i))
  
  #Plot.
  pltc <- ggplot(aes(LVX,LVY),data=clc) +
    geom_point(size=1) +
    geom_smooth(method=lm,se=F) +
    scale_x_continuous(limits=c(clc_xmin,clc_xmax)) +
    scale_y_continuous(limits=c(clc_ymin,clc_ymax)) +
    labs(x=clc_labs[1],y=clc_labs[2]) +
    theme(axis.text.x=element_text(size=clc_txtsize),
          axis.text.y=element_text(size=clc_txtsize))
  
  #Save.
  outfile <- paste0(outpath,'LV',as.character(i),'_corr.jpg')
  ggsave(outfile,pltc,units='in',width=4,height=3)
  
  #Set current bivariate relationship.
  c_XY <- clc
  
  #Generate bivariate QQ plot and save.
  outfile <- paste0(outpath,'LV',as.character(i),'_qq.jpg')
  jpeg(file=outfile)
  mvn(c_XY,multivariatePlot='qq')
  dev.off()
  
  #Do each test.
  for (ctest in testlabs) {

    #Whether indexing needs to be specific.
    if (grepl('mardia',ctest,fixed=T)) {
      
      #Skew or kurtosis value.
      if (ctest=='mardia_skew') {
        ctest_raw <- 'mardia'
        ctest_lab <- 'Mardia Skewness'
      } else if (ctest=='mardia_kurtosis') {
        ctest_raw <- 'mardia'
        ctest_lab <- 'Mardia Kurtosis'
      }
      
      #Do test.
      fitted <- mvn(c_XY,mvnTest=ctest_raw)$multivariateNormality
      fac <- fitted[fitted$Test==ctest_lab,'p value']
      cpr <- as.numeric(levels(fac)[fac])
    } else {
      
      #Do test.
      fitted <- mvn(c_XY,mvnTest=ctest)$multivariateNormality
      cpr <- fitted[1,'p value']
    }
    
    #Save it.
    pr_norm[paste0('LV',as.character(i)),ctest] <- cpr
  }
}

#Threshold p-values and save.
pr_norm_thres <- pr_norm
pr_norm_thres[pr_norm < 0.05] <- NA
pr_norm_thres <- signif(pr_norm_thres,2)
outfile <- paste0(outpath,'binormtest_p.csv')
write.csv(data.frame(pr_norm),outfile)
outfile <- paste0(outpath,'binormtest_pthres.csv')
write.csv(data.frame(pr_norm_thres),outfile)

#Save scores for flipped versions.
lx_scores <- matrix(NA,nrow=nrow(data_X),ncol=nsv)
rownames(lx_scores) <- rownames(data_X)
ly_scores <- lx_scores
for (i in 1:nsv) {
  
  #Extract scores.
  clx <- pls.res$TExPosition.Data$lx[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  if (yflip[i]) {
    clx <- -clx
    cly <- -cly
  }
  lx_scores[,i] <- clx
  ly_scores[,i] <- cly
}
colnames(lx_scores) <- paste0('LVX',as.character(1:nsv))
colnames(ly_scores) <- paste0('LVY',as.character(1:nsv))
outfile <- paste0(outpath,'LVX_scores.csv')
write.csv(lx_scores,outfile)
outfile <- paste0(outpath,'LVY_scores.csv')
write.csv(ly_scores,outfile)

# Bootstrap Stability -----------------------------------------------------

#Generate outpath.
outpath <- paste0(plscpath,'boot/')
dir.create(outpath,recursive=T)

#PLSC bootstrap feature stability.
bootres.plsc <- Boot4PLSC_F(data1,data2,
                            sublist,famlist,
                            eig=T,
                            nIter=as.numeric(iternum),
                            critical.value=as.numeric(bootcv),
                            nf2keep=nsv)

#Generate color for significant bootstrap ratios.
bootcolor.X <- ifelse(bootres.plsc[["bootRatiosSignificant.i"]],xsig.color,norm.color)
bootcolor.Y <- ifelse(bootres.plsc[["bootRatiosSignificant.j"]],ysig.color,norm.color)

#Plot X and Y loadings.
scaler <- 1
plot_list <- list()
pidx <- 1
lxtab <- matrix(NA,nrow=nsv,ncol=npred)
colnames(lxtab) <- predlabs
lytab <- matrix(NA,nrow=nsv,ncol=ncog)
colnames(lytab) <- coglabs
lxtab_full <- lxtab
lytab_full <- lytab
for (lidx in 1:nsv) {
  
  #Get significance label.
  if (permsig[lidx]) {
    addsig = '*'
  } else {
    addsig = ''
  }
  
  #Get and flip target Y.
  ly <- pls.res$TExPosition.Data$fj[,lidx]
  if (yflip[lidx]) {
    ly <- -ly
  }
  
  #Scale.
  ly <- ly*scaler
  
  #Plot.
  llab <- as.character(lidx)
  if (permsig[lidx]) {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=bootcolor.Y[,lidx],
                           font.size=3) 
  } else {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(norm.color,ncog),
                           font.size=3)
  }
  
  #Save.
  plot_list[[pidx]] <- plty
  pidx <- pidx + 1
  
  #Plot second version.
  collab <- ly
  collab[ly > 0] <- pos.color
  collab[ly < 0] <- neg.color
  collab[!(bootres.plsc[["bootRatiosSignificant.j"]][,lidx])] <- norm.color
  if (permsig[lidx]) {
    plty2 <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=collab,
                           font.size=3) 
  } else {
    plty2 <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(norm.color,ncog),
                           font.size=3)
  }
  outfile <- paste0(outpath,'LVY',llab,'_load.jpg')
  ggsave(outfile,plty2,units='in',width=2,height=3)
  
  #Get and flip target X.
  lx <- pls.res$TExPosition.Data$fi[,lidx] 
  if (yflip[lidx]) {
    lx <- -lx
  }
  
  #Scale.
  lx <- lx*scaler
  
  #Plot.
  if (permsig[lidx]) {
    pltx <- PrettyBarPlot2(lx,
                           threshold=0, 
                           color4bar=bootcolor.X[,lidx],
                           font.size=3)
  } else {
    pltx <- PrettyBarPlot2(lx,
                           threshold=0, 
                           color4bar=rep(norm.color,npred),
                           font.size=3)
  }  
  
  #Save.
  plot_list[[pidx]] <- pltx
  pidx <- pidx + 1
  
  #Save the loading tables.
  lx_thres <- lx
  lxtab_full[lidx,] <- lx_thres
  lx_thres[!bootres.plsc[["bootRatiosSignificant.i"]][,lidx]] <- NA
  lxtab[lidx,] <- lx_thres
  ly_thres <- ly
  lytab_full[lidx,] <- ly_thres
  ly_thres[!bootres.plsc[["bootRatiosSignificant.j"]][,lidx]] <- NA
  lytab[lidx,] <- ly_thres
  
  #Plot second version.
  for (didx in 1:ndyn) {
    
    #Extract.
    cdyn <- dynlabs[didx]
    dlx <- data.frame(t(lx)) %>%
      select(matches(cdyn))
    if (grepl('_',cdyn)) {
      namedyn <- gsub('_','',cdyn)
    } else {
      namedyn <- cdyn
    }
    dsig <- bootres.plsc[["bootRatiosSignificant.i"]][colnames(dlx),lidx]
    
    #Set up labels for plotting.
    if (namedyn == 'occur') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'dwell') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'numtrans') {
      inlabs <- paste0('Transition Number')
      nwidth <- 1
      nheight <- 3
    } else if (namedyn == 'trans') {
      inlabs <- gsub('trans_','',colnames(dlx))
      inlabs <- gsub('t','-',colnames(dlx))
      nwidth <- 10
      nheight <- 3
    } else if (namedyn == 'sim') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'jump') {
      inlabs <- gsub('jump_','',colnames(dlx))
      inlabs <- gsub('t','-',colnames(dlx))
      nwidth <- 10
      nheight <- 3
    }
    
    #Plot.
    collab <- dlx
    collab[dlx > 0] <- pos.color
    collab[dlx < 0] <- neg.color
    collab[!(dsig)] <- norm.color
    dlx <- unlist(as.vector(dlx))
    names(dlx) <- inlabs
    if (permsig[lidx]) {
      pltdx <- PrettyBarPlot2(dlx,
                              threshold=0, 
                              color4bar=collab,
                              font.size=3) 
    } else {
      pltdx <- PrettyBarPlot2(dlx,
                              threshold=0, 
                              color4bar=rep(norm.color,length(dlx)),
                              font.size=3)
    }
    outfile <- paste0(outpath,'LVX',llab,'_load_',namedyn,'.jpg')
    ggsave(outfile,pltdx,units='in',width=nwidth,height=nheight)
  }
}

#Format plots.
p1 <- plot_list[c(1,3,5)]
p2 <- plot_list[c(2,4,6)]

#Combine plots.
outfile <- paste0(outpath,'LVYall_load.jpg')
ggsave(outfile,arrangeGrob(grobs=p1,ncol=1),
       units='in',width=2,height=8)
outfile <- paste0(outpath,'LVXall_load.jpg')
ggsave(outfile,arrangeGrob(grobs=p2,ncol=1),
       units='in',width=30,height=8)

#Save tables.
outfile <- paste0(outpath,'LVXall_load.csv')
write.csv(lxtab,outfile)
outfile <- paste0(outpath,'LVYall_load.csv')
write.csv(lytab,outfile)
outfile <- paste0(outpath,'LVXall_load_full.csv')
write.csv(t(lxtab_full),outfile)
outfile <- paste0(outpath,'LVYall_load_full.csv')
write.csv(t(lytab_full),outfile)

#Plot X and Y bootstrap ratios.
scaler <- 1
plot_list <- list()
pidx <- 1
lxtab <- matrix(NA,nrow=nsv,ncol=npred)
colnames(lxtab) <- predlabs
lytab <- matrix(NA,nrow=nsv,ncol=ncog)
colnames(lytab) <- coglabs
lxtab_full <- lxtab
lytab_full <- lytab
for (lidx in 1:nsv) {
  
  #Get significance label.
  if (permsig[lidx]) {
    addsig = '*'
  } else {
    addsig = ''
  }
  
  #Get and flip target Y.
  ly <- bootres.plsc$bootRatios.j[,lidx]
  if (yflip[lidx]) {
    ly <- -ly
  }
  
  #Scale.
  ly <- ly*scaler
  
  #Plot.
  llab <- as.character(lidx)
  if (permsig[lidx]) {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=bootcolor.Y[,lidx],
                           font.size=3) 
  } else {
    plty <- PrettyBarPlot2(ly,
                           threshold=0, 
                           color4bar=rep(norm.color,ncog),
                           font.size=3)
  }
  
  #Save.
  plot_list[[pidx]] <- plty
  pidx <- pidx + 1
  
  #Plot second version.
  collab <- ly
  collab[ly > 0] <- pos.color
  collab[ly < 0] <- neg.color
  collab[!(bootres.plsc[["bootRatiosSignificant.j"]][,lidx])] <- norm.color
  if (permsig[lidx]) {
    plty2 <- PrettyBarPlot2(ly,
                            threshold=0, 
                            color4bar=collab,
                            font.size=3) 
  } else {
    plty2 <- PrettyBarPlot2(ly,
                            threshold=0, 
                            color4bar=rep(norm.color,ncog),
                            font.size=3)
  }
  outfile <- paste0(outpath,'LVY',llab,'_bootrat.jpg')
  ggsave(outfile,plty2,units='in',width=2,height=3)
  
  #Get and flip target X.
  lx <- bootres.plsc$bootRatios.i[,lidx]
  if (yflip[lidx]) {
    lx <- -lx
  }
  
  #Scale.
  lx <- lx*scaler
  
  #Plot.
  if (permsig[lidx]) {
    pltx <- PrettyBarPlot2(lx,
                           threshold=0, 
                           color4bar=bootcolor.X[,lidx],
                           font.size=3)
  } else {
    pltx <- PrettyBarPlot2(lx,
                           threshold=0, 
                           color4bar=rep(norm.color,npred),
                           font.size=3)
  }  
  
  #Save.
  plot_list[[pidx]] <- pltx
  pidx <- pidx + 1
  
  #Save the loading tables.
  lx_thres <- lx
  lxtab_full[lidx,] <- lx_thres
  lx_thres[!bootres.plsc[["bootRatiosSignificant.i"]][,lidx]] <- NA
  lxtab[lidx,] <- lx_thres
  ly_thres <- ly
  lytab_full[lidx,] <- ly_thres
  ly_thres[!bootres.plsc[["bootRatiosSignificant.j"]][,lidx]] <- NA
  lytab[lidx,] <- ly_thres
  
  #Plot second version.
  for (didx in 1:ndyn) {
    
    #Extract.
    cdyn <- dynlabs[didx]
    dlx <- data.frame(t(lx)) %>%
      select(matches(cdyn))
    if (grepl('_',cdyn)) {
      namedyn <- gsub('_','',cdyn)
    } else {
      namedyn <- cdyn
    }
    dsig <- bootres.plsc[["bootRatiosSignificant.i"]][colnames(dlx),lidx]
    
    #Set up labels for plotting.
    if (namedyn == 'occur') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'dwell') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'numtrans') {
      inlabs <- paste0('Transition Number')
      nwidth <- 1
      nheight <- 3
    } else if (namedyn == 'trans') {
      inlabs <- gsub('trans_','',colnames(dlx))
      inlabs <- gsub('t','-',colnames(dlx))
      nwidth <- 10
      nheight <- 3
    } else if (namedyn == 'sim') {
      inlabs <- paste0('State ',1:nk)
      nwidth <- 1.5
      nheight <- 3
    } else if (namedyn == 'jump') {
      inlabs <- gsub('jump_','',colnames(dlx))
      inlabs <- gsub('t','-',colnames(dlx))
      nwidth <- 10
      nheight <- 3
    }
    
    #Plot.
    collab <- dlx
    collab[dlx > 0] <- pos.color
    collab[dlx < 0] <- neg.color
    collab[!(dsig)] <- norm.color
    dlx <- unlist(as.vector(dlx))
    names(dlx) <- inlabs
    if (permsig[lidx]) {
      pltdx <- PrettyBarPlot2(dlx,
                              threshold=0, 
                              color4bar=collab,
                              font.size=3) 
    } else {
      pltdx <- PrettyBarPlot2(dlx,
                              threshold=0, 
                              color4bar=rep(norm.color,length(dlx)),
                              font.size=3)
    }
    outfile <- paste0(outpath,'LVX',llab,'_bootrat_',namedyn,'.jpg')
    ggsave(outfile,pltdx,units='in',width=nwidth,height=nheight)
  }
}

#Format plots.
p1 <- plot_list[c(1,3,5)]
p2 <- plot_list[c(2,4,6)]

#Combine plots.
outfile <- paste0(outpath,'LVYall_bootrat.jpg')
ggsave(outfile,arrangeGrob(grobs=p1,ncol=1),
       units='in',width=2,height=8)
outfile <- paste0(outpath,'LVXall_bootrat.jpg')
ggsave(outfile,arrangeGrob(grobs=p2,ncol=1),
       units='in',width=30,height=8)

#Save tables.
outfile <- paste0(outpath,'LVXall_bootrat.csv')
write.csv(lxtab,outfile)
outfile <- paste0(outpath,'LVYall_bootrat.csv')
write.csv(lytab,outfile)
outfile <- paste0(outpath,'LVXall_bootrat_full.csv')
write.csv(t(lxtab_full),outfile)
outfile <- paste0(outpath,'LVYall_bootrat_full.csv')
write.csv(t(lytab_full),outfile)

# Reproducibility ---------------------------------------------------------

#Create replication outpath.
outpath <- paste0(plscpath,'replicate/')
dir.create(outpath,recursive=T)

#Create a matrix for the SV, LVX, and LVY replicability values. Collect variance
#explained.
niter <- as.numeric(iternum)
maxsv <- min(ncol(data_X),ncol(data_Y))
svrep <- matrix(NA,nrow=niter,ncol=maxsv)
lvxrep <- matrix(NA,nrow=niter,ncol=maxsv)
lvyrep <- matrix(NA,nrow=niter,ncol=maxsv)
pls.var.x <- c()
pls.var.y <- c()

#For each of the splits.
set.seed(12345)
for (spidx in 1:niter) {
  
  #Extract splits.
  nfold <- 2
  both_idx <- cvF(sublist,famlist,nfold)
  idx_1 <- rownames(both_idx)[both_idx==1]
  idx_2 <- rownames(both_idx)[both_idx==2]
  x1 <- allmat[idx_1,predlabs]
  y1 <- allmat[idx_1,coglabs]
  a1 <- allmat[idx_1,]
  x2 <- allmat[idx_2,predlabs]
  y2 <- allmat[idx_2,coglabs]
  a2 <- allmat[idx_2,]
  
  #Residualize first set and extract coefficients, mean, and SD.
  coef_Y1 <- list()
  st_Y1 <- matrix(NA,2,ncog)
  rownames(st_Y1) <- c('mean','sd')
  colnames(st_Y1) <- coglabs
  coef_X1 <- list()
  st_X1 <- matrix(NA,2,npred)
  rownames(st_X1) <- c('mean','sd')
  colnames(st_X1) <- predlabs
  for (yidx in 1:ncog) {
    cY <- colnames(y1)[yidx]
    outlist <- resid_maker(cY,ycont_list,ycat_list,a1,'raw')
    cY_resid <- outlist[[1]]
    cY_coef <- outlist[[2]]
    y1[,yidx] <- cY_resid
    coef_Y1[[cY]] <- cY_coef
    st_Y1['mean',cY] <- mean(cY_resid[,cY])
    st_Y1['sd',cY] <- sd(cY_resid[,cY])
  }
  for (xidx in 1:npred) {
    cX <- colnames(x1)[xidx]
    outlist <- resid_maker(cX,xcont_list,xcat_list,a1,'raw')
    cX_resid <- outlist[[1]]
    cX_coef <- outlist[[2]]
    x1[,xidx] <- cX_resid
    coef_X1[[cX]] <- cX_coef
    st_X1['mean',cX] <- mean(cX_resid[,cX])
    st_X1['sd',cX] <- sd(cX_resid[,cX])
  }
  
  #Residualize and standardize second set using first set coefficients, mean, and SD.
  for (yidx in 1:ncog) {
    cY <- colnames(y2)[yidx]
    cY_coef <- coef_Y1[[cY]]
    outlist <- resid_maker(cY,ycont_list,ycat_list,a2,'raw',cY_coef)
    cY_resid <- outlist[[1]]
    y2[,yidx] <- cY_resid
  }
  for (xidx in 1:npred) {
    cX <- colnames(x2)[xidx]
    cX_coef <- coef_X1[[cX]]
    outlist <- resid_maker(cX,xcont_list,xcat_list,a2,'raw',cX_coef)
    cX_resid <- outlist[[1]]
    x2[,xidx] <- cX_resid
  }
  y2 <- sweep(y2,2,st_Y1['mean',],'-')
  y2 <- sweep(y2,2,st_Y1['sd',],'/')
  x2 <- sweep(x2,2,st_X1['mean',],'-')
  x2 <- sweep(x2,2,st_X1['sd',],'/')
  
  #Do PLS.
  pls1 <- tepPLS(x1,y1,
                 center1=T,scale1=T,
                 center2=T,scale2=T,
                 graphs=F)
  pls2 <- tepPLS(x2,y2,
                 center1=F,scale1=F,
                 center2=F,scale2=F,
                 graphs=F)
  
  #Calculate second correlation matrix.
  rxy2 <- cor(x2,y2)
  
  #Get first and second u and v.
  u1 <- pls1$TExPosition.Data$pdq$p
  v1 <- pls1$TExPosition.Data$pdq$q
  u2 <- pls2$TExPosition.Data$pdq$p
  v2 <- pls2$TExPosition.Data$pdq$q
  
  #Calculate SV reproducibility value using the second set Rxy and first set u and v.
  svrep_val <- t(u1)%*%rxy2%*%v1
  svrep_val <- diag(svrep_val)
  svrep[spidx,] <- svrep_val
  
  #Calculate LVX and LVY reproducibility value using the first and second set u and v.
  lvxrep_val <- t(u1)%*%u2
  lvxrep_val <- abs(diag(lvxrep_val))
  lvxrep[spidx,] <- lvxrep_val
  lvyrep_val <- t(v1)%*%v2
  lvyrep_val <- abs(diag(lvyrep_val))
  lvyrep[spidx,] <- lvyrep_val
  
  #Flip the sign if the correlation to the original singular vector is negative for every dimension.
  #Use just the Y singular vectors to initially check as X singular vectors will need the same flip.
  flip <- diag(cor(pls.res$TExPosition.Data$pdq$q,pls1$TExPosition.Data$pdq$q)) < 0
  pls1$TExPosition.Data$pdq$p[,flip] <- pls1$TExPosition.Data$pdq$p[,flip]*-1
  pls1$TExPosition.Data$pdq$q[,flip] <- pls1$TExPosition.Data$pdq$q[,flip]*-1
  
  #Project to test set which has been scaled using the train set by multiplying the test set 
  #X and Y by the weights represented by the training set X and Y singular vectors respectively,
  #generating estimated test set latent variables.
  train.res <- pls1$TExPosition.Data 
  lx2 <- as.matrix(x2) %*% train.res$pdq$p
  ly2 <- as.matrix(y2) %*% train.res$pdq$q
  
  #Calculate variance explained using the test set latent variables and save.
  lvx_var <- sum(cor(lx2,x2)^2)/(ncol(x2)*ncol(lx2))
  lvy_var <- sum(cor(ly2,y2)^2)/(ncol(y2)*ncol(ly2))
  pls.var.x <- c(pls.var.x,lvx_var)
  pls.var.y <- c(pls.var.y,lvy_var)
}

#Find the mean and std columnwise and z-score.
svrep_mean <- colMeans(svrep)
svrep_std <- apply(svrep,2,sd)
svrep_z <- svrep_mean / svrep_std
lvxrep_mean <- colMeans(lvxrep)
lvxrep_std <- apply(lvxrep,2,sd)
lvxrep_z <- lvxrep_mean / lvxrep_std
lvyrep_mean <- colMeans(lvyrep)
lvyrep_std <- apply(lvyrep,2,sd)
lvyrep_z <- lvyrep_mean / lvyrep_std

#Save each file.
outlabs <- c('svrep','svrep_z','lvxrep','lvxrep_z','lvyrep','lvyrep_z')
nout <- length(outlabs)
outvars <- list(svrep,svrep_z,lvxrep,lvxrep_z,lvyrep,lvyrep_z)
for (oidx in 1:nout) {
  
  #Extract.
  olab <- outlabs[[oidx]]
  ovar <- outvars[[oidx]]
  if (is.matrix(ovar)) {
    colnames(ovar) <- paste0('LV',1:maxsv)
  } else {
    ovar <- t(data.frame(ovar))
    colnames(ovar) <- paste0('LV',1:maxsv)
  }
  
  #Save.
  outfile <- paste0(outpath,olab,'.csv')
  write.csv(ovar,outfile,row.names=F)
}

#Average variance explained across splits and save.
lvx_var <- mean(pls.var.x)
lvy_var <- mean(pls.var.y)
outmat <- data.frame(c(lvx_var,lvy_var),row.names=c('X','Y'))
colnames(outmat) <- c('Variance Explained')
outfile <- paste0(outpath,'replicate_var_explained.csv')
write.csv(outmat,outfile)

#Create replication plot outpath.
plotpath <- paste0(plscpath,'replicate/plots/')
dir.create(plotpath,recursive=T)

#Label LVs.
lvlabs <- paste0('LV',1:maxsv)

#Plot SV replicability.
if (length(unique(svrep_z > 1.95)) == 2) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
} else if (unique(svrep_z > 1.95) == F) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('≤ 1.95')))
} else if (unique(svrep_z > 1.95) == T) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('> 1.95')))
}
svrep_color['LV'] <- lvlabs
plotmat <- data.frame(svrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(svrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
svrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
  geom_boxplot(outlier.size=0.1) +
  ggtitle('SV Reproducibility') +
  ylab('Reproducibility') +
  theme(axis.title.x = element_blank())
outfile <- paste0(plotpath,'svrep.jpg')
ggsave(outfile,svrep_plot,units='in')

#Plot LVX replicability.
if (length(unique(lvxrep_z > 1.95)) == 2) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
} else if (unique(svrep_z > 1.95) == F) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('≤ 1.95')))
} else if (unique(svrep_z > 1.95) == T) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('> 1.95')))
}
lvxrep_color['LV'] <- lvlabs
plotmat <- data.frame(lvxrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(lvxrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
lvxrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
  geom_boxplot(outlier.size=0.1) +
  ggtitle('LVX Reproducibility') +
  ylab('Reproducibility') +
  theme(axis.title.x = element_blank())
outfile <- paste0(plotpath,'lvxrep.jpg')
ggsave(outfile,lvxrep_plot,units='in')

#Plot LVY replicability.
if (length(unique(lvyrep_z > 1.95)) == 2) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
} else if (unique(svrep_z > 1.95) == F) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('≤ 1.95')))
} else if (unique(svrep_z > 1.95) == T) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('> 1.95')))
}
lvyrep_color['LV'] <- lvlabs
plotmat <- data.frame(lvyrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(lvyrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
lvyrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
  geom_boxplot(outlier.size=0.1) +
  ggtitle('LVY Reproducibility') +
  ylab('Reproducibility') +
  theme(axis.title.x = element_blank())
outfile <- paste0(plotpath,'lvyrep.jpg')
ggsave(outfile,lvyrep_plot,units='in')

#Isolate the significant LVs only.
siglv <- svlist
lvlabs <- paste0('LV',as.character(siglv))
svrep <- svrep[,siglv]
svrep_z <- svrep_z[siglv]
lvxrep <- lvxrep[,siglv]
lvxrep_z <- lvxrep_z[siglv]
lvyrep <- lvyrep[,siglv]
lvyrep_z <- lvyrep_z[siglv]
nsv <- length(siglv)

#Set summary text size.
sumtxt <- 10

#SV.
if (length(unique(svrep_z > 1.95)) == 2) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
  ncolors <- 'ft'
} else if (unique(svrep_z > 1.95) == F) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('≤ 1.95')))
  ncolors <- 'fone'
} else if (unique(svrep_z > 1.95) == T) {
  svrep_color <- data.frame(factor((svrep_z > 1.95),labels=c('> 1.95')))
  ncolors <- 'tone'
}
svrep_color['LV'] <- lvlabs
plotmat <- data.frame(svrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(svrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
if (ncolors == 'ft') {
  svrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    ggtitle('Singular Values') +
    ylab('Reproducibility') +
    theme(axis.title.x = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'fone') {
  svrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#F8766D')) +
    ggtitle('Singular Values') +
    ylab('Reproducibility') +
    theme(axis.title.x = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'tone') {
  svrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#00BFC4')) +
    ggtitle('Singular Values') +
    ylab('Reproducibility') +
    theme(axis.title.x = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
}

#LVX.
if (length(unique(lvxrep_z > 1.95)) == 2) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
  ncolors <- 'ft'
} else if (unique(lvxrep_z > 1.95) == F) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('≤ 1.95')))
  ncolors <- 'fone'
} else if (unique(lvxrep_z > 1.95) == T) {
  lvxrep_color <- data.frame(factor((lvxrep_z > 1.95),labels=c('> 1.95')))
  ncolors <- 'tone'
}
lvxrep_color['LV'] <- lvlabs
plotmat <- data.frame(lvxrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(lvxrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
if (ncolors == 'ft') {
  lvxrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    ggtitle('X Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'fone') {
  lvxrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#F8766D')) +
    ggtitle('X Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'tone') {
  lvxrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#00BFC4')) +
    ggtitle('X Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
}

#LVY.
if (length(unique(lvyrep_z > 1.95)) == 2) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('≤ 1.95', '> 1.95')))
  ncolors <- 'ft'
} else if (unique(lvyrep_z > 1.95) == F) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('≤ 1.95')))
  ncolors <- 'fone'
} else if (unique(lvyrep_z > 1.95) == T) {
  lvyrep_color <- data.frame(factor((lvyrep_z > 1.95),labels=c('> 1.95')))
  ncolors <- 'tone'
}
lvyrep_color['LV'] <- lvlabs
plotmat <- data.frame(lvyrep)
colnames(plotmat) <- lvlabs
plotmat <- stack(plotmat)
colnames(plotmat) <- c('Reproducibility','LV')
plotmat <- plotmat %>%
  full_join(lvyrep_color)
plotmat['LV'] <- factor(plotmat$LV,levels=lvlabs)
colnames(plotmat) <- c('Reproducibility','LV','Critical')
if (ncolors == 'ft') {
  lvyrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    ggtitle('Y Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'fone') {
  lvyrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#F8766D')) +
    ggtitle('Y Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
} else if (ncolors == 'tone') {
  lvyrep_plot <- ggplot(plotmat,aes(x=LV,y=Reproducibility,fill=Critical)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c('#00BFC4')) +
    ggtitle('Y Singular Vectors') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          axis.text.x=element_text(size=sumtxt),
          axis.text.y=element_text(size=sumtxt))
}  
(svrep_plot + lvyrep_plot + lvxrep_plot)
ggsave(paste0(plotpath,'repsummary.jpg'),units='in',width=5.25,height=3)
dev.off()

# K-fold CV ---------------------------------------------------------------

#Create cross-validation outpath.
outpath <- paste0(plscpath,'crossval/')
dir.create(outpath,recursive=T)

#Investigate results across for k-folds.
nfold <- 10
seedval <- 12345
set.seed(seedval)

#Collect fold outputs.
pls.kfcv <- list(lx.hat = matrix(NA, 
                                 nrow = nrow(pls.res$TExPosition.Data$lx),
                                 ncol = ncol(pls.res$TExPosition.Data$lx),
                                 dimnames = dimnames(pls.res$TExPosition.Data$lx)),
                 ly.hat = matrix(NA, 
                                 nrow = nrow(pls.res$TExPosition.Data$ly),
                                 ncol = ncol(pls.res$TExPosition.Data$ly),
                                 dimnames = dimnames(pls.res$TExPosition.Data$ly)),
                 lxly.cor = list(),
                 p = list(),
                 q = list(),
                 ci = list(),
                 cj = list(),
                 Dv = list(),
                 eig = list())
pls.var.x <- c()
pls.var.y <- c()

#Extract test fold IDs.
both_idx <- cvF(sublist,famlist,nfold)

#For each fold.
for (foidx in 1:nfold) {
  
  #Extract splits.
  idx_1 <- rownames(both_idx)[both_idx!=foidx]
  idx_2 <- rownames(both_idx)[both_idx==foidx]
  x1 <- allmat[idx_1,predlabs]
  y1 <- allmat[idx_1,coglabs]
  a1 <- allmat[idx_1,]
  x2 <- allmat[idx_2,predlabs]
  y2 <- allmat[idx_2,coglabs]
  a2 <- allmat[idx_2,]
  cat('Train:',length(idx_1),'\n')
  cat('Test:',length(idx_2),'\n')
  
  #Residualize first set and extract coefficients, mean, and SD.
  coef_Y1 <- list()
  st_Y1 <- matrix(NA,2,ncog)
  rownames(st_Y1) <- c('mean','sd')
  colnames(st_Y1) <- coglabs
  coef_X1 <- list()
  st_X1 <- matrix(NA,2,npred)
  rownames(st_X1) <- c('mean','sd')
  colnames(st_X1) <- predlabs
  for (yidx in 1:ncog) {
    cY <- colnames(y1)[yidx]
    outlist <- resid_maker(cY,ycont_list,ycat_list,a1,'raw')
    cY_resid <- outlist[[1]]
    cY_coef <- outlist[[2]]
    y1[,yidx] <- cY_resid
    coef_Y1[[cY]] <- cY_coef
    st_Y1['mean',cY] <- mean(cY_resid[,cY])
    st_Y1['sd',cY] <- sd(cY_resid[,cY])
  }
  for (xidx in 1:npred) {
    cX <- colnames(x1)[xidx]
    outlist <- resid_maker(cX,xcont_list,xcat_list,a1,'raw')
    cX_resid <- outlist[[1]]
    cX_coef <- outlist[[2]]
    x1[,xidx] <- cX_resid
    coef_X1[[cX]] <- cX_coef
    st_X1['mean',cX] <- mean(cX_resid[,cX])
    st_X1['sd',cX] <- sd(cX_resid[,cX])
  }
  
  #Residualize and standardize second set using first set coefficients, mean, and SD.
  for (yidx in 1:ncog) {
    cY <- colnames(y2)[yidx]
    cY_coef <- coef_Y1[[cY]]
    outlist <- resid_maker(cY,ycont_list,ycat_list,a2,'raw',cY_coef)
    cY_resid <- outlist[[1]]
    y2[,yidx] <- cY_resid
  }
  for (xidx in 1:npred) {
    cX <- colnames(x2)[xidx]
    cX_coef <- coef_X1[[cX]]
    outlist <- resid_maker(cX,xcont_list,xcat_list,a2,'raw',cX_coef)
    cX_resid <- outlist[[1]]
    x2[,xidx] <- cX_resid
  }
  y2 <- sweep(y2,2,st_Y1['mean',],'-')
  y2 <- sweep(y2,2,st_Y1['sd',],'/')
  x2 <- sweep(x2,2,st_X1['mean',],'-')
  x2 <- sweep(x2,2,st_X1['sd',],'/')
  
  #Run PLSC with train set and scale with its own parameters.
  train.pls <- tepPLS(x1,y1,
                      center1=T,scale1=T,
                      center2=T,scale2=T,
                      graphs=F)
  
  #Flip the sign if the correlation to the original singular vector is negative for every dimension.
  #Use just the Y singular vectors to initially check as X singular vectors will need the same flip.
  flip <- diag(cor(pls.res$TExPosition.Data$pdq$q,train.pls$TExPosition.Data$pdq$q)) < 0
  train.pls$TExPosition.Data$pdq$p[,flip] <- train.pls$TExPosition.Data$pdq$p[,flip]*-1
  train.pls$TExPosition.Data$pdq$q[,flip] <- train.pls$TExPosition.Data$pdq$q[,flip]*-1
  
  #Project to test set which has been scaled using the train set by multiplying the test set 
  #X and Y by the weights represented by the training set X and Y singular vectors respectively,
  #generating estimated test set latent variables.
  train.res <- train.pls$TExPosition.Data 
  lx2 <- as.matrix(x2) %*% train.res$pdq$p
  ly2 <- as.matrix(y2) %*% train.res$pdq$q
  
  #Calculate variance explained using the test set latent variables and save.
  lvx_var <- sum(cor(lx2,x2)^2)/(ncol(x2)*ncol(lx2))
  lvy_var <- sum(cor(ly2,y2)^2)/(ncol(y2)*ncol(ly2))
  pls.var.x <- c(pls.var.x,lvx_var)
  pls.var.y <- c(pls.var.y,lvy_var)
  
  #Save results from train set.
  pls.kfcv$p[[foidx]] <- train.pls$TExPosition.Data$pdq$p
  pls.kfcv$q[[foidx]] <- train.pls$TExPosition.Data$pdq$q
  pls.kfcv$ci[[foidx]] <- train.pls$TExPosition.Data$pdq$p^2
  pls.kfcv$cj[[foidx]] <- train.pls$TExPosition.Data$pdq$q^2
  pls.kfcv$Dv[[foidx]] <- train.pls$TExPosition.Data$pdq$Dv
  pls.kfcv$eig[[foidx]] <- train.pls$TExPosition.Data$pdq$eigs
  
  #Save results from test set, placing values for the test individuals back into 
  #the full matrix by ID.
  pls.kfcv$lx.hat[rownames(x2),] <- lx2
  pls.kfcv$ly.hat[rownames(y2),] <- ly2
}

#Collect outputs for each dimension based on the test values. Calculate the 
#correlation between the estimated test set latent variables of X and Y in order
#to generate the overall test set estimated relationship. Calculate the correlation
#between the original latent variables of X and Y and the estimated test set latent
#variables of X and Y in order to see how consistent they are.
pls.kfcv$lxly.cor[['cor.lxhatlyhat']] <- diag(cor(pls.kfcv$lx.hat,pls.kfcv$ly.hat,method='spearman'))
pls.kfcv$lxly.cor[['cor.lxlxhat']] <- diag(cor(pls.res$TExPosition.Data$lx,pls.kfcv$lx.hat,method='spearman'))
pls.kfcv$lxly.cor[['cor.lylyhat']] <- diag(cor(pls.res$TExPosition.Data$ly,pls.kfcv$ly.hat,method='spearman'))

#Average variance explained across folds and save.
lvx_var <- mean(pls.var.x)
lvy_var <- mean(pls.var.y)
outmat <- data.frame(c(lvx_var,lvy_var),row.names=c('X','Y'))
colnames(outmat) <- c('Variance Explained')
outfile <- paste0(outpath,'cv_var_explained.csv')
write.csv(outmat,outfile)

#Save all of these outputs.
cvout <- rbind(pls.kfcv$lxly.cor[['cor.lxhatlyhat']],
               pls.kfcv$lxly.cor[['cor.lxlxhat']],
               pls.kfcv$lxly.cor[['cor.lylyhat']])
rownames(cvout) <- c('lxhat_lyhat',
                     'lx_lxhat',
                     'ly_lyhat')
outfile <- paste0(outpath,'cvout.csv')
write.csv(cvout,outfile)

#Produce max and min values.
lxh_lim <- c()
lx_lim <- c()
lyh_lim <- c()
ly_lim <- c()
for (i in svlist) {
  
  #Extract latent variables.
  clxh <- pls.kfcv$lx.hat[,i]
  clx <- pls.res$TExPosition.Data$lx[,i]
  clyh <- pls.kfcv$ly.hat[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  
  #Extract max and min.
  lxh_lim <- c(lxh_lim,max(clxh),min(clxh))
  lx_lim <- c(lx_lim,max(clx),min(clx))
  lyh_lim <- c(lyh_lim,max(clyh),min(clyh))
  ly_lim <- c(ly_lim,max(cly),min(cly))
}
lxh_max <- max(lxh_lim)
lxh_min <- min(lxh_lim)
lx_max <- max(lx_lim)
lx_min <- min(lx_lim)
lyh_max <- max(lyh_lim)
lyh_min <- min(lyh_lim)
ly_max <- max(ly_lim)
ly_min <- min(ly_lim)
lxh_diff <- lxh_max - lxh_min
lx_diff <- lx_max - lx_min
lyh_diff <- lyh_max - lyh_min
ly_diff <- ly_max - ly_min

#Plot scatterplots and do bivariate normality tests.
testlabs <- c("mardia_skew","mardia_kurtosis","hz", "royston", "dh","energy")
ntest <- length(testlabs)
pr_norm_lxhlyh <- matrix(NA,nrow=nsv,ncol=ntest)
rownames(pr_norm_lxhlyh) <- paste0('LV',as.character(svlist))
colnames(pr_norm_lxhlyh) <- testlabs
pr_norm_lxlxh <- pr_norm_lxhlyh
pr_norm_lylyh <- pr_norm_lxhlyh
clc_txtsize <- 11
for (i in 1:nsv) {
  
  #Extract latent variables.
  clxh <- pls.kfcv$lx.hat[,i]
  clx <- pls.res$TExPosition.Data$lx[,i]
  clyh <- pls.kfcv$ly.hat[,i]
  cly <- pls.res$TExPosition.Data$ly[,i]
  
  #Combine into comparisons.
  lxhlyh <- data.frame('LVX_H'=clxh,'LVY_H'=clyh)
  lxhlyh_labs <- paste0(colnames(lxhlyh),as.character(i))
  lxlxh <- data.frame('LVX'=clx,'LVX_H'=clxh)
  lxlxh_labs <- paste0(colnames(lxlxh),as.character(i))
  lylyh <- data.frame('LVY'=cly,'LVY_H'=clyh)
  lylyh_labs <- paste0(colnames(lylyh),as.character(i))
  
  #Plot.
  xdiff <- lxh_diff*0.01
  ydiff <- lyh_diff*0.01
  plt_lxhlyh <- ggplot(aes(clxh,clyh),data=lxhlyh) +
    geom_point(size=1) +
    geom_smooth(method=lm,se=F) +
    scale_x_continuous(limits=c(lxh_min-xdiff,lxh_max+xdiff)) +
    scale_y_continuous(limits=c(lyh_min-ydiff,lyh_max+ydiff)) +
    labs(x=lxhlyh_labs[1],y=lxhlyh_labs[2]) +
    theme(axis.text.x=element_text(size=clc_txtsize),
          axis.text.y=element_text(size=clc_txtsize))
  xdiff <- lx_diff*0.01
  ydiff <- lxh_diff*0.01
  plt_lxlxh <- ggplot(aes(clx,clxh),data=lxlxh) +
    geom_point(size=1) +
    geom_smooth(method=lm,se=F) +
    scale_x_continuous(limits=c(lx_min-xdiff,lx_max+xdiff)) +
    scale_y_continuous(limits=c(lxh_min-ydiff,lxh_max+ydiff)) +
    labs(x=lxlxh_labs[1],y=lxlxh_labs[2]) +
    theme(axis.text.x=element_text(size=clc_txtsize),
          axis.text.y=element_text(size=clc_txtsize))
  xdiff <- ly_diff*0.01
  ydiff <- lyh_diff*0.01
  plt_lylyh <- ggplot(aes(cly,clyh),data=lylyh) +
    geom_point(size=1) +
    geom_smooth(method=lm,se=F) +
    scale_x_continuous(limits=c(ly_min-xdiff,ly_max+xdiff)) +
    scale_y_continuous(limits=c(lyh_min-ydiff,lyh_max+ydiff)) +
    labs(x=lylyh_labs[1],y=lylyh_labs[2]) +
    theme(axis.text.x=element_text(size=clc_txtsize),
          axis.text.y=element_text(size=clc_txtsize))
  
  #Save.
  outfile <- paste0(outpath,'LV',as.character(i),'_LXHLYH.jpg')
  ggsave(outfile,plt_lxhlyh,units='in',width=4,height=3)
  outfile <- paste0(outpath,'LV',as.character(i),'_LXLXH.jpg')
  ggsave(outfile,plt_lxlxh,units='in',width=4,height=3)
  outfile <- paste0(outpath,'LV',as.character(i),'_LYLYH.jpg')
  ggsave(outfile,plt_lylyh,units='in',width=4,height=3)
  
  #Set current bivariate relationship.
  c_XY <- lxhlyh
  c_XY_lab <- 'LXHLYH'
  
  #Generate bivariate QQ plot and save.
  outfile <- paste0(outpath,'LV',as.character(i),'_',c_XY_lab,'_qq.jpg')
  jpeg(file=outfile)
  mvn(c_XY,multivariatePlot='qq')
  dev.off()
  
  #Do each test.
  for (ctest in testlabs) {
    
    #Whether indexing needs to be specific.
    if (grepl('mardia',ctest,fixed=T)) {
      
      #Skew or kurtosis value.
      if (ctest=='mardia_skew') {
        ctest_raw <- 'mardia'
        ctest_lab <- 'Mardia Skewness'
      } else if (ctest=='mardia_kurtosis') {
        ctest_raw <- 'mardia'
        ctest_lab <- 'Mardia Kurtosis'
      }
      
      #Do test.
      fitted <- mvn(c_XY,mvnTest=ctest_raw)$multivariateNormality
      fac <- fitted[fitted$Test==ctest_lab,'p value']
      cpr <- as.numeric(levels(fac)[fac])
    } else {
      
      #Do test.
      fitted <- mvn(c_XY,mvnTest=ctest)$multivariateNormality
      cpr <- fitted[1,'p value']
    }
    
    #Save it.
    pr_norm_lxhlyh[paste0('LV',as.character(i)),ctest] <- cpr
  }
  
  #Set current bivariate relationship.
  c_XY <- lxlxh
  c_XY_lab <- 'LXLXH'
  
  #Generate bivariate QQ plot and save.
  outfile <- paste0(outpath,'LV',as.character(i),'_',c_XY_lab,'_qq.jpg')
  jpeg(file=outfile)
  mvn(c_XY,multivariatePlot='qq')
  dev.off()
  
  #Do each test.
  for (ctest in testlabs) {
    
    #Whether indexing needs to be specific.
    if (grepl('mardia',ctest,fixed=T)) {
      
      #Skew or kurtosis value.
      if (ctest=='mardia_skew') {
        ctest_raw <- 'mardia'
        ctest_lab <- 'Mardia Skewness'
      } else if (ctest=='mardia_kurtosis') {
        ctest_raw <- 'mardia'
        ctest_lab <- 'Mardia Kurtosis'
      }
      
      #Do test.
      fitted <- mvn(c_XY,mvnTest=ctest_raw)$multivariateNormality
      fac <- fitted[fitted$Test==ctest_lab,'p value']
      cpr <- as.numeric(levels(fac)[fac])
    } else {
      
      #Do test.
      fitted <- mvn(c_XY,mvnTest=ctest)$multivariateNormality
      cpr <- fitted[1,'p value']
    }
    
    #Save it.
    pr_norm_lxlxh[paste0('LV',as.character(i)),ctest] <- cpr
  }
  
  #Set current bivariate relationship.
  c_XY <- lylyh
  c_XY_lab <- 'LYLYH'
  
  #Generate bivariate QQ plot and save.
  outfile <- paste0(outpath,'LV',as.character(i),'_',c_XY_lab,'_qq.jpg')
  jpeg(file=outfile)
  mvn(c_XY,multivariatePlot='qq')
  dev.off()
  
  #Do each test.
  for (ctest in testlabs) {
    
    #Whether indexing needs to be specific.
    if (grepl('mardia',ctest,fixed=T)) {
      
      #Skew or kurtosis value.
      if (ctest=='mardia_skew') {
        ctest_raw <- 'mardia'
        ctest_lab <- 'Mardia Skewness'
      } else if (ctest=='mardia_kurtosis') {
        ctest_raw <- 'mardia'
        ctest_lab <- 'Mardia Kurtosis'
      }
      
      #Do test.
      fitted <- mvn(c_XY,mvnTest=ctest_raw)$multivariateNormality
      fac <- fitted[fitted$Test==ctest_lab,'p value']
      cpr <- as.numeric(levels(fac)[fac])
    } else {
      
      #Do test.
      fitted <- mvn(c_XY,mvnTest=ctest)$multivariateNormality
      cpr <- fitted[1,'p value']
    }
    
    #Save it.
    pr_norm_lylyh[paste0('LV',as.character(i)),ctest] <- cpr
  }
}

#Threshold p-values and save.
outlab <- 'LXHLYH'
pr_norm <- pr_norm_lxhlyh
pr_norm_thres <- pr_norm
pr_norm_thres[pr_norm < 0.05] <- NA
pr_norm_thres <- signif(pr_norm_thres,2)
outfile <- paste0(outpath,'binormtest_',outlab,'_p.csv')
write.csv(data.frame(pr_norm),outfile)
outfile <- paste0(outpath,'binormtest_',outlab,'_pthres.csv')
write.csv(data.frame(pr_norm_thres),outfile)

#Threshold p-values and save.
outlab <- 'LXLXH'
pr_norm <- pr_norm_lxlxh
pr_norm_thres <- pr_norm
pr_norm_thres[pr_norm < 0.05] <- NA
pr_norm_thres <- signif(pr_norm_thres,2)
outfile <- paste0(outpath,'binormtest_',outlab,'_p.csv')
write.csv(data.frame(pr_norm),outfile)
outfile <- paste0(outpath,'binormtest_',outlab,'_pthres.csv')
write.csv(data.frame(pr_norm_thres),outfile)

#Threshold p-values and save.
outlab <- 'LYLYH'
pr_norm <- pr_norm_lylyh
pr_norm_thres <- pr_norm
pr_norm_thres[pr_norm < 0.05] <- NA
pr_norm_thres <- signif(pr_norm_thres,2)
outfile <- paste0(outpath,'binormtest_',outlab,'_p.csv')
write.csv(data.frame(pr_norm),outfile)
outfile <- paste0(outpath,'binormtest_',outlab,'_pthres.csv')
write.csv(data.frame(pr_norm_thres),outfile)
