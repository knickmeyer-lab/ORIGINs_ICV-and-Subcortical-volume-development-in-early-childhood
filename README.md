# ORIGINs_ICV-and-Subcortical-volume-development-in-early-childhood

This repository contains the scripts that were used to do the analysis for the paper "Mapping Subcortical Brain Development and Cognition in Infancy and Early Childhood: A Global, Multi-Cohort Study "

The data is obtained from 8 cohorts. The data for 3 of the cohorts (UNC Early Brain Development in twins, UNC/UMN Baby Connectome Project,Longitudinal MRI Study of Infants at Risk for Autism (IBIS)) are available upon request for access from NDA. The other cohorts data are available upon request to the parent cohorts on IRB approval.

#########################################
# Functions for non-linear mixed models #
#########################################

fit_nlme <- function(formula_basic, formula_random, covariates, par_start, het_var, data, max_iter){
  # Grouping data by ID
  data_grouped <- groupedData(formula_basic, data = data)
  # Define group of heterogeneous variances
  eval(parse(text = paste0('wts <- varIdent(form = ~ 1 | ', het_var, ')')))
  s0 <- rep(0, length(covariates))
  # List of starting parameters
  start_list <- list(fixed = c(par_start[1], s0, par_start[2], s0, par_start[3]))
  # Define Asymptotic Regression Model
  model_formula <- as.formula(paste0(all.vars(formula_basic)[1], ' ~ SSasymp(', all.vars(formula_basic)[2], ', asym, r0, lrc)'))
  # Formula for fixed effects
  fixed_list <- list(as.formula(paste('asym ~ ', paste(covariates, collapse = '+'))), 
                     as.formula(paste('r0 ~ ', paste(covariates, collapse = '+'))), lrc ~ 1)
  # Try to fit model
  tmp <- nlme(model = model_formula, 
              data = data_grouped,
              fixed = fixed_list,
              random = formula_random,
              start = start_list,
              control = lmeControl(msMaxIter = max_iter),
              weights = wts,
              na.action = na.exclude)
  # Save model used
  tmp$call[["model"]] <- model_formula
  tmp$fixed_list <- fixed_list
  tmp$formula_random <- formula_random
  
  return(tmp)
}

# Function to obtain model predictions
predict_nlme <- function(model, newdata, level1 = 'ID'){
  require(nlme)
  newdata <- as.data.frame(newdata)
  mod_fix <- fixed.effects(model)
  mod_ran <- random.effects(model)
  x_name <- as.character(as.formula(as.character(model$call)[2])[[3]][2])
  x <- newdata[, x_name]
  asym_fix <- mod_fix[grep('asym', names(mod_fix))]
  asym_fix_names <- sapply(strsplit(names(asym_fix[-1]), '\\.'), function(z)z[2])
  for (i in asym_fix_names) {
    if(!i %in% colnames(newdata)) newdata[, i] <- 0
    newdata[, i] <- as.numeric(as.character(newdata[, i]))
  }
  asym_fix_pred <- asym_fix[1] + as.matrix(newdata[, asym_fix_names]) %*% as.vector(asym_fix[-1])
  r0_fix <- mod_fix[grep('r0', names(mod_fix))]
  r0_fix_names <- sapply(strsplit(names(r0_fix[-1]), '\\.'), function(z)z[2])
  for (i in r0_fix_names) newdata[, i] <- as.numeric(as.character(newdata[, i]))
  r0_fix_pred <- r0_fix[1] + as.matrix(newdata[, r0_fix_names]) %*% as.vector(r0_fix[-1])
  
  ID_col <- grep(level1, colnames(newdata), value = T)
  if(length(ID_col) > 0){
    asym_ran <- mod_ran[grep('asym', names(mod_ran))]
    asym_ran_pred <- asym_ran[as.character(newdata[,ID_col]),]
    r0_ran <- mod_ran[grep('r0', names(mod_ran))]
    r0_ran_pred <- r0_ran[as.character(newdata[, ID_col]),]
    asym_pred <- asym_fix_pred + asym_ran_pred
    r0_pred <- r0_fix_pred + r0_ran_pred
  } else {
    asym_pred <- asym_fix_pred
    r0_pred <- r0_fix_pred
  }
  
  SSasymp(input = x, Asym = asym_pred, R0 = r0_pred, lrc = mod_fix[grep('lrc', names(mod_fix))])
}


###############################
# Fit non-linear mixed models #
###############################

# load data and nlme package
load('volume_data.RData')
library(nlme)

# Starting initial parameters
par_start_list <- list('TLM' = c(13744.64506, 7977.90075, -6.27325),
                       'CDT' = c(6865.672587, 3532.31516, -6.23229),
                       'PTM' = c(8721.183910, 4576.25337, -6.34319),
                       'PLD' = c(2091.695183, 1094.77481, -6.08433),
                       'HPS' = c(5014.750573, 2204.84316, -5.78780),
                       'AGD' = c(2062.760786, 423.066260, -5.385317),
                       'ICV' = c(1.337962e+06, 3.900078e+05, -5.720917e+00))

# Basic formula
formula_basic_list <- list('TLM' = TLM ~ VAGE | ID, 
                           'CDT' = CDT ~ VAGE | ID,
                           'PTM' = PTM ~ VAGE | ID,
                           'PLD' = PLD ~ VAGE | ID,
                           'HPS' = HPS ~ VAGE | ID,
                           'AGD' = AGD ~ VAGE | ID,
                           'ICV' = ICV ~ VAGE | ID)

# string of covariates
covariates <- c('PRET', 'LIF', 'SEXM', 'LME', 'LBW')

# for loop to store models by volume
vol_models <- list()
for (i in 1:length(par_start_list)){
  vol_models[[i]] <- fit_nlme(formula_basic = formula_basic_list[[i]],
                                 formula_random = asym + r0 ~ 1,
                                 covariates = covariates,
                                 par_start = par_start_list[[i]],
                                 het_var = 'group_CA',
                                 data = vdata,
                                 max_iter = 1000)
  
  names(vol_models)[i] <- names(par_start_list)[i]
  print(names(par_start_list)[i])
}


#########################
# LRT test by covariate #
#########################

# Initial matrix
lrt_pvals <- matrix(NA, nrow = 7, ncol = 5)
colnames(lrt_pvals) <- covariates
rownames(lrt_pvals) <- names(par_start_list)
# For loop for LRT
for (j in 1:ncol(lrt_pvals)){
  for(i in 1:nrow(lrt_pvals)){
      tmp <- fit_nlme(formula_basic = formula_basic_list[[i]],
                         formula_random = asym + r0 ~ 1,
                         covariates = covariates[-j],
                         par_start = par_start_list[[i]],
                         het_var = 'group_CA',
                         data = vdata,
                         max_iter = 1000, 
                         attempts = 100)
      
    lrt_pvals[i, j] <- anova(vol_models[[i]], tmp)[2, 9]
    print(paste(value = lrt_pvals[i, j], i, j))
  }
}


###############################
# Fit linear mixed models #
###############################

# load data and lme4 package
load('mul_data.RData')
library(lme4)

# Define names of Muller scores
mul_raw_scores <- colnames(mul_raw)[9:13]

# List of covariates with cohort effect
covariates <- c('LBW', 'PRET', 'LME', 'SEXM', 'LIF', 'CHT')

# Fit mixed model for each score in a for loop
mul_raw_models <- list()
for (i in 1:length(mul_raw_scores)){
  mul_formula <- as.formula(paste(mul_raw_scores[i], '~ VAGE +', paste(covariates, collapse = '+'), '+ (1|ID)'))
  mul_raw_models[[i]] <- lmer(mul_formula, data = mul_raw, REML = F)
}
names(mul_raw_models) <- mul_raw_scores


#########################
# LRT test by covariate #
#########################

# Initial matrix
mul_raw_lrt <- matrix(NA, nrow = length(mul_raw_scores), ncol = length(covariates))
colnames(mul_raw_lrt) <- covariates
rownames(mul_raw_lrt) <- mul_raw_scores

# For loop for LRT
for (j in 1:length(mul_raw_scores)){
  for(i in 1:length(covariates)){
    tmp_formula <- as.formula(paste(mul_raw_scores[j], '~ VAGE +', paste(covariates[-i], collapse = '+'), ' + (1|ID)'))
    tmp <- lmer(tmp_formula, data = mul_raw, REML = F)
    mul_raw_lrt[j, i] <- anova(mul_raw_models[[j]], tmp)[2, 8]
  }
}


############################
# Obtain model predictions #
############################

# Create template for predictions based on original data
col_vol <- c('ID', 'CHT', 'LBW', 'PRET', 'LME', 'SEXM', 'LIF')
data_pred <- vdata[!duplicated(vdata[, col_vol]), col_vol]
CHT_means <- table(data_pred$CHT) / nrow(data_pred)

# Age of prediction (days after birth)
j = 730

# Volume predictions
for (i in 1:length(vol_models)) 
  data_pred[, paste0(names(vol_models)[i], '_vol')] <- predict_nlme(vol_models[[i]], cbind(data_pred, VAGE = j), level1 = 'ID')
  

# Mullen predictions
for (i in 1:length(mul_raw_models)) {
  if (!isSingular(mul_raw_models[[i]])) {
    coefs <- coef(mul_raw_models[[i]])$ID
    col_name <- paste0(names(mul_raw_models)[i], '_mul')
    data_pred[, col_name] <- NA
    
    for (rowi in 1:nrow(data_pred)) {
      coef_row <- rownames(coefs) %in% as.character(data_pred[rowi,'ID'])
      
      if (any(coef_row)) {
        betas <- coefs[coef_row,]
        x <- c(1, j, as.numeric(data_pred[rowi, c('LBW', 'PRET', 'LME', 'SEXM', 'LIF')]), 
               CHT_means[gsub('CHT','', grep('CHT', names(betas), value = T))])
        data_pred[rowi, col_name] <- sum(x * betas)
      }
    }
  }
}


# Mediation Analyis
# load data and mediation package

load(data_pred)
library(mediation)

b <-lm (vol~SEXM+LBW+LME+LIF+CHT, data = data_pred)
summary(b) 

c <-lm(mul_score~vol+SEXM+LBW+LME+LIF+CHT, data = data_pred)
summary(c)

model <- mediate(b, c, treat=(covariate), mediator=(vol), boot= TRUE,sims=10000)
summary(model)
