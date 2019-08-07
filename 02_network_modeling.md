Network modeling
================

These analyses were used to generate the network plots in figure 6, as
well
as

# Analysis of CST Occurrence vs Immunological Parameters - quasi-Poisson Duration

Iterates through every pairwise combination of CST and immunological
parameter (either IST or metacluster abundance at one of the three
timepoints) and tests for a significant associations between the immune
parameter and the number of days a subject spends in the CST.

``` r
#Start time
start.time <- Sys.time()

library(readr)

type <- "duration"

sites = c("REC", "NAS")

for (site in sites) {
  
  #Read in metadata and CST duration data
  mapping <- read_tsv(file.path(refined, sprintf("%s_Surv_Mapping.txt", site)))
  
  raw_table <- read_tsv(file.path(refined, sprintf("%s_Duration_Input.txt", site)))
  
  mapping_length <- ncol(mapping)
  
  table_length <- ncol(raw_table)
  
  #Iterate through CSTs
  for (cst_i in 3:table_length) {
    
    tmp.time <- Sys.time()
    
    cst <- colnames(raw_table[,cst_i])
    
    cst_in <- data.frame(raw_table[, 1:2], raw_table[, cst_i])
    
    #Lists to store results
    var_names <- list()
    model_pvals <- list()
    k <- 1
    
    voi_names <- list()
    var_pvals <- list()
    var_coeffs <- list()
    vars <- list()
    renamed_vars <- list()
    j <- 1
    
    #Iterate through immune parameters
    for (var_i in 4:mapping_length) {
      
      var <- colnames(mapping[,var_i])
      
      mapping_in <- data.frame(mapping[, 1:3], mapping[, var_i])
      
      fac = FALSE
      
      #Join metadata with CST duration data on subject
      working_table <- merge(mapping_in, cst_in, by = "Subject")
      
      working_table <- working_table[complete.cases(working_table[,]), ]
      
      #If the immmune parameter is IST (vs metacluster abundance), remove any ISTs that occur less than ten times
      if (typeof(working_table[[var]]) == "character") {
        
        fac = TRUE
        
        tmp_table <- as.data.frame(table(working_table[[var]]))
        
        rare_cats <- list()
        n <- 1
        
        for (cat_i in 1:nrow(tmp_table)) {
          
          if (tmp_table[[2]][cat_i] < 10) {
            
            rare_cats[[n]] <- levels(tmp_table[[1]])[cat_i]
            n <- n + 1
            
          }
          
        }
        
        working_table <- working_table[!(working_table[[var]] %in% rare_cats), ]
        
      }
      
      #If the current CST occurs in <10% of the remaining subjects or all but one IST has been filtered out, skip this comparison
      if ((sum(working_table[[cst]] == 0) > nrow(working_table)*0.9) || ((fac) && (nrow(table(working_table[[var]])) < 2))) {
        
        next
        
      }
      
      working_table$cst_freq <- working_table[[cst]]
      
      #Convert the immune parameter to a factor if it's an IST, otherwise center and scale the metacluster abundance
      if (fac) {
        
        working_table$voi <- factor(working_table[[var]])
        
      } else {
        
        working_table$voi <- c(scale(working_table[[var]], center = TRUE, scale = TRUE))
        
      }
      
      #Center GA on 37 weeks (term)
      working_table$GAB <- ((working_table$gaBirth)/37) - 1
      working_table$ldays = log(working_table$Total)
      
      #Fit the model with and without the immune parameter (voi = Variable of Interest = the immune parameter)
      fit.full <- glm(cst_freq ~ MOD + GAB + offset(ldays) + voi, data = working_table, family = quasipoisson)
      fit.null <- glm(cst_freq ~ MOD + GAB + offset(ldays), data = working_table, family = quasipoisson)
      
      #Test for significance
      model_test <- anova(fit.null, fit.full, test = "Chisq")
      
      var_names[[k]] <- var
      model_pvals[[k]] <- model_test$`Pr(>Chi)`[2]
      k <- k + 1

      #Write the model summary
      if(! dir.exists(pth <- file.path(refined,  sprintf("%s_results", type))))
        dir.create(pth)
      
      cat("Current Variable: ", file = file.path(refined,  sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat(var, file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat('\n', file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      capture.output(summary(fit.full), file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat('\n###\n###\n###\n\n', file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      
      
      #If the immune parameter was significant, get the coefficients and p-values
      if (model_test$`Pr(>Chi)`[2] < 0.05) {
        
        s <- summary(fit.full)
        l <- grep("voi", rownames(s$coefficients))
        for (i in l) {
          
          voi_names[[j]] <- var
          vars[[j]] <- rownames(s$coefficients)[i]
          #Because ISTs were renamed after they were defined, provide the alternate name
          renamed_vars[[j]] <- switch(rownames(s$coefficients)[i], 'voi'='voi', 'voiICS_1'='ICS_8', 'voiICS_2'='ICS_6', 'voiICS_3'='ICS_4', 'voiICS_4'='ICS_2', 'voiICS_5'='ICS_3', 'voiICS_6'='ICS_1', 'voiICS_7'='ICS_5', 'voiICS_8'='ICS_7', 'voiTPHE_1'='TPHE_7',   'voiTPHE_2'='TPHE_2',   'voiTPHE_3'='TPHE_1',   'voiTPHE_4'='TPHE_3',   'voiTPHE_5'='TPHE_5',   'voiTPHE_6'='TPHE_4',   'voiTPHE_7'='TPHE_6')
          var_coeffs[[j]] <- s$coefficients[i,1]
          var_pvals[[j]] <- s$coefficients[i,4]
          j <- j + 1
          
        }
        
      }
        
    }
    
    #Multiple test correction, model level
    adj_model_pvals <- p.adjust(model_pvals, method = "fdr")
    #Assemble model level results
    results <- cbind(var_names, model_pvals, adj_model_pvals)
    #Write model level results
    write.csv(results, file.path(refined, sprintf("%s_results/%s_model_pvals.csv", type, cst)))
    
    #Do the same thing for variables
    adj_var_pvals <- p.adjust(var_pvals, method = "fdr")
    term_results <- cbind(voi_names, vars, renamed_vars, var_coeffs, var_pvals, adj_var_pvals)
    
    write.csv(term_results, file.path(refined,  sprintf("%s_results/%s_term_pvals_and_betas.csv", type, cst)))
    
    cat("\nTime taken to complete last CST: ")
    cat(difftime(Sys.time(), tmp.time, units = "mins"))
    
  }
  
}
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   TPHE_IST_Birth = col_character(),
    ##   TPHE_IST_Disch = col_character(),
    ##   TPHE_IST_1YR = col_character(),
    ##   ICS_IST_Birth = col_character(),
    ##   ICS_IST_Disch = col_character(),
    ##   ICS_IST_1YR = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   Total = col_integer(),
    ##   REC_1 = col_integer(),
    ##   REC_2 = col_integer(),
    ##   REC_3 = col_integer(),
    ##   REC_4 = col_integer(),
    ##   REC_5 = col_integer(),
    ##   REC_6 = col_integer(),
    ##   REC_7 = col_integer(),
    ##   REC_8 = col_integer(),
    ##   REC_9 = col_integer(),
    ##   REC_10 = col_integer(),
    ##   REC_11 = col_integer(),
    ##   REC_12 = col_integer(),
    ##   REC_13 = col_integer()
    ## )

    ## 
    ## Time taken to complete last CST: 0.05015535
    ## Time taken to complete last CST: 0.07282457
    ## Time taken to complete last CST: 0.1226355
    ## Time taken to complete last CST: 0.08694412
    ## Time taken to complete last CST: 0.05218252
    ## Time taken to complete last CST: 0.09292652
    ## Time taken to complete last CST: 0.06812287
    ## Time taken to complete last CST: 0.05140285
    ## Time taken to complete last CST: 0.0447908
    ## Time taken to complete last CST: 0.004313219
    ## Time taken to complete last CST: 0.0472049
    ## Time taken to complete last CST: 0.0488318
    ## Time taken to complete last CST: 0.04824179

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   TPHE_IST_Birth = col_character(),
    ##   TPHE_IST_Disch = col_character(),
    ##   TPHE_IST_1YR = col_character(),
    ##   ICS_IST_Birth = col_character(),
    ##   ICS_IST_Disch = col_character(),
    ##   ICS_IST_1YR = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   Total = col_integer(),
    ##   NAS_1 = col_integer(),
    ##   NAS_2 = col_integer(),
    ##   NAS_3 = col_integer(),
    ##   NAS_4 = col_integer(),
    ##   NAS_5 = col_integer(),
    ##   NAS_6 = col_integer(),
    ##   NAS_7 = col_integer(),
    ##   NAS_8 = col_integer(),
    ##   NAS_9 = col_integer(),
    ##   NAS_10 = col_integer(),
    ##   NAS_11 = col_integer(),
    ##   NAS_12 = col_integer(),
    ##   NAS_13 = col_integer()
    ## )

    ## 
    ## Time taken to complete last CST: 0.04688052
    ## Time taken to complete last CST: 5.356234
    ## Time taken to complete last CST: 0.0702706
    ## Time taken to complete last CST: 0.05398927
    ## Time taken to complete last CST: 0.0526436
    ## Time taken to complete last CST: 0.0510429
    ## Time taken to complete last CST: 0.04509755
    ## Time taken to complete last CST: 0.04682993
    ## Time taken to complete last CST: 0.04853682
    ## Time taken to complete last CST: 0.0488218
    ## Time taken to complete last CST: 0.04807422
    ## Time taken to complete last CST: 0.04590672
    ## Time taken to complete last CST: 0.0472614

``` r
#Print total runtime
cat("\nTotal runtime:")
```

    ## 
    ## Total runtime:

``` r
cat(difftime(Sys.time(), start.time, units = "hours"))
```

    ## 0.1126649

``` r
cat("\n\n\n")
```

# Analysis of CST Occurrence vs Immunological Parameters - Logistic Occurrence Probability

Iterates through every pairwise combination of CST and immunological
parameter (either IST or metacluster abundance at one of the three
timepoints) and tests for a significant associations between the immune
parameter and whether or not the CST occurs at all within a subject.

Same as the quasi-Poisson duration above, just fitting a different
model.

``` r
#Start time
start.time <- Sys.time()

library(readr)

type <- "logit"

sites = c("REC", "NAS")

for (site in sites) {
  #Read in metadata and CST occurrence data
  mapping <- read_tsv(file.path(refined, sprintf("%s_Surv_Mapping.txt", site)))
  
  raw_table <- read_tsv(file.path(refined, sprintf("%s_Logit_Input.txt", site)))
  
  mapping_length <- ncol(mapping)
  
  table_length <- ncol(raw_table)
  
  for (cst_i in 3:table_length) {
    
    tmp.time <- Sys.time()
    
    cst <- colnames(raw_table[,cst_i])
    
    cst_in <- data.frame(raw_table[, 1:2], raw_table[, cst_i])
    
    #Lists to store results
    var_names <- list()
    model_pvals <- list()
    k <- 1
    
    voi_names <- list()
    var_pvals <- list()
    var_coeffs <- list()
    vars <- list()
    j <- 1
    
    for (var_i in 4:mapping_length) {
      
      var <- colnames(mapping[,var_i])
      
      mapping_in <- data.frame(mapping[, 1:3], mapping[, var_i])
      
      fac = FALSE
      
      #Merge metadata and CST occurrence data on subject
      working_table <- merge(mapping_in, cst_in, by = "Subject")
      
      working_table <- working_table[complete.cases(working_table[,]), ]
      
      #If the immune parameter is IST, remove any ISTs that occur less than ten times at the given time point.
      if (typeof(working_table[[var]]) == "character") {
        
        fac = TRUE
        
        tmp_table <- as.data.frame(table(working_table[[var]]))
        
        rare_cats <- list()
        n <- 1
        
        for (cat_i in 1:nrow(tmp_table)) {
          
          if (tmp_table[[2]][cat_i] < 10) {
            
            rare_cats[[n]] <- levels(tmp_table[[1]])[cat_i]
            n <- n + 1
            
          }
          
        }
        
        working_table <- working_table[!(working_table[[var]] %in% rare_cats), ]
        
      }
      
      #If the CST occurs in less than 10% of remaining subjects, skip this comparison
      if ((sum(working_table[[cst]] == "Yes") < nrow(working_table)*0.1) || ((fac) && (nrow(table(working_table[[var]])) < 2))) {
        
        next
        
      }
      
      working_table$cst_yn <- factor(working_table[[cst]])
      
      #Convert the immune parameter to a factor if it's an IST, otherwise center and scale the metacluster abundance
      if (fac) {
        
        working_table$voi <- factor(working_table[[var]])
        
      } else {
        
        working_table$voi <- c(scale(working_table[[var]], center = TRUE, scale = TRUE))
        
      }
      
      #Center and scale the number of observations per subjects and GA at birth
      working_table$Obs <- c(scale(working_table$Total, center = TRUE, scale = TRUE))
      working_table$GAB <- ((working_table$gaBirth)/37) - 1
      
      #Fit the model with and without the immune parameter
      fit.full <- glm(cst_yn ~ MOD + GAB + Obs + voi, data = working_table, family = binomial)
      fit.null <- glm(cst_yn ~ MOD + GAB + Obs, data = working_table, family = binomial)
      
      #Test for significance
      model_test <- anova(fit.null, fit.full, test = "Chisq")
      
      var_names[[k]] <- var
      model_pvals[[k]] <- model_test$`Pr(>Chi)`[2]
      k <- k + 1

      if(! dir.exists(pth <- file.path(refined,  sprintf("%s_results", type))))
        dir.create(pth)
      
      #Write the model summary
     cat("Current Variable: ", file = file.path(refined,  sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
     cat(var, file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
cat('\n', file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      capture.output(summary(fit.full), file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat('\n###\n###\n###\n\n', file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      
      #If it's significant, store additional information
      if (model_test$`Pr(>Chi)`[2] < 0.05) {
        
        s <- summary(fit.full)
        l <- grep("voi", rownames(s$coefficients))
        for (i in l) {
          
          voi_names[[j]] <- var
          vars[[j]] <- rownames(s$coefficients)[i]
          var_coeffs[[j]] <- s$coefficients[i,1]
          var_pvals[[j]] <- s$coefficients[i,4]
          j <- j + 1
          
        }
        
      }
        
    }
    
    #For the model itself and for the variables of interest, perform multiple test correction and save the results.
    adj_model_pvals <- p.adjust(model_pvals, method = "fdr")
    results <- cbind(var_names, model_pvals, adj_model_pvals)
    
    write.csv(results, file.path(refined, sprintf("%s_results/%s_model_pvals.csv", type, cst)))
    
    adj_var_pvals <- p.adjust(var_pvals, method = "fdr")
    term_results <- cbind(voi_names, vars, var_coeffs, var_pvals, adj_var_pvals)
    
    write.csv(term_results, file.path(refined, sprintf("%s_results/%s_term_pvals_and_betas.csv", type, cst)))
    
    cat("\nTime taken to complete last CST: ")
    cat(difftime(Sys.time(), tmp.time, units = "mins"))
    
  }
  
}
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   TPHE_IST_Birth = col_character(),
    ##   TPHE_IST_Disch = col_character(),
    ##   TPHE_IST_1YR = col_character(),
    ##   ICS_IST_Birth = col_character(),
    ##   ICS_IST_Disch = col_character(),
    ##   ICS_IST_1YR = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   Total = col_integer(),
    ##   REC_2 = col_character(),
    ##   REC_3 = col_character(),
    ##   REC_9 = col_character(),
    ##   REC_1 = col_character(),
    ##   REC_6 = col_character(),
    ##   REC_8 = col_character(),
    ##   REC_11 = col_character(),
    ##   REC_7 = col_character(),
    ##   REC_4 = col_character(),
    ##   REC_5 = col_character(),
    ##   REC_12 = col_character(),
    ##   REC_13 = col_character(),
    ##   REC_10 = col_character()
    ## )

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.0389856

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03748385
    ## Time taken to complete last CST: 0.03800955

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03767507

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.0363434
    ## Time taken to complete last CST: 0.03532879

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.0376293

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03729373

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03815925
    ## Time taken to complete last CST: 0.03909737

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03633742

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03683447
    ## Time taken to complete last CST: 0.003948534

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   TPHE_IST_Birth = col_character(),
    ##   TPHE_IST_Disch = col_character(),
    ##   TPHE_IST_1YR = col_character(),
    ##   ICS_IST_Birth = col_character(),
    ##   ICS_IST_Disch = col_character(),
    ##   ICS_IST_1YR = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   Total = col_integer(),
    ##   NAS_1 = col_character(),
    ##   NAS_7 = col_character(),
    ##   NAS_4 = col_character(),
    ##   NAS_5 = col_character(),
    ##   NAS_11 = col_character(),
    ##   NAS_6 = col_character(),
    ##   NAS_9 = col_character(),
    ##   NAS_13 = col_character(),
    ##   NAS_2 = col_character(),
    ##   NAS_8 = col_character(),
    ##   NAS_12 = col_character(),
    ##   NAS_3 = col_character(),
    ##   NAS_10 = col_character()
    ## )

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.0402618

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03910108
    ## Time taken to complete last CST: 0.03772487
    ## Time taken to complete last CST: 0.03833652
    ## Time taken to complete last CST: 0.0356443
    ## Time taken to complete last CST: 0.0353639

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03589335

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03476913
    ## Time taken to complete last CST: 0.0387526

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.0381929
    ## Time taken to complete last CST: 0.03593445

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    
    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03656717

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## 
    ## Time taken to complete last CST: 0.03734653

``` r
#Print total runtime
cat("\nTotal runtime:")
```

    ## 
    ## Total runtime:

``` r
cat(difftime(Sys.time(), start.time, units = "hours"))
```

    ## 0.01567103

``` r
cat("\n\n\n")
```

The warnings we see about 0/1 fitted values indicate that we need to use
LRT tests, and be careful with coefficient estimates at the
boundary.

# Analysis of CST Occurrence vs Immunological Parameters - Survival Time to Occurrence

Iterates through every pairwise combination of CST and immunological
parameter (either IST or metacluster abundance at one of the three
timepoints) and tests for a significant associations between the immune
parameter and how low it takes a subject to initially transition into a
CST (or out of a CST, in the cases of REC 1 and NAS 1)

Similar as the quasi-Poisson duration and logistic above, except fitting
a different model and not comparing a full and null model to assess
significance. The p-value of the term of interest in the full model only
is used to assess significance.

``` r
#Start time
start.time <- Sys.time()

library(survival)
library(icenReg)
```

    ## Warning: package 'icenReg' was built under R version 3.5.2

    ## Loading required package: Rcpp

    ## Warning: package 'Rcpp' was built under R version 3.5.2

    ## Loading required package: coda

    ## Warning: package 'coda' was built under R version 3.5.2

``` r
library(readr)

type <- "surv"

sites = c("REC", "NAS")

for (site in sites) {
  #Read in metadata and CST occurrence interval information
  mapping <- read_tsv(file.path(refined, sprintf("%s_Surv_Mapping.txt", site)))
  
  raw_table <- read_tsv(file.path(refined, sprintf("%s_Surv_Input.txt", site)))
  
  mapping_length <- ncol(mapping)
  
  table_length <- ncol(raw_table)
  
  for (cst_i in 3:table_length) {
    
    tmp.time <- Sys.time()
    
    cst <- colnames(raw_table[,cst_i])
    
    cst_in <- data.frame(raw_table[, 1:2], raw_table[, cst_i])
    
    cst_in <- cst_in[complete.cases(cst_in[,]),]
    #Lists to store results
    var_names <- list()
    voi_names <- list()
    var_coeffs <- list()
    var_exp_coeffs <- list()
    variable_pvals <- list()
    k <- 1
    
    for (var_i in 4:mapping_length) {
      
      var <- colnames(mapping[,var_i])
      
      mapping_in <- data.frame(mapping[, 1:3], mapping[, var_i])
      
      fac = FALSE
      #Join metadata and CST occurrence info by Subject
      working_table <- merge(mapping_in, cst_in, by = "Subject")
      
      working_table <- working_table[complete.cases(working_table[,]), ]
      
      #If the immune parameter is IST, remove any ISTs that occur less than ten times at the given time point.
      if (typeof(working_table[[var]]) == "character") {
        
        fac = TRUE
        
        tmp_table <- as.data.frame(table(working_table[[var]]))
        
        rare_cats <- list()
        n <- 1
        
        for (cat_i in 1:nrow(tmp_table)) {
          
          if (tmp_table[[2]][cat_i] < 10) {
            
            rare_cats[[n]] <- levels(tmp_table[[1]])[cat_i]
            n <- n + 1
            
          }
          
        }
        
        working_table <- working_table[!(working_table[[var]] %in% rare_cats), ]
        
      }
      
      #If fewer than 20 subjects remain having entered the CST, skip this comparison
      if ((nrow(working_table) < 20) || ((fac) && (nrow(table(working_table[[var]])) < 2))) {
        
        next
        
      }
      
      working_table$cst_obs <- working_table[[cst]]
      
      #Convert the immune parameter to a factor if it's an IST, otherwise center and scale the metacluster abundance
      if (fac) {
        
        working_table$voi <- factor(working_table[[var]])
        
      } else {
        
        working_table$voi <- c(scale(working_table[[var]], center = TRUE, scale = TRUE))
        
      }
      
      #Center and scale GA at birth
      working_table$GAB <- ((working_table$gaBirth)/37) - 1
      
      #Define the interval in which initial transition into a CST must have occurred
      cst_surv <- Surv(working_table$PrevSampleDOL, working_table$cst_obs, type = "interval2")
      
      #Fit the model
      fit <- ic_par(cst_surv ~ MOD + GAB + voi, data = working_table, model = "aft", dist = "loglogistic")
      
      #Collect resulting parameters
      fit_summary <- summary(fit)
      
      for (p in 5:nrow(fit_summary$summaryParameters)) {
        
        var_names[[k]] <- var
        voi_names[[k]] <- rownames(fit_summary$summaryParameters)[p]
        var_coeffs[[k]] <- fit_summary$summaryParameters[p, 1]
        var_exp_coeffs[[k]] <- fit_summary$summaryParameters[p, 2]
        variable_pvals[[k]] <- fit_summary$summaryParameters[p, 5]
        k <- k + 1
        
      }
    
      #Write the model summary
      if(! dir.exists(pth <- file.path(refined,  sprintf("%s_results", type))))
  dir.create(pth)
      
      cat("Current Variable: ", file = file.path(refined,  sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
     cat(var, file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat('\n', file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      capture.output(summary(fit), file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat('\n###\n###\n###\n\n', file = file.path(refined, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
    }
    
    #Perform multiple test correction and write the results
    adj_variable_pvals <- p.adjust(variable_pvals, method = "fdr")
    results <- cbind(var_names, voi_names, var_coeffs, var_exp_coeffs, variable_pvals, adj_variable_pvals)
    
    write.csv(results, file.path(refined, sprintf("%s_results/%s_model_pvals.csv", type, cst)))
    
    cat("\nTime taken to complete last CST: ")
    cat(difftime(Sys.time(), tmp.time, units = "mins"))
    
  }
  
}
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   TPHE_IST_Birth = col_character(),
    ##   TPHE_IST_Disch = col_character(),
    ##   TPHE_IST_1YR = col_character(),
    ##   ICS_IST_Birth = col_character(),
    ##   ICS_IST_Disch = col_character(),
    ##   ICS_IST_1YR = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   PrevSampleDOL = col_integer(),
    ##   REC_2 = col_integer(),
    ##   REC_3 = col_integer(),
    ##   REC_9 = col_integer(),
    ##   REC_1 = col_integer(),
    ##   REC_6 = col_integer(),
    ##   REC_8 = col_integer(),
    ##   REC_11 = col_integer(),
    ##   REC_7 = col_integer(),
    ##   REC_4 = col_integer(),
    ##   REC_5 = col_integer(),
    ##   REC_12 = col_integer(),
    ##   REC_13 = col_integer(),
    ##   REC_10 = col_integer()
    ## )

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced
    
    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## 
    ## Time taken to complete last CST: 0.03383714
    ## Time taken to complete last CST: 0.03054432
    ## Time taken to complete last CST: 0.03068927
    ## Time taken to complete last CST: 0.04457994
    ## Time taken to complete last CST: 0.02803228
    ## Time taken to complete last CST: 0.02785612
    ## Time taken to complete last CST: 0.02729083
    ## Time taken to complete last CST: 0.02734911
    ## Time taken to complete last CST: 0.02776447
    ## Time taken to complete last CST: 0.02740963
    ## Time taken to complete last CST: 0.02640962
    ## Time taken to complete last CST: 0.02767125
    ## Time taken to complete last CST: 0.003057337

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   TPHE_IST_Birth = col_character(),
    ##   TPHE_IST_Disch = col_character(),
    ##   TPHE_IST_1YR = col_character(),
    ##   ICS_IST_Birth = col_character(),
    ##   ICS_IST_Disch = col_character(),
    ##   ICS_IST_1YR = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   PrevSampleDOL = col_integer(),
    ##   NAS_1 = col_integer(),
    ##   NAS_7 = col_integer(),
    ##   NAS_4 = col_integer(),
    ##   NAS_5 = col_integer(),
    ##   NAS_11 = col_integer(),
    ##   NAS_6 = col_integer(),
    ##   NAS_9 = col_integer(),
    ##   NAS_13 = col_integer(),
    ##   NAS_2 = col_integer(),
    ##   NAS_8 = col_integer(),
    ##   NAS_12 = col_integer(),
    ##   NAS_3 = col_integer(),
    ##   NAS_10 = col_integer()
    ## )

    ## 
    ## Time taken to complete last CST: 0.02896803
    ## Time taken to complete last CST: 0.03063723
    ## Time taken to complete last CST: 0.03217995
    ## Time taken to complete last CST: 0.0308751
    ## Time taken to complete last CST: 0.02984157
    ## Time taken to complete last CST: 0.02787715
    ## Time taken to complete last CST: 0.02841572
    ## Time taken to complete last CST: 0.02923811
    ## Time taken to complete last CST: 0.02755238
    ## Time taken to complete last CST: 0.02678901
    ## Time taken to complete last CST: 0.02796894
    ## Time taken to complete last CST: 0.0214157
    ## Time taken to complete last CST: 0.02630523

``` r
#Print total runtime
cat("\nTotal runtime:")
```

    ## 
    ## Total runtime:

``` r
cat(difftime(Sys.time(), start.time, units = "hours"))
```

    ## 0.01250218

``` r
cat("\n\n\n")
```

# Targeted analysis of Alloiococcus abundance, Tphe5, and acute illness

``` r
library("readr")
library("tidyverse")
```

    ## ── Attaching packages ─────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.0.0     ✔ purrr   0.2.5
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
    ## ✔ tidyr   0.8.1     ✔ stringr 1.3.1
    ## ✔ ggplot2 3.0.0     ✔ forcats 0.3.0

    ## Warning: package 'tibble' was built under R version 3.5.2

    ## Warning: package 'dplyr' was built under R version 3.5.2

    ## ── Conflicts ────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(lme4)
```

    ## Warning: package 'lme4' was built under R version 3.5.2

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

``` r
library(geepack)

#Read in mapping file with metadata including Alloiococcus abundance, Tphe5 at birth or discharge, and acute illness
md.nas <- read_delim(file.path(refined, "NAS_Focused_Mapping.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   Subject = col_character(),
    ##   DOL = col_integer(),
    ##   MOD = col_character(),
    ##   gaBirth = col_double(),
    ##   Reads = col_integer(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character(),
    ##   TPHE_5 = col_character(),
    ##   Alloiococcus = col_double()
    ## )

``` r
rns <- md.nas[[1]]
md.nas <- md.nas[ , 2:ncol(md.nas)]
rownames(md.nas) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
#Turn it into a dataframe
md.nas <- data.frame(md.nas)
#Extract only post-discharge samples
md.post <- md.nas[which(md.nas$PostInitialDischarge == "Yes"), ]

#Center and scale numerical variables and convert categorical variables to factors
md.post$Subject <- factor(md.post$Subject)
md.post$GAB <- c(scale(md.post$gaBirth, center = TRUE, scale = TRUE))
md.post$nDOL <- c(scale(md.post$DOL, center = TRUE, scale = TRUE))
md.post$MOD <- factor(md.post$MOD)
md.post$IllnessVisit <- factor(md.post$IllnessVisit)
md.post$TPHE_5 <- factor(md.post$TPHE_5) #This is a binary variable with an affirmative value if the subject exhibited TPHE 5 at either birth or discharge

md.post$allo_freq <- round(md.post$Reads*md.post$Alloiococcus) #Convert Alloiococcus relative abundance to counts

#Test Alloiococcus and TPHE_5 as predictors of acute illness, both by themselves and jointly, controlling for confounders.
summary(glmer(IllnessVisit ~ Alloiococcus + nDOL + GAB + MOD + (1|Subject), data = md.post, family = binomial))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: IllnessVisit ~ Alloiococcus + nDOL + GAB + MOD + (1 | Subject)
    ##    Data: md.post
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    695.7    731.6   -340.9    681.7     1229 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -0.9847 -0.3009 -0.2192 -0.1586  6.0123 
    ## 
    ## Random effects:
    ##  Groups  Name        Variance Std.Dev.
    ##  Subject (Intercept) 0.9655   0.9826  
    ## Number of obs: 1236, groups:  Subject, 141
    ## 
    ## Fixed effects:
    ##                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -2.62185    0.25210 -10.400  < 2e-16 ***
    ## Alloiococcus      -3.52729    1.11728  -3.157 0.001594 ** 
    ## nDOL               0.41943    0.11638   3.604 0.000313 ***
    ## GAB                0.25975    0.15816   1.642 0.100532    
    ## MODVaginal_Breech -0.78757    1.25550  -0.627 0.530463    
    ## MODVaginal_Vertex  0.08336    0.30131   0.277 0.782053    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Allccc nDOL   GAB    MODV_B
    ## Alloiococcs -0.149                            
    ## nDOL        -0.078 -0.055                     
    ## GAB         -0.035  0.069  0.227              
    ## MODVgnl_Brc -0.136 -0.008  0.009  0.027       
    ## MODVgnl_Vrt -0.598 -0.020 -0.014 -0.154  0.109

``` r
summary(glmer(IllnessVisit ~ TPHE_5 + nDOL + GAB + MOD + (1|Subject), data = md.post, family = binomial))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: IllnessVisit ~ TPHE_5 + nDOL + GAB + MOD + (1 | Subject)
    ##    Data: md.post
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    709.0    744.8   -347.5    695.0     1229 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -0.9920 -0.3026 -0.2217 -0.1687  6.2358 
    ## 
    ## Random effects:
    ##  Groups  Name        Variance Std.Dev.
    ##  Subject (Intercept) 0.9307   0.9647  
    ## Number of obs: 1236, groups:  Subject, 141
    ## 
    ## Fixed effects:
    ##                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -2.95798    0.26197 -11.291  < 2e-16 ***
    ## TPHE_5Yes          0.54736    0.35835   1.527 0.126645    
    ## nDOL               0.41220    0.11557   3.567 0.000362 ***
    ## GAB                0.31411    0.15600   2.014 0.044056 *  
    ## MODVaginal_Breech -0.93397    1.21827  -0.767 0.443297    
    ## MODVaginal_Vertex  0.07927    0.29793   0.266 0.790185    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) TPHE_5 nDOL   GAB    MODV_B
    ## TPHE_5Yes   -0.342                            
    ## nDOL        -0.098  0.008                     
    ## GAB         -0.009 -0.080  0.252              
    ## MODVgnl_Brc -0.126 -0.029  0.005  0.039       
    ## MODVgnl_Vrt -0.587  0.076 -0.016 -0.153  0.107

``` r
summary(glmer(IllnessVisit ~ Alloiococcus + TPHE_5 + nDOL + GAB + MOD + (1|Subject), data = md.post, family = binomial))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: IllnessVisit ~ Alloiococcus + TPHE_5 + nDOL + GAB + MOD + (1 |  
    ##     Subject)
    ##    Data: md.post
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    696.8    737.7   -340.4    680.8     1228 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0085 -0.3039 -0.2184 -0.1567  6.0978 
    ## 
    ## Random effects:
    ##  Groups  Name        Variance Std.Dev.
    ##  Subject (Intercept) 0.9324   0.9656  
    ## Number of obs: 1236, groups:  Subject, 141
    ## 
    ## Fixed effects:
    ##                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)        -2.6988     0.2666 -10.123  < 2e-16 ***
    ## Alloiococcus       -3.4039     1.1158  -3.051 0.002284 ** 
    ## TPHE_5Yes           0.3582     0.3609   0.993 0.320888    
    ## nDOL                0.4203     0.1164   3.610 0.000307 ***
    ## GAB                 0.2489     0.1573   1.582 0.113557    
    ## MODVaginal_Breech  -0.8171     1.2431  -0.657 0.510984    
    ## MODVaginal_Vertex   0.1034     0.2997   0.345 0.730081    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Allccc TPHE_5 nDOL   GAB    MODV_B
    ## Alloiococcs -0.168                                   
    ## TPHE_5Yes   -0.353  0.105                            
    ## nDOL        -0.075 -0.052  0.005                     
    ## GAB         -0.015  0.062 -0.064  0.228              
    ## MODVgnl_Brc -0.119 -0.012 -0.027  0.008  0.030       
    ## MODVgnl_Vrt -0.583 -0.014  0.075 -0.014 -0.155  0.108

``` r
#Test TPHE_5 and acute illness as predictors of Alloiococcus abundance, both by themselves and jointly, controlling for confounders.
summary(geeglm(allo_freq ~ DOL*gaBirth + MOD + TPHE_5, id = Subject, corstr = "exchangeable", family = poisson(link = 'log'), data = md.post, offset = log(Reads)))
```

    ## 
    ## Call:
    ## geeglm(formula = allo_freq ~ DOL * gaBirth + MOD + TPHE_5, family = poisson(link = "log"), 
    ##     data = md.post, offset = log(Reads), id = Subject, corstr = "exchangeable")
    ## 
    ##  Coefficients:
    ##                     Estimate    Std.err   Wald Pr(>|W|)    
    ## (Intercept)       -1.7865106  1.4120761  1.601   0.2058    
    ## DOL                0.0060183  0.0040405  2.219   0.1364    
    ## gaBirth           -0.0123747  0.0395382  0.098   0.7543    
    ## MODVaginal_Breech  1.0138070  0.3999270  6.426   0.0112 *  
    ## MODVaginal_Vertex  0.1036130  0.2249404  0.212   0.6451    
    ## TPHE_5Yes         -1.9505040  0.3689272 27.952 1.24e-07 ***
    ## DOL:gaBirth       -0.0001946  0.0001174  2.749   0.0973 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Estimated Scale Parameters:
    ##             Estimate Std.err
    ## (Intercept)    15666    2282
    ## 
    ## Correlation: Structure = exchangeable  Link = identity 
    ## 
    ## Estimated Correlation Parameters:
    ##       Estimate Std.err
    ## alpha   0.2465 0.04085
    ## Number of clusters:   141   Maximum cluster size: 16

``` r
summary(geeglm(allo_freq ~ DOL*gaBirth + MOD + IllnessVisit, id = Subject, corstr = "exchangeable", family = poisson(link = 'log'), data = md.post, offset = log(Reads)))
```

    ## 
    ## Call:
    ## geeglm(formula = allo_freq ~ DOL * gaBirth + MOD + IllnessVisit, 
    ##     family = poisson(link = "log"), data = md.post, offset = log(Reads), 
    ##     id = Subject, corstr = "exchangeable")
    ## 
    ##  Coefficients:
    ##                    Estimate   Std.err  Wald Pr(>|W|)    
    ## (Intercept)       -1.380401  1.424449  0.94    0.333    
    ## DOL                0.004885  0.003962  1.52    0.218    
    ## gaBirth           -0.030411  0.039890  0.58    0.446    
    ## MODVaginal_Breech  0.838724  0.482997  3.02    0.082 .  
    ## MODVaginal_Vertex  0.206798  0.244865  0.71    0.398    
    ## IllnessVisitYes   -0.907277  0.215164 17.78  2.5e-05 ***
    ## DOL:gaBirth       -0.000155  0.000115  1.79    0.181    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Estimated Scale Parameters:
    ##             Estimate Std.err
    ## (Intercept)    16563    2520
    ## 
    ## Correlation: Structure = exchangeable  Link = identity 
    ## 
    ## Estimated Correlation Parameters:
    ##       Estimate Std.err
    ## alpha    0.287  0.0488
    ## Number of clusters:   141   Maximum cluster size: 16

``` r
summary(geeglm(allo_freq ~ DOL*gaBirth + MOD + TPHE_5 + IllnessVisit, id = Subject, corstr = "exchangeable", family = poisson(link = 'log'), data = md.post, offset = log(Reads)))
```

    ## 
    ## Call:
    ## geeglm(formula = allo_freq ~ DOL * gaBirth + MOD + TPHE_5 + IllnessVisit, 
    ##     family = poisson(link = "log"), data = md.post, offset = log(Reads), 
    ##     id = Subject, corstr = "exchangeable")
    ## 
    ##  Coefficients:
    ##                    Estimate   Std.err  Wald Pr(>|W|)    
    ## (Intercept)       -1.729464  1.378299  1.57    0.210    
    ## DOL                0.005705  0.003906  2.13    0.144    
    ## gaBirth           -0.014884  0.038727  0.15    0.701    
    ## MODVaginal_Breech  1.005486  0.393484  6.53    0.011 *  
    ## MODVaginal_Vertex  0.173316  0.226981  0.58    0.445    
    ## TPHE_5Yes         -1.908447  0.410353 21.63  3.3e-06 ***
    ## IllnessVisitYes   -0.905316  0.215667 17.62  2.7e-05 ***
    ## DOL:gaBirth       -0.000180  0.000114  2.49    0.115    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Estimated Scale Parameters:
    ##             Estimate Std.err
    ## (Intercept)    15115    2158
    ## 
    ## Correlation: Structure = exchangeable  Link = identity 
    ## 
    ## Estimated Correlation Parameters:
    ##       Estimate Std.err
    ## alpha    0.254  0.0427
    ## Number of clusters:   141   Maximum cluster size: 16
