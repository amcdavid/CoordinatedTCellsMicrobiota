Network modeling
================

These analyses were used to generate the network plots in figure 6, as
well as

# Analysis of CST Occurrence vs Immunological Parameters - quasi-Poisson Duration

Iterates through every pairwise combination of CST and immunological
parameter (either IST or metacluster abundance at one of the three
timepoints) and tests for a significant associations between the immune
parameter in adjusted and unadjusted models.

``` r
to_fit = expand.grid(site = c('REC', 'NAS'), #microbiome sampling site
                     model_type = c('logit', 'duration', 'time_to_event'), #type of association. logit = odds of CST ever, duration = number of days in CST, time_to_event = rate of first entry to CST
                     specification = c('mod_only', 
                                       'baseline', # MOD, GAB
                                       'full'),# full = MOD, GAB, ABX, MILK
                     stringsAsFactors = FALSE
                     ) 

logit =  function(form_, data_) arm::bayesglm(form_, family = 'binomial', data = data_)
duration = function(form_, data_) arm::bayesglm(form_, family = 'quasipoisson', data = data_)
time_to_event = function(form_, data_){
  data_ = dplyr::filter(data_, !is.na(response))
  data_$response = Surv(data_$PrevSampleDOL, data_$response, type = "interval2")
  res = ic_par(form_, data = data_, model = "aft", dist = "loglogistic")
  res
}

knitr::kable(to_fit)
```

| site | model\_type     | specification |
|:-----|:----------------|:--------------|
| REC  | logit           | mod\_only     |
| NAS  | logit           | mod\_only     |
| REC  | duration        | mod\_only     |
| NAS  | duration        | mod\_only     |
| REC  | time\_to\_event | mod\_only     |
| NAS  | time\_to\_event | mod\_only     |
| REC  | logit           | baseline      |
| NAS  | logit           | baseline      |
| REC  | duration        | baseline      |
| NAS  | duration        | baseline      |
| REC  | time\_to\_event | baseline      |
| NAS  | time\_to\_event | baseline      |
| REC  | logit           | full          |
| NAS  | logit           | full          |
| REC  | duration        | full          |
| NAS  | duration        | full          |
| REC  | time\_to\_event | full          |
| NAS  | time\_to\_event | full          |

Scenarios to fit.

``` r
nabx_polish = read_csv('data/antibiotic_exposure.csv') %>% 
  tidyr::pivot_wider(Subject, names_from = 'discharge', values_from = 'Number of systemic antibiotic') %>% 
  rename(abx_hospital = 'FALSE', abx_discharge = 'TRUE') %>%
  mutate(abx_hospital = (abx_hospital - mean(abx_hospital))/sd(abx_hospital), abx_discharge = abx_discharge- mean(abx_discharge))
```

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   discharge = col_logical(),
    ##   group = col_character(),
    ##   `Number of systemic antibiotic` = col_double()
    ## )

``` r
milk_months = read_csv(file.path('intermediates', 'milk_subject.csv')) %>% select(Subject, milk_months)
```

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   months_surveyed = col_double(),
    ##   milk_months = col_double(),
    ##   any_milk_quarter = col_double()
    ## )

``` r
hospital_humilk = read_csv('data/milk_hospital.csv') %>% rename(milk_perinatal = `Any Human Milk Perinatal`)
```

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   `Any Human Milk Perinatal` = col_logical()
    ## )

``` r
milk = full_join(milk_months, hospital_humilk) %>% mutate(milk_months = milk_months - mean(milk_months, na.rm = TRUE))
```

    ## Joining, by = "Subject"

Set up confounders (Antibiotics and Milk). Both are centered. Hospital
antibiotics (days) are Z-scored so represent SD increases.

``` r
network_results_path = file.path(refined, 'network_results')
if(!dir.exists(network_results_path)) dir.create(network_results_path)

fit_scenario_result = list()
for(i in seq_len(nrow(to_fit))){
  this_fit = to_fit[i,]
  ist_extra_tbl <- read_tsv(file.path(refined, sprintf("%s_Surv_Mapping.txt", this_fit$site))) %>% mutate(gaBirth = gaBirth - 37)
  
  # Get response data and
  if(this_fit$model_type == 'logit'){ # Logistic
     microbe_tbl <- read_tsv(file.path(refined, sprintf("%s_Logit_Input.txt", this_fit$site))) %>%
       mutate_at(.vars = vars(-Subject, -Total), .funs = function(x) x=='Yes')
  } else if(this_fit$model_type == 'duration'){  # Duration
    microbe_tbl =  read_tsv(file.path(refined, sprintf("%s_Duration_Input.txt", this_fit$site)))
    
  } else if(this_fit$model_type == 'time_to_event'){ #Survival
    microbe_tbl = read_tsv(file.path(refined, sprintf("%s_Surv_Input.txt",  this_fit$site)))
  }
  
  joined_tbl = left_join(microbe_tbl, ist_extra_tbl, by = 'Subject')
  regr_fcn = get(this_fit$model_type)
  
  # mod_only
  this_formula = formula(response ~  MOD)
  if(this_fit$specification == 'full'){
    this_formula = update(this_formula, . ~ . +  milk_perinatal + milk_months + abx_hospital + abx_discharge + gaBirth)
  } else if(this_fit$specification == 'baseline'){
     this_formula = update(this_formula, . ~ . + gaBirth)
  }
  
  # Total used as offset
  if(this_fit$model_type %in% c('logit', 'duration')){
    this_formula = update(this_formula, . ~. + offset(log(Total)))
    
  } else{
    ist_extra_tbl = left_join(microbe_tbl[c('PrevSampleDOL', 'Subject')], ist_extra_tbl, by = 'Subject')
    microbe_tbl = select(microbe_tbl, -PrevSampleDOL)
  }
  
  ist_fields = ist_extra_tbl %>% select(contains('TPHE'), contains('CD4'), contains('CD8'), contains('ICS')) %>% names()
  microbe_fields = microbe_tbl  %>% select(contains('REC'), contains('NAS')) %>% names()
  joined_tbl2 = joined_tbl %>%
    left_join(milk) %>%
    left_join(nabx_polish)
  
  message("Fitting ", paste0(this_fit, collapse = "-"))
  fit_scenario_result[[i]] = pairwise_regression(joined_tbl2, ist_fields = ist_fields, microbe_fields = microbe_fields, null_formula = this_formula, regression_fcn = regr_fcn)
  result_file = stringr::str_c(do.call(paste, c(list(sep = '_'), this_fit)), "_results.csv")
   
  write_csv(fit_scenario_result[[i]], file.path(network_results_path, result_file))
 
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
    ##   Total = col_double(),
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

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-logit-mod_only

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
    ##   Total = col_double(),
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

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-logit-mod_only

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
    ##   Total = col_double(),
    ##   REC_1 = col_double(),
    ##   REC_2 = col_double(),
    ##   REC_3 = col_double(),
    ##   REC_4 = col_double(),
    ##   REC_5 = col_double(),
    ##   REC_6 = col_double(),
    ##   REC_7 = col_double(),
    ##   REC_8 = col_double(),
    ##   REC_9 = col_double(),
    ##   REC_10 = col_double(),
    ##   REC_11 = col_double(),
    ##   REC_12 = col_double(),
    ##   REC_13 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-duration-mod_only

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
    ##   Total = col_double(),
    ##   NAS_1 = col_double(),
    ##   NAS_2 = col_double(),
    ##   NAS_3 = col_double(),
    ##   NAS_4 = col_double(),
    ##   NAS_5 = col_double(),
    ##   NAS_6 = col_double(),
    ##   NAS_7 = col_double(),
    ##   NAS_8 = col_double(),
    ##   NAS_9 = col_double(),
    ##   NAS_10 = col_double(),
    ##   NAS_11 = col_double(),
    ##   NAS_12 = col_double(),
    ##   NAS_13 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-duration-mod_only

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
    ##   PrevSampleDOL = col_double(),
    ##   REC_2 = col_double(),
    ##   REC_3 = col_double(),
    ##   REC_9 = col_double(),
    ##   REC_1 = col_double(),
    ##   REC_6 = col_double(),
    ##   REC_8 = col_double(),
    ##   REC_11 = col_double(),
    ##   REC_7 = col_double(),
    ##   REC_4 = col_double(),
    ##   REC_5 = col_double(),
    ##   REC_12 = col_double(),
    ##   REC_13 = col_double(),
    ##   REC_10 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-time_to_event-mod_only

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

    ## Error in solve.default(fit$hessian) : 
    ##   Lapack routine dgesv: system is exactly singular: U[1,1] = 0

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

    ## Error in solve.default(fit$hessian) : 
    ##   Lapack routine dgesv: system is exactly singular: U[4,4] = 0

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Error in solve.default(fit$hessian) : 
    ##   Lapack routine dgesv: system is exactly singular: U[1,1] = 0

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Warning in sqrt(diag(fit$var)): NaNs produced

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
    ##   PrevSampleDOL = col_double(),
    ##   NAS_1 = col_double(),
    ##   NAS_7 = col_double(),
    ##   NAS_4 = col_double(),
    ##   NAS_5 = col_double(),
    ##   NAS_11 = col_double(),
    ##   NAS_6 = col_double(),
    ##   NAS_9 = col_double(),
    ##   NAS_13 = col_double(),
    ##   NAS_2 = col_double(),
    ##   NAS_8 = col_double(),
    ##   NAS_12 = col_double(),
    ##   NAS_3 = col_double(),
    ##   NAS_10 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-time_to_event-mod_only

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
    ##   Total = col_double(),
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

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-logit-baseline

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
    ##   Total = col_double(),
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

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-logit-baseline

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
    ##   Total = col_double(),
    ##   REC_1 = col_double(),
    ##   REC_2 = col_double(),
    ##   REC_3 = col_double(),
    ##   REC_4 = col_double(),
    ##   REC_5 = col_double(),
    ##   REC_6 = col_double(),
    ##   REC_7 = col_double(),
    ##   REC_8 = col_double(),
    ##   REC_9 = col_double(),
    ##   REC_10 = col_double(),
    ##   REC_11 = col_double(),
    ##   REC_12 = col_double(),
    ##   REC_13 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-duration-baseline

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
    ##   Total = col_double(),
    ##   NAS_1 = col_double(),
    ##   NAS_2 = col_double(),
    ##   NAS_3 = col_double(),
    ##   NAS_4 = col_double(),
    ##   NAS_5 = col_double(),
    ##   NAS_6 = col_double(),
    ##   NAS_7 = col_double(),
    ##   NAS_8 = col_double(),
    ##   NAS_9 = col_double(),
    ##   NAS_10 = col_double(),
    ##   NAS_11 = col_double(),
    ##   NAS_12 = col_double(),
    ##   NAS_13 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-duration-baseline

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
    ##   PrevSampleDOL = col_double(),
    ##   REC_2 = col_double(),
    ##   REC_3 = col_double(),
    ##   REC_9 = col_double(),
    ##   REC_1 = col_double(),
    ##   REC_6 = col_double(),
    ##   REC_8 = col_double(),
    ##   REC_11 = col_double(),
    ##   REC_7 = col_double(),
    ##   REC_4 = col_double(),
    ##   REC_5 = col_double(),
    ##   REC_12 = col_double(),
    ##   REC_13 = col_double(),
    ##   REC_10 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-time_to_event-baseline

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

    ## Error in solve.default(fit$hessian) : 
    ##   Lapack routine dgesv: system is exactly singular: U[1,1] = 0

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Error in solve.default(fit$hessian) : 
    ##   Lapack routine dgesv: system is exactly singular: U[1,1] = 0

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Error in solve.default(fit$hessian) : 
    ##   Lapack routine dgesv: system is exactly singular: U[1,1] = 0

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Warning in sqrt(diag(fit$var)): NaNs produced

    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level

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
    ##   PrevSampleDOL = col_double(),
    ##   NAS_1 = col_double(),
    ##   NAS_7 = col_double(),
    ##   NAS_4 = col_double(),
    ##   NAS_5 = col_double(),
    ##   NAS_11 = col_double(),
    ##   NAS_6 = col_double(),
    ##   NAS_9 = col_double(),
    ##   NAS_13 = col_double(),
    ##   NAS_2 = col_double(),
    ##   NAS_8 = col_double(),
    ##   NAS_12 = col_double(),
    ##   NAS_3 = col_double(),
    ##   NAS_10 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-time_to_event-baseline

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
    ##   Total = col_double(),
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

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-logit-full

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
    ##   Total = col_double(),
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

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-logit-full

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
    ##   Total = col_double(),
    ##   REC_1 = col_double(),
    ##   REC_2 = col_double(),
    ##   REC_3 = col_double(),
    ##   REC_4 = col_double(),
    ##   REC_5 = col_double(),
    ##   REC_6 = col_double(),
    ##   REC_7 = col_double(),
    ##   REC_8 = col_double(),
    ##   REC_9 = col_double(),
    ##   REC_10 = col_double(),
    ##   REC_11 = col_double(),
    ##   REC_12 = col_double(),
    ##   REC_13 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-duration-full

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
    ##   Total = col_double(),
    ##   NAS_1 = col_double(),
    ##   NAS_2 = col_double(),
    ##   NAS_3 = col_double(),
    ##   NAS_4 = col_double(),
    ##   NAS_5 = col_double(),
    ##   NAS_6 = col_double(),
    ##   NAS_7 = col_double(),
    ##   NAS_8 = col_double(),
    ##   NAS_9 = col_double(),
    ##   NAS_10 = col_double(),
    ##   NAS_11 = col_double(),
    ##   NAS_12 = col_double(),
    ##   NAS_13 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-duration-full

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
    ##   PrevSampleDOL = col_double(),
    ##   REC_2 = col_double(),
    ##   REC_3 = col_double(),
    ##   REC_9 = col_double(),
    ##   REC_1 = col_double(),
    ##   REC_6 = col_double(),
    ##   REC_8 = col_double(),
    ##   REC_11 = col_double(),
    ##   REC_7 = col_double(),
    ##   REC_4 = col_double(),
    ##   REC_5 = col_double(),
    ##   REC_12 = col_double(),
    ##   REC_13 = col_double(),
    ##   REC_10 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting REC-time_to_event-full

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

    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level
    ## Error in ic_par(form_, data = data_, model = "aft", dist = "loglogistic") : 
    ##   covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level

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
    ##   PrevSampleDOL = col_double(),
    ##   NAS_1 = col_double(),
    ##   NAS_7 = col_double(),
    ##   NAS_4 = col_double(),
    ##   NAS_5 = col_double(),
    ##   NAS_11 = col_double(),
    ##   NAS_6 = col_double(),
    ##   NAS_9 = col_double(),
    ##   NAS_13 = col_double(),
    ##   NAS_2 = col_double(),
    ##   NAS_8 = col_double(),
    ##   NAS_12 = col_double(),
    ##   NAS_3 = col_double(),
    ##   NAS_10 = col_double()
    ## )

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

    ## Fitting NAS-time_to_event-full

    ## Warning in sqrt(diag(fit$var)): NaNs produced

Fit all the scenarios.

# Top interactions

``` r
clamp = function (x, modulus = 5) {
    x[x < -modulus] = -modulus
    x[x > modulus] = modulus
    x
}

sign_max = function(x) {
  sign(x[which.max(abs(x))])
}


library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 3.5.2

``` r
all_fits = bind_rows(fit_scenario_result, .id = 'scenario') %>% left_join(to_fit %>% mutate(scenario = as.character(seq_along(site))))
```

    ## Joining, by = "scenario"

``` r
all_fits  = all_fits %>% mutate(voi_full = str_replace(voi, 'TPHE_(C[1-9]+)', 'TPHE_CD4_\\1'),
                                voi_full = str_replace(voi_full, '(?<!TPHE)(_CD[48]_C[1-9]+)', '_ICS\\1'))

per_type_fdr = filter(all_fits, term == 'anova')  %>% group_by(model_type, specification)%>% mutate(fdr = p.adjust(p.value, 'fdr'))

per_type_sign = all_fits %>%  filter(stringr::str_detect(term, stringr::fixed('voi'))) %>% group_by(model_type, specification, voi, response) %>% summarize(sign_max = sign_max(estimate))

all_fits = all_fits %>%left_join(per_type_sign) %>% left_join(per_type_fdr[c('fdr', intersect(names(per_type_sign), names(per_type_fdr)))])
```

    ## Joining, by = c("response", "voi", "model_type", "specification")

    ## Joining, by = c("response", "voi", "model_type", "specification")

``` r
top = all_fits %>% filter(term == 'anova') %>% group_by(site, model_type) %>% arrange(p.value) %>% do(head(., n = 1))

top_coefs = semi_join(all_fits, top[c('scenario', 'voi', 'response')]) %>% anti_join(tibble(term = c('anova', '(Intercept)')))
```

    ## Joining, by = c("scenario", "response", "voi")

    ## Joining, by = "term"

``` r
ggplot(top_coefs, aes(x = term, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, color = clamp(-log10(p.value)))) + geom_pointrange() + facet_wrap(~model_type  + voi + response, scales = 'free') + coord_flip()
```

    ## Warning: Removed 1 rows containing missing values (geom_segment).

![](02_network_modeling_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
of_interest = tibble(response = c("NAS_8", 'REC_10'), voi = c('TPHE_IST_Disch', 'ICS_IST_1YR'))

write_csv(semi_join(all_fits, of_interest), path = 'intermediates/selected_network_effects.csv')
```

    ## Joining, by = c("response", "voi")

[Coefficients / etc for NAS\_8 and
REC\_10](intermediates/selected_network_effects.csv)

[Full results](intermediates/network_results)

# Network plot figures

``` r
enmatrix = function(x, rownames_from){
  y = x[,setdiff(names(x), rownames_from)]
  y = as.matrix(y)
  rownames(y) = x[[rownames_from]]
  y[is.na(y)] = 0
  y
}

per_type_fdr = left_join(per_type_fdr, per_type_sign) %>% mutate(signed_pval = clamp(-log10(p.value), 8) * sign_max, zeroed_pval = signed_pval * (fdr< .1))
```

    ## Joining, by = c("response", "voi", "model_type", "specification")

``` r
per_type_nest = per_type_fdr %>% select(voi_full, response, zeroed_pval, model_type, specification) %>% tidyr::nest(data = c(zeroed_pval, voi_full, response))

per_type_nest$adj_matrix = purrr::map(per_type_nest$data, ~ tidyr::pivot_wider(.x, voi_full, names_from = response, values_from = zeroed_pval) %>% enmatrix('voi_full'))

per_type_nest %>% rowwise() %>% mutate(n_edges = sum(abs(adj_matrix)>0)) %>% select(model_type, specification, n_edges)
```

    ## Source: local data frame [9 x 3]
    ## Groups: <by row>
    ## 
    ## # A tibble: 9 x 3
    ##   model_type    specification n_edges
    ##   <chr>         <chr>           <int>
    ## 1 logit         mod_only          100
    ## 2 duration      mod_only          166
    ## 3 time_to_event mod_only           39
    ## 4 logit         baseline            0
    ## 5 duration      baseline           41
    ## 6 time_to_event baseline           44
    ## 7 logit         full                0
    ## 8 duration      full               49
    ## 9 time_to_event full               15

``` r
mc_suffix = dplyr::tibble(Family = c('TPHE4', 'TPHE8', 'ICS4', 'ICS8'), 
                   suffix = c('t4', 't8', 'i4', 'i8'),
                   voi3_prefix = c('TPHE_CD4', 'TPHE_CD8', 'ICS_CD4', 'ICS_CD8'))
# Descriptive names for T cell subpop
metacluster_rn = readr::read_csv('intermediates/Metacluster Identities.csv') %>% select(-X5) %>% mutate(Family = Family %>% toupper(), Family = str_replace_all(Family, ' ', '')) %>% left_join(mc_suffix) %>% mutate(marker = str_c('Meta.Cluster_', Cluster, '_', suffix), voi3 = str_c(voi3_prefix, '_C', Cluster), Identity = make.unique(Identity))
```

    ## Warning: Missing column names filled in: 'X5' [5]

    ## Parsed with column specification:
    ## cols(
    ##   Family = col_character(),
    ##   Cluster = col_double(),
    ##   Category = col_character(),
    ##   Identity = col_character(),
    ##   X5 = col_logical()
    ## )

    ## Joining, by = "Family"

``` r
#Bipartite Graph of CST vs Immunome Results

#Load packages (Probably don't even need half of them - I tried a lot of things before I settled on this relatively simple solution.)
library(ggplot2)
library(network)
```

    ## network: Classes for Relational Data
    ## Version 1.16.1 created on 2020-10-06.
    ## copyright (c) 2005, Carter T. Butts, University of California-Irvine
    ##                     Mark S. Handcock, University of California -- Los Angeles
    ##                     David R. Hunter, Penn State University
    ##                     Martina Morris, University of Washington
    ##                     Skye Bender-deMoll, University of Washington
    ##  For citation information, type citation("network").
    ##  Type help("network-package") to get started.

``` r
library(GGally)
```

    ## 
    ## Attaching package: 'GGally'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     nasa

``` r
library(RColorBrewer)
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 3.5.2

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ tibble  2.1.3     ✓ purrr   0.3.3
    ## ✓ tidyr   1.0.2     ✓ forcats 0.5.0

    ## Warning: package 'tibble' was built under R version 3.5.2

    ## Warning: package 'tidyr' was built under R version 3.5.2

    ## Warning: package 'purrr' was built under R version 3.5.2

    ## Warning: package 'forcats' was built under R version 3.5.2

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
#Define constants

set.seed(11)

rec_csts = c("REC_4", "REC_1", "REC_2", "REC_9", "REC_10", "REC_5", "REC_8", "REC_6", "REC_3", "REC_13", "REC_7", "REC_11", "REC_12")
nas_csts = c("NAS_4", "NAS_1", "NAS_2", "NAS_9", "NAS_10", "NAS_5", "NAS_8", "NAS_6", "NAS_3", "NAS_13", "NAS_7", "NAS_11", "NAS_12")
#tests = c("logit", "surv", "binom", "alt_binom")
# tests = c("logit", "alt_surv")

to_plot = tibble(
  response = c('REC_10', 'NAS_3', 'NAS_8', 'NAS_9'),
  model_type = 'duration',
  specification = 'full')
  
to_plot_response = semi_join(per_type_nest, to_plot) %>%  select(-adj_matrix) %>% unnest(cols = c(data)) %>% ungroup() %>% 
  inner_join(to_plot) %>% filter(!is.na(zeroed_pval) & abs(zeroed_pval) > 0) 
```

    ## Joining, by = c("model_type", "specification")

    ## Joining, by = c("model_type", "specification", "response")

``` r
#Loop through the sites and make a figure for each as the target site (i.e. predicted by the CSTs of the other two sites)
for (test in seq_len(nrow(per_type_nest))){
  
  #Read in the site-specific adjacency matrix specifying associations between taxa and the CSTs of the other body sites
  adj.matrix = per_type_nest$adj_matrix[[test]]
    
  #Use that adjacency matrix to construct a network object that is bipartite
  net.obj <- network(adj.matrix, matrix.type = "bipartite", ignore.eval = FALSE, names.eval = "weights")
  
  #Color the CST vertices according to the body site they represent
  network::set.vertex.attribute(net.obj, "color", ifelse(net.obj %v% "vertex.names" %in% rec_csts, "red", ifelse(net.obj %v% "vertex.names" %in% nas_csts, "blue", ifelse(grepl("TPHE", net.obj %v% "vertex.names"), "orange", "green"))))
  
  #Shape the vertices according to the visit they represent
  network::set.vertex.attribute(net.obj, "shape", ifelse(grepl("Birth", net.obj %v% "vertex.names"), 19, ifelse(grepl("Disch", net.obj %v% "vertex.names"), 17, ifelse(grepl("OneYear", net.obj %v% "vertex.names"), 15, 18))))
  
  #Make non-significant associations invisible
  network::set.edge.attribute(net.obj, "alpha", ifelse(net.obj %e% "weight" == 0.0, 0, 1))
  
  if(sum(network::has.edges(net.obj)) == 0) next
  connected.net = network::get.inducedSubgraph(net.obj, which(network::has.edges(net.obj)))
  
  #Make a color gradient that will be used to color the edges
  rbPalFun0 <- colorRamp(c('darkblue', 'grey70',  'red'), space = 'Lab')
  
  #Map edge weights (signed log10 FDR adjusted p-vals) to colors in the gradient
  weights = network::get.edge.attribute(net.obj, "weights")
  rbPalFun = function(x){
    sx = x*.5/max(abs(x)) + .5
    #sx = (x-max(abs(weights))/(max(weights)-min(weights)))
    rgb(rbPalFun0(sx)/255)
  }
  ecols <- rbPalFun(weights)
  
  print(ggnet2(connected.net, color = "color", shape = "shape", edge.color = ecols, label = TRUE, label.size = 3, fontface = "bold", edge.size = 1,layout.par = list(niter = 1000))  + 
          ggtitle(paste0(per_type_nest[test,'model_type'], ':', per_type_nest[test,'specification'])))
  
  if(per_type_nest[test,'model_type'] == 'duration' && per_type_nest[test, 'specification'] == 'full'){
    to_plot_v = which((connected.net %v% "vertex.names") %in% to_plot$response)
    to_plot_neigh =  do.call(c, lapply(to_plot_v, network::get.neighborhood, x = connected.net))
    #subnet = network::get.inducedSubgraph(net.obj,  to_plot_v, alters = setdiff(seq_along(net.obj %v% "vertex.names"), to_plot_v))
    subnet = network::get.inducedSubgraph(connected.net,  union(to_plot_v, to_plot_neigh))
      ecols_sub <- rbPalFun(subnet %eattr% "weights")
  }

}
```

![](02_network_modeling_files/figure-gfm/networks-1.png)<!-- -->![](02_network_modeling_files/figure-gfm/networks-2.png)<!-- -->![](02_network_modeling_files/figure-gfm/networks-3.png)<!-- -->![](02_network_modeling_files/figure-gfm/networks-4.png)<!-- -->![](02_network_modeling_files/figure-gfm/networks-5.png)<!-- -->![](02_network_modeling_files/figure-gfm/networks-6.png)<!-- -->![](02_network_modeling_files/figure-gfm/networks-7.png)<!-- -->

``` r
#ggplot2.multiplot(nas_graph, rec_graph, thr_graph, cols = 2)

#ggsave(paste(workDir, sprintf("%s_graph.pdf", site), sep = "/"), width = 15, height = 15, units = "in")
```

``` r
 print(ggnet2(subnet, color = "color", shape = "shape", edge.color = ecols_sub, label = TRUE, label.size = 3, fontface = "bold", edge.size = 1,layout.par = list(niter = 1000))  + 
          ggtitle('duration:full (selected nodes)'))
```

![](02_network_modeling_files/figure-gfm/core_network-1.png)<!-- -->

``` r
cleanup_coefs = . %>%
  filter(str_detect(term, 'voi')) %>%
  mutate(is_anova = !str_detect(term, '^voi$'),
    voi2 = ifelse(is_anova, str_replace(term, 'voi', ''), voi_full),
         tp = str_extract(voi, 'OneYear|Disch|Birth'),
         voi3 = str_replace(voi2, '_?(OneYear(_)?|Disch(arge)?(_)?|Birth(_)?)', ''),
         cut_fdr = cut(-log10(ifelse(is.na(fdr), 1, fdr)), breaks = c(0, 1, 2, Inf), include.lowest = TRUE)
         ) %>%
  left_join(metacluster_rn) %>% mutate(metacluster_nm = coalesce(Identity, voi3))

to_plot_coefse_all = semi_join(all_fits, to_plot_response, by = c('model_type', 'specification', 'response')) %>% semi_join(to_plot_response, by = 'voi_full') %>% cleanup_coefs
```

    ## Joining, by = "voi3"

``` r
to_plot_coefse_sub = semi_join(all_fits, to_plot_response) %>% cleanup_coefs %>% mutate(metacluster_nm = factor(metacluster_nm, levels = unique(metacluster_nm[order(voi3)])))
```

    ## Joining, by = c("response", "model_type", "specification", "voi_full")
    ## Joining, by = "voi3"

``` r
plt = ggplot(to_plot_coefse_all, aes(y = clamp(estimate, 2), ymin = estimate-1.96*std.error, ymax = estimate + std.error*1.96, x = voi3, color = tp))+ facet_grid(is_anova ~ response,  space = 'free_y', scales = 'free', ) + xlab('') + ylab('Log Fold Change') + geom_pointrange(position = position_dodge(width = .4)) + geom_hline(yintercept = 0, lty = 2) + theme_minimal()  + coord_flip(ylim = c(-2, 2)) + scale_color_discrete('Time point')

plt + aes(alpha = cut_fdr) + scale_alpha_manual(values = c('[0,1]'=.3, '(1,2]'=.7, '(2,Inf]'=1))
```

![](02_network_modeling_files/figure-gfm/core_network-2.png)<!-- -->

``` r
plt %+% to_plot_coefse_sub
```

![](02_network_modeling_files/figure-gfm/core_network-3.png)<!-- -->

``` r
(plt  %+% to_plot_coefse_sub) + aes(x = metacluster_nm)
```

![](02_network_modeling_files/figure-gfm/core_networks_mc-1.png)<!-- -->

``` r
fsom_expr = readRDS('flowcytometry/intermediates/metacluster_mfi.rds')
descaled = bind_rows(fsom_expr$descaled, .id = 'population') %>% pivot_longer(cols = c(-cell_clustering, -population)) %>% filter(!is.na(value))
```

# Targeted analysis of Alloiococcus abundance, Tphe5, and acute illness

``` r
library("readr")

library(lme4)
```

    ## Warning: package 'lme4' was built under R version 3.5.2

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
library(geepack)
```

    ## Warning: package 'geepack' was built under R version 3.5.2

``` r
#Read in mapping file with metadata including Alloiococcus abundance, Tphe5 at birth or discharge, and acute illness
md.nas <- read_delim(file.path(refined, "NAS_Focused_Mapping.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   Subject = col_character(),
    ##   DOL = col_double(),
    ##   MOD = col_character(),
    ##   gaBirth = col_double(),
    ##   Reads = col_double(),
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

    ## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
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

    ## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
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

    ## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
    ##  Family: binomial  ( logit )
    ## Formula: IllnessVisit ~ Alloiococcus + TPHE_5 + nDOL + GAB + MOD + (1 |      Subject)
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
    ## Correlation structure = exchangeable 
    ## Estimated Scale Parameters:
    ## 
    ##             Estimate Std.err
    ## (Intercept)    15666    2282
    ##   Link = identity 
    ## 
    ## Estimated Correlation Parameters:
    ##       Estimate Std.err
    ## alpha   0.2465 0.04085
    ## Number of clusters:   141  Maximum cluster size: 16

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
    ## Correlation structure = exchangeable 
    ## Estimated Scale Parameters:
    ## 
    ##             Estimate Std.err
    ## (Intercept)    16563    2520
    ##   Link = identity 
    ## 
    ## Estimated Correlation Parameters:
    ##       Estimate Std.err
    ## alpha    0.287  0.0488
    ## Number of clusters:   141  Maximum cluster size: 16

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
    ## Correlation structure = exchangeable 
    ## Estimated Scale Parameters:
    ## 
    ##             Estimate Std.err
    ## (Intercept)    15115    2158
    ##   Link = identity 
    ## 
    ## Estimated Correlation Parameters:
    ##       Estimate Std.err
    ## alpha    0.254  0.0427
    ## Number of clusters:   141  Maximum cluster size: 16
