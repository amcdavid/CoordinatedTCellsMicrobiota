library(readr)
library(dplyr)
library(broom)
library(survival)
library(icenReg)

refined = 'intermediates'


# Take a big ol' table
# names of response columns and predictor columns

# Iterate through pairs of columns in ist_tbl and microbe_tbl, adjusting for extra_tbl
# All keyed by Subject
# null_formula should include all confounders
# regression_fcn(formula, data) returns an object that permits an anova test and tidy method
# result_id?
pairwise_regression = function(extra_tbl, ist_tbl, microbe_tbl, null_formula, regression_fcn, result_id){
  predict_tbl = inner_join(extra_tbl, ist_tbl, by = 'Subject')
  response_tbl = left_join(predict_tbl['Subject'], microbe_tbl, by = 'Subject')
  #predict_tbl = left_join(response_tbl['Subject'], predict_tbl, by = 'Subject')
  # response_tbl / predict_tbl should be in lockstep.
  ist_tbl = select(ist_tbl, -Subject)
  response_tbl = select(response_tbl, -Subject)
  
  
  table_length <- ncol(response_tbl) 
  full_formula = update(null_formula, . ~ . + voi)
  
  
  if(! dir.exists(pth <- file.path(refined,  sprintf("%s_results", result_id))))
    dir.create(pth)
  
  all_results = list()
  i = 1
  
  warn_trapped_fcn = purrr::quietly(regression_fcn)
  for( cst_i in seq_len(table_length)){
    response_nm = names(response_tbl)[cst_i]
    for( voi_i in names(ist_tbl)){
      this_data = tibble(response = response_tbl[[cst_i]], voi = predict_tbl[[voi_i]]) %>% 
        bind_cols(predict_tbl[names(extra_tbl)]) %>% na.omit()
      if(is.numeric(this_data$voi)){
        this_data$voi = scale(this_data$voi)
      }
      try({
      fit.null = warn_trapped_fcn(null_formula, this_data)
      fit.full = warn_trapped_fcn(full_formula, this_data)
      model_test = anova(fit.null$result, fit.full$result, test = "Chisq")
      full_results = tidy(fit.full$result) %>% mutate(response = response_nm,
                                               voi = voi_i, produced_warning = length(fit.full$warnings) > 0)
      anova_p = model_test$`Pr(>Chi)`[2]
      anova_result = tibble(term = 'anova', p.value = anova_p, response = response_nm, voi = voi_i)
      all_results[[i]] = bind_rows(full_results, anova_result)
      })
      i = i + 1
    }
  }
 bind_rows(all_results )
}

to_fit = expand.grid(site = c('REC', 'NAS'), #microbiome sampling site
                     model_type = c('logit', 'duration', 'time_to_event'), #type of association. logit = odds of CST ever, duration = number of days in CST, time_to_event = rate of first entry to CST
                     specification = c('baseline', # baseline = MOD, GAB,
                                       'full'),# full = MOD, GAB, ABX, MILK
                     stringsAsFactors = FALSE
                     ) 

logit =  function(form_, data_) glm(form_, family = 'binomial', data = data_)
duration = function(form_, data_) glm(form_, family = 'quasipoisson', data = data_)
time_to_event = function(form_, data_){
  data_ = filter(data_, !is.na(response))
  response = Surv(data_$PrevSampleDOL, data_$response, type = "interval2")
  data_$response = NULL
  ic_par(form_, data = data_, model = "aft", dist = "loglogistic")
}

nabx_polish = read_csv('data/antibiotic_exposure.csv') %>% 
  tidyr::pivot_wider(Subject, names_from = 'discharge', values_from = 'Number of systemic antibiotic') %>% 
  rename(abx_hospital = 'FALSE', abx_discharge = 'TRUE')


milk_months = read_csv(file.path('intermediates', 'milk_subject.csv')) %>% select(Subject, milk_months)

hospital_humilk = read_csv('data/milk_hospital.csv') %>% rename(milk_perinatal = `Any Human Milk Perinatal`)
milk = full_join(milk_months, hospital_humilk)

fit_scenario_result = list()
for(i in seq_len(nrow(to_fit))){
  this_fit = to_fit[i,]
  ist_extra_tbl <- read_tsv(file.path(refined, sprintf("%s_Surv_Mapping.txt", this_fit$site)))
  
  # Get response data and
  if(this_fit$model_type == 'logit'){ # Logistic
     microbe_tbl <- read_tsv(file.path(refined, sprintf("%s_Logit_Input.txt", this_fit$site))) %>%
       mutate_at(.vars = vars(-Subject, -Total), .funs = function(x) x=='Yes')
  } else if(this_fit$model_type == 'duration'){  # Duration
    microbe_tbl =  read_tsv(file.path(refined, sprintf("%s_Duration_Input.txt", this_fit$site)))
  } else if(this_fit$model_type == 'time_to_event'){ #Survival
    microbe_tbl = read_tsv(file.path(refined, sprintf("%s_Surv_Input.txt",  this_fit$site)))
  }
  
  regr_fcn = get(this_fit$model_type)
  
  # baseline
  this_formula = formula(response ~  MOD + gaBirth)
  if(this_fit$specification == 'full'){
    this_formula = update(this_formula, . ~ . +  milk_perinatal + milk_months + abx_hospital + abx_discharge)
  }
  
  # Total used as offset
  if(this_fit$model_type %in% c('logit', 'duration')){
    ist_extra_tbl = left_join(ist_extra_tbl, microbe_tbl[c('Total', 'Subject')], by = 'Subject')
    microbe_tbl = select(microbe_tbl, -Total)
    this_formula = update(this_formula, . ~. + offset(log(Total)))
  } else{
    ist_extra_tbl = left_join(microbe_tbl[c('PrevSampleDOL', 'Subject')], ist_extra_tbl, by = 'Subject')
    microbe_tbl = select(microbe_tbl, -PrevSampleDOL)
  }
  
  ist_tbl = ist_extra_tbl %>% select(contains('TPHE'), contains('CD4'), contains('CD8'), contains('ICS'), 'Subject')
  extra_tbl = ist_extra_tbl[,c('Subject', setdiff(names(ist_extra_tbl), names(ist_tbl)))] %>%
    left_join(milk) %>%
    left_join(nabx_polish)
  
  fit_scenario_result[[i]] = pairwise_regression(extra_tbl, ist_tbl, microbe_tbl, null_formula = this_formula, regression_fcn = regr_fcn, do.call(paste, c(list(sep = '_'), this_fit)))

  
 
}



# Total (duration) should go into extra_tbl, and be an offset.
