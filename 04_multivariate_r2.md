Multivariate *R*<sup>2</sup> Calculations
================
Andrew McDavid
07/06/2020

``` r
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE)
knitr::opts_chunk$set(dev = c('png', 'pdf'))
library(broom)
library(tidyverse)
library(ggbeeswarm)
library(car)

subject = read_csv('data/subject_covariates.csv') %>% transmute(ispreterm = preterm_weeks >= 0, Subject, preterm_weeks)
```

    ## Parsed with column specification:
    ## cols(
    ##   Gender = col_character(),
    ##   Race = col_character(),
    ##   `Birth Season` = col_character(),
    ##   preterm_weeks = col_double(),
    ##   auc14 = col_double(),
    ##   PRD = col_character(),
    ##   preg_antibiotics = col_character(),
    ##   mode_delivery = col_character(),
    ##   cchorio = col_character(),
    ##   preg_membrane_18hr = col_character(),
    ##   birth_wt_gms = col_double(),
    ##   `cmv test` = col_character(),
    ##   Subject = col_character()
    ## )

``` r
all_feats = read_csv('intermediates/all_tcell_features.csv') %>% left_join(subject)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Subject = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Joining, by = "Subject"

``` r
all_mb = read_csv('intermediates/microbiome_joined.csv') %>% left_join(subject)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Subject = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Joining, by = "Subject"

``` r
timeline = read_csv('data/subject_timeline.csv')
```

    ## Parsed with column specification:
    ## cols(
    ##   `Sequence Num` = col_double(),
    ##   DOL = col_double(),
    ##   cga = col_double(),
    ##   Subject = col_character()
    ## )

# Multivariate *R*<sup>2</sup> associations (figure 1)

``` r
loglik_mlm = function(object){
   res = object$residuals
  sig = cov(res)
  p = object$rank
  n = nrow(res)
  sum(mvtnorm::dmvnorm(res, mean = rep(0, ncol(res)), sigma = sig, log = TRUE))
}

aic_mlm = function(object){
   res = object$residuals
    p = object$rank
    N = nrow(res)
  -2*loglik_mlm(object) + 2*p
}


bic_mlm = function(object){
   res = object$residuals
    p = object$rank
    N = nrow(res)
  -2*loglik_mlm(object) + p*log(nrow(res))
}


as_char_form = function(x) paste0(as.character(x)[c(2, 1, 3:length(x))], collapse = ' ')

estimate_mlm_r2 = function(lhs_data, rhs_data, formula_, null_form = . ~ 1, response_label, predictor_label, adjusted = NA){
  # Drop NAs and constant columns
  lhs_data = as.matrix(lhs_data)
  #rhs_data = as.matrix(rhs_data)
  colVar = apply(lhs_data, 2, var, na.rm = TRUE)
  lhs_data = lhs_data[,colVar > 1e-4,drop = FALSE]
  colVar = apply(rhs_data, 2, var, na.rm = TRUE)
  rhs_data = rhs_data[,colVar > 1e-4, drop = FALSE]
  good_right = rowSums(is.na(rhs_data)) == 0
  good_left = rowSums(is.na(lhs_data)) == 0
  good_good = good_right & good_left
  lhs_data = lhs_data[good_good,,drop = FALSE]
  rhs_data = rhs_data[good_good,,drop = FALSE]
  
  # Left hand side is compositional, take ILR
  lhs_data_comp = compositions::ilr(compositions::clo(lhs_data))
  formula_ = as.formula(formula_)
  environment(formula_) = environment()
  full = lm(formula_, data = rhs_data)
  null = update(full, null_form)
  empty = update(full, . ~ 1)
  adjr2 = 1 - (sum(full$residuals^2)/full$df.residual) /  (sum(null$residuals^2)/null$df.residual)
  total_r2 = 1- (sum(full$residuals^2)/full$df.residual) /  (sum(empty$residuals^2)/empty$df.residual)
  anova_ = anova(full, null, test = 'Wilks', tol = 1e-6)
  stats = tibble(total_r2 = total_r2, adjr2 = adjr2, pval = anova_[2,'Pr(>F)'], response_label = response_label, predictor_label = predictor_label, n = min(nrow(rhs_data), nrow(lhs_data)), p_rhs = ncol(rhs_data), p_lhs = ncol(lhs_data), adjusted = adjusted, null_form = as_char_form(null_form), full_form = as_char_form(formula_), aic_full = aic_mlm(full), bic_full = bic_mlm(full))
  stats
}

push = function(x, y){
  if(missing(y)){
    x = list(x)
    } else{
      if(!is.list(x)) stop('Not list')
      x[[length(x)+1]] = y
    }
  x
}

all_mb_dol = all_mb %>% mutate(log2DOL = log2(DOL + 1 ))
rec_mb = all_mb_dol %>% select(ends_with('.rec')) %>% as.matrix()


all_feats_dol = all_feats %>% left_join(timeline[c('Subject',  'Sequence Num', 'DOL')], by = c('Subject',  'Sequence Num')) %>% mutate(log2DOL = log2(`DOL` + 1))
```

## Setup T cell with MB

Get Visits 1, 7, 19 from timeline, and their DOL

Merge onto MB timelines

``` r
closest_1_7_19 = timeline %>% filter(`Sequence Num` %in% c(1, 7, 19)) %>% select(`Sequence Num`, `Subject`, DOL_t = DOL) %>%  left_join(all_mb %>% select(`Subject`, DOL_m = DOL), by = 'Subject') 

closest_visit = closest_1_7_19 %>% mutate(diff = DOL_m -DOL_t)%>% group_by(`Subject`, `Sequence Num`) %>% arrange(abs(diff)) %>% mutate(rank_diff = seq_along(diff)) %>% filter(rank_diff == 1, abs(diff) < 28)

closest_visit %>% group_by(`Subject`) %>% summarize(n())
```

    ## # A tibble: 215 x 2
    ##    Subject `n()`
    ##    <chr>   <int>
    ##  1 C01D8       3
    ##  2 C0427       2
    ##  3 C04D3       2
    ##  4 C04F0       3
    ##  5 C0522       3
    ##  6 C054B       3
    ##  7 C0796       3
    ##  8 C08F3       2
    ##  9 C09B7       3
    ## 10 C09C2       3
    ## # … with 205 more rows

``` r
all_mb_closest = inner_join(all_mb[!duplicated(all_mb[c('Subject', 'DOL')]),], 
                            closest_visit %>% select(`Sequence Num.T` = `Sequence Num`, `Subject`, DOL = DOL_m))

mb_tcell_closest = inner_join(all_feats %>% rename(`Sequence Num.T` = `Sequence Num`), select(all_mb_closest, -`Sequence Num`, -ispreterm, -cga), by = c('Subject', 'Sequence Num.T'), suffix = c('.T', '.B'))

mb_tcell_closest %>% group_by(`Subject`) %>% summarize(n())
```

    ## # A tibble: 147 x 2
    ##    Subject `n()`
    ##    <chr>   <int>
    ##  1 C01D8       3
    ##  2 C0427       1
    ##  3 C04D3       2
    ##  4 C0522       3
    ##  5 C08F3       1
    ##  6 C09B7       3
    ##  7 C09C2       3
    ##  8 C0C1C       3
    ##  9 C0E4F       3
    ## 10 C1114       2
    ## # … with 137 more rows

``` r
pred_table = tibble(covariates = c('Term', 'DOL', 'PMA', 'T cell', 'NAS', 'REC', 'T cell', 'NAS', 'REC', 'PMA', 'DOL'),
      adjusted = c("no", "Term", "Term", "Term", "Term", "Term", 'PMA', 'PMA', 'PMA', 'no', 'PMA'))

r2_res = list()
for(i in seq_len(nrow(pred_table))){
  for(resp in c('T cell', 'NAS', 'REC')){
    pred = pred_table[[i, 'covariates']]
    adjusted = pred_table[[i, 'adjusted']]
    if(pred == resp) next
    if( (pred  == 'T cell' && resp %in% c('NAS', 'REC')) || (resp == 'T cell' && pred %in% c('NAS', 'REC'))) {
      data = mb_tcell_closest
    } else if (resp == 'T cell') {
      data = all_feats_dol
    } else{
      data = all_mb_dol
    }
    
    null_form = . ~  ispreterm
    if(pred == 'Term'){
      form_ = "lhs_data_comp ~ ispreterm"
      null_form = . ~ 1
      rhs_data = data['ispreterm']
    } else if(pred == 'DOL'){
      form_ = "lhs_data_comp ~ log2DOL"
      rhs_data = data[c('ispreterm', 'log2DOL')]
    } else if( pred == 'PMA'){
      form_ = "lhs_data_comp ~ cga"
      rhs_data = data[c('ispreterm', 'cga')]
    } else if (pred == 'T cell'){
      form_ = "lhs_data_comp ~ ."
      rhs_data = data %>% select(starts_with('Meta.Cluster'), ispreterm)
    } else if( pred == 'NAS'){
      form_ = "lhs_data_comp ~ ."
      rhs_data = data %>% select(ends_with('.nas'), ispreterm)
    } else if (pred == 'REC'){
      form_ = "lhs_data_comp ~ ."
      rhs_data = data %>% select(ends_with('.rec'), ispreterm)
    }
    
    # always adjust for preterm
    if(adjusted == 'Term' || adjusted == 'PMA'){
      form_ = str_c(form_, ' + ispreterm')
    }
    
    if(adjusted == 'GAB'){
      form_ = str_c(form_, ' + preterm_weeks')
    }
    
    if(adjusted == 'PMA'){
      form_ = str_c(form_, ' + cga')
      rhs_data = cbind(rhs_data, data['cga'])
      null_form = . ~ ispreterm + cga
    }
    
    if(resp == 'T cell'){
      lhs_data = data %>% select(starts_with('Meta.Cluster'))
    } else if(resp == 'NAS'){
      lhs_data = data %>% select(ends_with('.nas'))
    } else if( resp == 'REC'){
      lhs_data = data %>% select(ends_with('.rec'))
    }

    r2_res = push(r2_res, estimate_mlm_r2(lhs_data, rhs_data, form_, null_form, response_label = resp, predictor_label = pred, adjusted = adjusted))
  }
}
```

``` r
r2_res = bind_rows(r2_res) 

r2_res_graph = filter(r2_res, adjusted %in% c('no', 'Term', 'PMA'), !(adjusted == 'no' && predictor_label == "PMA"),
                      !(adjusted == 'PMA' && predictor_label == 'DOL')) %>% 
  mutate(predictor_label2 = factor(predictor_label, levels = rev(c('Term', 'DOL', 'PMA', 'T cell', 'NAS', 'REC'))),
         response_label = factor(response_label, levels = c('T cell', 'NAS', 'REC')), 
         pma_adjusted = factor(ifelse(adjusted == 'PMA', 'PMA', 'no'), levels = c('PMA', 'no')),
         pstar = case_when(pval < 1e-20 ~ '**',  pval < 1e-4 ~ '*', TRUE ~ ''))

#p_rhs = r2_res %>% arrange(predictor_label2, p_rhs) %>% split(f = .$predictor_label2) %>% map_dfr(~ .x[[1,]])


ggplot(r2_res_graph, aes(y = adjr2, x = predictor_label2, fill = pma_adjusted)) + geom_col(position = 'dodge') + coord_flip() + facet_grid(~response_label)+ geom_text(aes(y = adjr2 + .01, label = pstar), position = position_dodge(width = 1), size = 6) + ylab('Adjusted R2') + xlab('Predictor(s)') + scale_y_continuous(limits = c(0, .3), breaks = c(0, .1, .2)) + theme(legend.position = 'bottom') + scale_fill_discrete('Adjusted?', direction = 1) + theme_minimal()
```

![](04_multivariate_r2_files/figure-gfm/plot_mv_r2-1.png)<!-- -->

``` r
r2_res %>% dplyr::select(predictor_label, response_label, everything()) %>% mutate(adjr2 = round(adjr2, 3)) %>% write_csv(path = 'intermediates/di_results/r2_supp_table.csv')
```

## AIC (Smaller is better)

``` r
dplyr::filter(r2_res, predictor_label %in% c('Term', 'DOL', 'PMA')) %>% select(predictor_label, response_label, aic_full, adjusted) %>% tidyr::spread(predictor_label, aic_full) %>% knitr::kable()
```

| response\_label | adjusted |       DOL |       PMA |      Term |
|:----------------|:---------|----------:|----------:|----------:|
| NAS             | no       |        NA | 212265.38 | 213630.42 |
| NAS             | PMA      | 211872.17 |        NA |        NA |
| NAS             | Term     | 212418.66 | 212067.19 |        NA |
| REC             | no       |        NA | 599526.07 | 600603.30 |
| REC             | PMA      | 598630.16 |        NA |        NA |
| REC             | Term     | 599222.27 | 599005.28 |        NA |
| T cell          | no       |        NA |  10976.49 |  11465.76 |
| T cell          | PMA      |  10589.97 |        NA |        NA |
| T cell          | Term     |  10888.47 |  10773.92 |        NA |

### BIC (Smaller is better)

``` r
bic_r2 = dplyr::filter(r2_res, predictor_label %in% c('DOL', 'PMA'), adjusted %in% c('Term', 'PMA')) %>% 
  dplyr::select(predictor_label, response_label, bic_full, adjusted, adjr2, total_r2) %>% 
  dplyr::group_by(response_label) %>% 
  dplyr::mutate(log_norm = matrixStats::logSumExp(bic_full)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(`log10(p) vs worst model` = (bic_full-log_norm)*log10(exp(1)))

knitr::kable(bic_r2)
```

| predictor\_label | response\_label | bic\_full | adjusted |     adjr2 | total\_r2 | log\_norm | log10(p) vs worst model |
|:-----------------|:----------------|----------:|:---------|----------:|----------:|----------:|------------------------:|
| DOL              | T cell          |  10900.33 | Term     | 0.1657228 | 0.1843943 |  10900.33 |                 0.00000 |
| DOL              | NAS             | 212436.32 | Term     | 0.0550058 | 0.0745323 | 212436.32 |                 0.00000 |
| DOL              | REC             | 599240.26 | Term     | 0.0282133 | 0.0360097 | 599240.26 |                 0.00000 |
| PMA              | T cell          |  10785.78 | Term     | 0.1788257 | 0.1972038 |  10900.33 |               -49.75249 |
| PMA              | NAS             | 212084.85 | Term     | 0.0660061 | 0.0853053 | 212436.32 |              -152.64342 |
| PMA              | REC             | 599023.27 | Term     | 0.0306001 | 0.0383773 | 599240.26 |               -94.23715 |
| DOL              | T cell          |  10605.78 | PMA      | 0.0167490 | 0.2106498 |  10900.33 |              -127.92417 |
| DOL              | NAS             | 211895.71 | PMA      | 0.0065262 | 0.0912748 | 212436.32 |              -234.78286 |
| DOL              | REC             | 598654.15 | PMA      | 0.0038198 | 0.0420505 | 599240.26 |              -254.54379 |

``` r
readr::write_csv(bic_r2, 'intermediates/bic_supp_table.csv')
```

## 
