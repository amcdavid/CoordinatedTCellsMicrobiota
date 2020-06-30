CoordinatedTCellsMicrobiota
================

Data and analysis supporting “Aberrant newborn T cell and microbiota
developmental trajectories predict respiratory compromise during
infancy.” [Preprint](https://www.biorxiv.org/content/10.1101/736090v2)

## Contents

  - data: minimally processed data – microbiome OTU counts, T cell
    subpopulation abundances, subject covariates, etc
  - intermediates: derived tables and output from scripts

To reproduce some results, data will need to be downloaded [from
dbGAP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001347.v2.p1),
as some data posted here has been jittered, quantized or redacted to
reduce the risk of reidentification of subjects. The `Alias` key will
join subject-level tables to data in dbGAP.

## Scripts

  - 01\_make\_cst: Run DMN model on OTU and flowsom metaclusters and
    find best-fitting model
  - 01\_table1\_demographics: Summarize some demographic variables to
    generate portions of “table 1” in the paper.
  - 02\_network\_modeling: Test for pairwise associations between CST,
    IST and flowsom metacluster
  - 02\_ist\_cst\_results: Heatmaps of IST/CST composition, associations
    with gestational age and postmenstrual age
  - 03\_development\_index: prediction of PMA from microbiome and T cell
    subpopulations. Developmental index. Prediction of PRD. Associations
    of T cells subpopulation with PMA. Broken stick regression models.
    IST associations with perinatal exposures.
  - flowcytometry: scripts for processing T cell flow cytometry
    (TPHE/ICS) to produce flowsom “metacluster” populations. FCS files
    would need to be [downloaded from
    IMMPORT](https://browser.immport.org/browser?path=SDY1302) to run
    these scripts.

These don’t necessarily need to be run in order, as their inputs have
been saved under `intermediates`.

To install all dependencies run `devtools::install(CoordinatedTDeps, lib
= <some new directory>)` where you may wish to set `<some new
directory>` to be a path besides your default `.libPaths()` to keep
these packages isolated.

Run using R 3.5.1 and Bioconductor 3.7, and YMMV
    otherwise.

# Colophone

``` r
devtools::load_all('CoordinatedTDeps')
```

    ## Loading CoordinatedTDeps

``` r
devtools::session_info()
```

    ## Session info ------------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.5.1 (2018-07-02)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/New_York            
    ##  date     2020-06-29

    ## Packages ----------------------------------------------------------------------

    ##  package              * version    date       source                         
    ##  abind                  1.4-5      2016-07-21 CRAN (R 3.5.0)                 
    ##  assertthat             0.2.1      2019-03-21 cran (@0.2.1)                  
    ##  backports              1.1.5      2019-10-02 cran (@1.1.5)                  
    ##  base                 * 3.5.1      2018-07-05 local                          
    ##  BBmisc                 1.11       2017-03-10 CRAN (R 3.5.0)                 
    ##  beeswarm               0.2.3      2016-04-25 CRAN (R 3.5.0)                 
    ##  Biobase                2.40.0     2018-05-01 Bioconductor                   
    ##  BiocGenerics           0.26.0     2018-05-01 Bioconductor                   
    ##  boot                   1.3-20     2017-08-06 CRAN (R 3.5.1)                 
    ##  broom                  0.5.5      2020-02-29 cran (@0.5.5)                  
    ##  car                    3.0-7      2020-03-11 cran (@3.0-7)                  
    ##  carData                3.0-3      2019-11-16 cran (@3.0-3)                  
    ##  cellranger             1.1.0      2016-07-27 CRAN (R 3.5.0)                 
    ##  checkmate              2.0.0      2020-02-06 cran (@2.0.0)                  
    ##  cluster                2.0.7-1    2018-04-13 CRAN (R 3.5.0)                 
    ##  coda                   0.19-3     2019-07-05 CRAN (R 3.5.2)                 
    ##  codetools              0.2-16     2018-12-24 CRAN (R 3.5.2)                 
    ##  colorspace             1.3-2      2016-12-14 CRAN (R 3.5.0)                 
    ##  commonmark             1.5        2018-04-28 CRAN (R 3.5.0)                 
    ##  compiler               3.5.1      2018-07-05 local                          
    ##  ConsensusClusterPlus   1.44.0     2018-05-01 Bioconductor                   
    ##  CoordinatedTDeps     * 0.0.0.9000 <NA>       local                          
    ##  corpcor                1.6.9      2017-04-01 CRAN (R 3.5.0)                 
    ##  crayon                 1.3.4      2017-09-16 CRAN (R 3.5.0)                 
    ##  curl                   4.3        2019-12-02 cran (@4.3)                    
    ##  data.table             1.12.8     2019-12-09 cran (@1.12.8)                 
    ##  datasets             * 3.5.1      2018-07-05 local                          
    ##  DBI                    1.1.0      2019-12-15 cran (@1.1.0)                  
    ##  dbplyr                 1.4.2      2019-06-17 cran (@1.4.2)                  
    ##  DEoptimR               1.0-8      2016-11-19 cran (@1.0-8)                  
    ##  devtools               1.13.6     2018-06-27 CRAN (R 3.5.0)                 
    ##  digest                 0.6.25     2020-02-23 cran (@0.6.25)                 
    ##  DirichletMultinomial   1.22.0     2018-05-01 Bioconductor                   
    ##  dplyr                  0.8.5      2020-03-07 cran (@0.8.5)                  
    ##  evaluate               0.14       2019-05-28 cran (@0.14)                   
    ##  farver                 2.0.3      2020-01-16 cran (@2.0.3)                  
    ##  fastmap                1.0.1      2019-10-08 cran (@1.0.1)                  
    ##  fastmatch              1.1-0      2017-01-28 CRAN (R 3.5.0)                 
    ##  fda                    2.4.8.1    2020-02-19 CRAN (R 3.5.2)                 
    ##  flowCore               1.46.2     2018-09-12 Bioconductor                   
    ##  FlowSOM                1.12.0     2018-05-01 Bioconductor                   
    ##  flowStats              3.38.0     2018-05-01 Bioconductor                   
    ##  flowViz                1.44.0     2018-05-01 Bioconductor                   
    ##  flowWorkspace          3.28.2     2018-09-13 Bioconductor                   
    ##  forcats                0.5.0      2020-03-01 cran (@0.5.0)                  
    ##  foreach                1.4.8      2020-02-09 cran (@1.4.8)                  
    ##  foreign                0.8-71     2018-07-20 CRAN (R 3.5.0)                 
    ##  fs                     1.3.2      2020-03-05 cran (@1.3.2)                  
    ##  geepack                1.3-1      2019-12-13 cran (@1.3-1)                  
    ##  generics               0.0.2      2018-11-29 cran (@0.0.2)                  
    ##  gganimate              1.0.5      2020-02-09 cran (@1.0.5)                  
    ##  ggbeeswarm             0.6.0      2017-08-07 CRAN (R 3.5.0)                 
    ##  ggplot2                3.3.0      2020-03-05 cran (@3.3.0)                  
    ##  ggpubr                 0.2.5      2020-02-13 cran (@0.2.5)                  
    ##  ggrepel                0.8.2      2020-03-08 cran (@0.8.2)                  
    ##  ggsignif               0.6.0      2019-08-08 CRAN (R 3.5.2)                 
    ##  glue                   1.3.2      2020-03-12 cran (@1.3.2)                  
    ##  graph                  1.58.0     2018-05-01 Bioconductor                   
    ##  graphics             * 3.5.1      2018-07-05 local                          
    ##  grDevices            * 3.5.1      2018-07-05 local                          
    ##  grid                   3.5.1      2018-07-05 local                          
    ##  gridExtra              2.3        2017-09-09 CRAN (R 3.5.0)                 
    ##  gtable                 0.3.0      2019-03-25 cran (@0.3.0)                  
    ##  haven                  2.2.0      2019-11-08 cran (@2.2.0)                  
    ##  hexbin                 1.28.1     2020-02-03 cran (@1.28.1)                 
    ##  hms                    0.5.3      2020-01-08 cran (@0.5.3)                  
    ##  htmltools              0.4.0      2019-10-04 cran (@0.4.0)                  
    ##  htmlwidgets            1.5.1      2019-10-08 cran (@1.5.1)                  
    ##  httpuv                 1.5.2      2019-09-11 cran (@1.5.2)                  
    ##  httr                   1.4.1      2019-08-05 cran (@1.4.1)                  
    ##  icenReg                2.0.13     2019-12-16 cran (@2.0.13)                 
    ##  IDPmisc                1.1.20     2020-01-21 CRAN (R 3.5.2)                 
    ##  igraph                 1.2.5      2020-03-19 cran (@1.2.5)                  
    ##  IRanges                2.14.11    2018-08-24 Bioconductor                   
    ##  iterators              1.0.12     2019-07-26 cran (@1.0.12)                 
    ##  jsonlite               1.6.1      2020-02-02 cran (@1.6.1)                  
    ##  KernSmooth             2.23-15    2015-06-29 CRAN (R 3.5.1)                 
    ##  km.ci                  0.5-2      2009-08-30 CRAN (R 3.5.0)                 
    ##  KMsurv                 0.1-5      2012-12-03 CRAN (R 3.5.0)                 
    ##  knitr                  1.28       2020-02-06 cran (@1.28)                   
    ##  ks                     1.11.7     2020-02-11 CRAN (R 3.5.2)                 
    ##  later                  1.0.0      2019-10-04 cran (@1.0.0)                  
    ##  lattice                0.20-38    2018-11-04 CRAN (R 3.5.0)                 
    ##  latticeExtra           0.6-28     2016-02-09 cran (@0.6-28)                 
    ##  lazyeval               0.2.2      2019-03-15 cran (@0.2.2)                  
    ##  lifecycle              0.2.0      2020-03-06 cran (@0.2.0)                  
    ##  lme4                   1.1-21     2019-03-05 CRAN (R 3.5.2)                 
    ##  lubridate              1.7.4      2018-04-11 CRAN (R 3.5.0)                 
    ##  magrittr               1.5        2014-11-22 CRAN (R 3.5.0)                 
    ##  MASS                   7.3-50     2018-04-30 CRAN (R 3.5.1)                 
    ##  Matrix                 1.2-14     2018-04-13 CRAN (R 3.5.1)                 
    ##  matrixStats            0.54.0     2018-07-23 CRAN (R 3.5.0)                 
    ##  mclust                 5.4.1      2018-06-27 CRAN (R 3.5.0)                 
    ##  memoise                1.1.0      2018-03-13 Github (hadley/memoise@611cfad)
    ##  methods              * 3.5.1      2018-07-05 local                          
    ##  mgcv                   1.8-24     2018-06-23 CRAN (R 3.5.1)                 
    ##  mime                   0.9        2020-02-04 cran (@0.9)                    
    ##  minqa                  1.2.4      2014-10-09 CRAN (R 3.5.0)                 
    ##  mlr                    2.17.0     2020-01-10 cran (@2.17.0)                 
    ##  modelr                 0.1.6      2020-02-22 cran (@0.1.6)                  
    ##  munsell                0.5.0      2018-06-12 CRAN (R 3.5.0)                 
    ##  mvtnorm                1.0-8      2018-05-31 CRAN (R 3.5.0)                 
    ##  ncdfFlow               2.26.0     2018-05-01 Bioconductor                   
    ##  nlme                   3.1-137    2018-04-07 CRAN (R 3.5.1)                 
    ##  nloptr                 1.0.4      2017-08-22 CRAN (R 3.5.0)                 
    ##  openxlsx               4.1.4      2019-12-06 cran (@4.1.4)                  
    ##  parallel               3.5.1      2018-07-05 local                          
    ##  parallelMap            1.4        2019-05-17 cran (@1.4)                    
    ##  ParamHelpers           1.13       2019-12-02 cran (@1.13)                   
    ##  pcaPP                  1.9-73     2018-01-14 CRAN (R 3.5.0)                 
    ##  pheatmap               1.0.12     2019-01-04 CRAN (R 3.5.2)                 
    ##  pillar                 1.4.3      2019-12-20 cran (@1.4.3)                  
    ##  pkgconfig              2.0.3      2019-09-22 cran (@2.0.3)                  
    ##  plotly                 4.9.2      2020-02-12 cran (@4.9.2)                  
    ##  prettyunits            1.1.1      2020-01-24 cran (@1.1.1)                  
    ##  progress               1.2.2      2019-05-16 cran (@1.2.2)                  
    ##  promises               1.1.0      2019-10-04 cran (@1.1.0)                  
    ##  purrr                  0.3.3      2019-10-18 cran (@0.3.3)                  
    ##  R6                     2.4.1      2019-11-12 cran (@2.4.1)                  
    ##  RColorBrewer           1.1-2      2014-12-07 CRAN (R 3.5.0)                 
    ##  Rcpp                   1.0.4      2020-03-17 cran (@1.0.4)                  
    ##  readr                  1.3.1      2018-12-21 cran (@1.3.1)                  
    ##  readxl                 1.3.1      2019-03-13 cran (@1.3.1)                  
    ##  reprex                 0.3.0      2019-05-16 cran (@0.3.0)                  
    ##  Rgraphviz              2.24.0     2018-05-01 Bioconductor                   
    ##  rio                    0.5.16     2018-11-26 cran (@0.5.16)                 
    ##  rlang                  0.4.5      2020-03-01 cran (@0.4.5)                  
    ##  rmarkdown              2.1        2020-01-20 cran (@2.1)                    
    ##  robustbase             0.93-2     2018-07-27 CRAN (R 3.5.0)                 
    ##  roxygen2               6.1.0      2018-07-27 CRAN (R 3.5.0)                 
    ##  rrcov                  1.5-2      2020-01-16 CRAN (R 3.5.2)                 
    ##  Rtsne                  0.15       2018-11-10 cran (@0.15)                   
    ##  rvest                  0.3.5      2019-11-08 cran (@0.3.5)                  
    ##  S4Vectors              0.18.3     2018-06-08 Bioconductor                   
    ##  scales                 1.1.0      2019-11-18 cran (@1.1.0)                  
    ##  shiny                  1.4.0.2    2020-03-13 cran (@1.4.0.2)                
    ##  slam                   0.1-47     2019-12-21 cran (@0.1-47)                 
    ##  splines                3.5.1      2018-07-05 local                          
    ##  stats                * 3.5.1      2018-07-05 local                          
    ##  stats4                 3.5.1      2018-07-05 local                          
    ##  stringi                1.4.6      2020-02-17 cran (@1.4.6)                  
    ##  stringr                1.4.0      2019-02-10 cran (@1.4.0)                  
    ##  survival               2.42-6     2018-07-13 CRAN (R 3.5.0)                 
    ##  survminer              0.4.6      2019-09-03 cran (@0.4.6)                  
    ##  survMisc               0.5.5      2018-07-05 CRAN (R 3.5.0)                 
    ##  tibble                 2.1.3      2019-06-06 CRAN (R 3.5.2)                 
    ##  tidyr                  1.0.2      2020-01-24 cran (@1.0.2)                  
    ##  tidyselect             1.0.0      2020-01-27 cran (@1.0.0)                  
    ##  tidyverse              1.3.0      2019-11-21 cran (@1.3.0)                  
    ##  tools                  3.5.1      2018-07-05 local                          
    ##  tsne                   0.1-3      2016-07-15 cran (@0.1-3)                  
    ##  tweenr                 1.0.1      2018-12-14 CRAN (R 3.5.0)                 
    ##  utils                * 3.5.1      2018-07-05 local                          
    ##  vctrs                  0.2.4      2020-03-10 cran (@0.2.4)                  
    ##  vipor                  0.4.5      2017-03-22 CRAN (R 3.5.0)                 
    ##  viridis                0.5.1      2018-03-29 CRAN (R 3.5.0)                 
    ##  viridisLite            0.3.0      2018-02-01 CRAN (R 3.5.0)                 
    ##  withr                  2.1.2      2018-03-15 CRAN (R 3.5.0)                 
    ##  xfun                   0.12       2020-01-13 cran (@0.12)                   
    ##  XML                    3.99-0.3   2020-01-20 cran (@3.99-0.)                
    ##  xml2                   1.2.5      2020-03-11 cran (@1.2.5)                  
    ##  xtable                 1.8-4      2019-04-21 cran (@1.8-4)                  
    ##  yaml                   2.2.1      2020-02-01 cran (@2.2.1)                  
    ##  zip                    2.0.4      2019-09-01 cran (@2.0.4)                  
    ##  zlibbioc               1.26.0     2018-05-01 Bioconductor                   
    ##  zoo                    1.8-7      2020-01-10 cran (@1.8-7)
