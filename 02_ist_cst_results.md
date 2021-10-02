IST/CST Results
================

``` r
#Figure 3 IST Composition Heatmaps

library("readr")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("viridis")
library("ggplot2")
library("igraph")
library("slam")
library("scales")

#Name of IST/cluster column in metadata
tphe.cc <- "TPHE.IST"
ics.cc <- "ICS.IST"

#Read in sample metadata
pth = file.path(refined, 'dmn')
tphe.md <- read_delim(file.path(pth, "tphe_md.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   `TPHE IST` = col_character()
    ## )

``` r
ics.md <- read_delim(file.path(pth, "ics_md.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   `ICS IST` = col_character()
    ## )

``` r
rns <- tphe.md[[1]]
tphe.md <- tphe.md[ , 2:ncol(tphe.md)]
rownames(tphe.md) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
rns <- ics.md[[1]]
ics.md <- ics.md[ , 2:ncol(ics.md)]
rownames(ics.md) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
#Read in sample composition tables
tphe.cd4 <- read_delim(file.path(pth, "tphe_cd4_comp.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Metacluster = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
tphe.cd8 <- read_delim(file.path(pth, "tphe_cd8_comp.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Metacluster = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
ics.cd4 <- read_delim(file.path(pth, "ics_cd4_comp.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Metacluster = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
ics.cd8 <- read_delim(file.path(pth, "ics_cd8_comp.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Metacluster = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
rns <- tphe.cd4[[1]]
tphe.cd4 <- tphe.cd4[ , 2:ncol(tphe.cd4)]
rownames(tphe.cd4) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
rns <- tphe.cd8[[1]]
tphe.cd8 <- tphe.cd8[ , 2:ncol(tphe.cd8)]
rownames(tphe.cd8) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
rns <- ics.cd4[[1]]
ics.cd4 <- ics.cd4[ , 2:ncol(ics.cd4)]
rownames(ics.cd4) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
rns <- ics.cd8[[1]]
ics.cd8 <- ics.cd8[ , 2:ncol(ics.cd8)]
rownames(ics.cd8) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
#Match metadata and composition data
tphe.md <- data.frame(tphe.md)
ics.md <- data.frame(ics.md)
tphe.cd4 <- data.matrix(tphe.cd4)
tphe.cd8 <- data.matrix(tphe.cd8)
ics.cd4 <- data.matrix(ics.cd4)
ics.cd8 <- data.matrix(ics.cd8)
tphe.cd4 <- tphe.cd4[ , (colnames(tphe.cd4) %in% rownames(tphe.md))]
tphe.cd8 <- tphe.cd8[ , (colnames(tphe.cd8) %in% rownames(tphe.md))]
ics.cd4 <- ics.cd4[ , (colnames(ics.cd4) %in% rownames(ics.md))]
ics.cd8 <- ics.cd8[ , (colnames(ics.cd8) %in% rownames(ics.md))]

#Construct properly formatted annotations, color scheme, and composition matrices for use in pheatmap
tphe.cd4.anno <- tphe.md[colnames(tphe.cd4), tphe.cc, drop=FALSE]
tphe.cd8.anno <- tphe.md[colnames(tphe.cd8), tphe.cc, drop=FALSE]
ics.cd4.anno <- ics.md[colnames(ics.cd4), ics.cc, drop=FALSE]
ics.cd8.anno <- ics.md[colnames(ics.cd8), ics.cc, drop=FALSE]
tphe.cd4.anno[[tphe.cc]] <- factor(tphe.cd4.anno[[tphe.cc]])
tphe.cd8.anno[[tphe.cc]] <- factor(tphe.cd8.anno[[tphe.cc]])
ics.cd4.anno[[ics.cc]] <- factor(ics.cd4.anno[[ics.cc]])
ics.cd8.anno[[ics.cc]] <- factor(ics.cd8.anno[[ics.cc]])
tphe.cd4.anno <- tphe.cd4.anno[order(tphe.cd4.anno[[tphe.cc]]), , drop=FALSE]
tphe.cd8.anno <- tphe.cd8.anno[order(tphe.cd8.anno[[tphe.cc]]), , drop=FALSE]
ics.cd4.anno <- ics.cd4.anno[order(ics.cd4.anno[[ics.cc]]), , drop=FALSE]
ics.cd8.anno <- ics.cd8.anno[order(ics.cd8.anno[[ics.cc]]), , drop=FALSE]
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100)
tphe.cd4.mat <- tphe.cd4[, rownames(tphe.cd4.anno), drop=FALSE]
tphe.cd8.mat <- tphe.cd8[, rownames(tphe.cd8.anno), drop=FALSE]
ics.cd4.mat <- ics.cd4[, rownames(ics.cd4.anno), drop=FALSE]
ics.cd8.mat <- ics.cd8[, rownames(ics.cd8.anno), drop=FALSE]
tphe.cd4.mat <- t(apply(tphe.cd4.mat, 1L, scales::rescale))
tphe.cd8.mat <- t(apply(tphe.cd8.mat, 1L, scales::rescale))
ics.cd4.mat <- t(apply(ics.cd4.mat, 1L, scales::rescale))
ics.cd8.mat <- t(apply(ics.cd8.mat, 1L, scales::rescale))

#Make the heatmaps
pheatmap(mat = tphe.cd4.mat, color = colors, annotation_col = tphe.cd4.anno, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = cumsum(unname(table(tphe.cd4.anno[[tphe.cc]]))))
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
pheatmap(mat = ics.cd4.mat, color = colors, annotation_col = ics.cd4.anno, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = cumsum(unname(table(ics.cd4.anno[[ics.cc]]))))
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
pheatmap(mat = tphe.cd8.mat, color = colors, annotation_col = tphe.cd8.anno, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = cumsum(unname(table(tphe.cd8.anno[[tphe.cc]]))))
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
pheatmap(mat = ics.cd8.mat, color = colors, annotation_col = ics.cd8.anno, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = cumsum(unname(table(ics.cd8.anno[[ics.cc]]))))
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->

``` r
#Figure 3 IST Occurence Over PMA

library(readr)
library(ggplot2)
library(tidyverse)
library(ggbeeswarm)




# To get CGA at timepoints
timeline = read_csv(file.path('data', 'subject_timeline.csv')) %>% mutate(SampleID = str_c(Subject, '_',  `Sequence Num`))
```

    ## Parsed with column specification:
    ## cols(
    ##   `Sequence Num` = col_double(),
    ##   DOL = col_double(),
    ##   cga = col_double(),
    ##   Subject = col_character()
    ## )

``` r
subject = read_csv(file.path('data', 'subject_covariates.csv')) %>% mutate(GAB = 37 - preterm_weeks, BirthCohort = ifelse(preterm_weeks > 0, 'Pre-term', 'Full-term'))
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
timeline = left_join(timeline, subject)
```

    ## Joining, by = "Subject"

``` r
# recode by mean CGA
recode_ist = function(tab){
  s = tab %>% group_by(IST) %>% summarize(mcga = mean(CGA)) %>% arrange(mcga)
  s = s %>% mutate(oIST = fct_reorder(factor(str_c('oIST_', seq_len(nrow(.)))), mcga), oIST_num = as.numeric(oIST))
  left_join(tab, s)
}

# ICS by Subject and timepoint
pth = file.path(refined, 'dmn')
tphe <- read_delim(file.path(pth, "tphe_md.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600) %>% 
  mutate(IST = gsub('TPHE IST ([0-9])', 'TPHE_\\1', `TPHE IST`)) %>% left_join(timeline)%>% select(Subject, CGA = cga, IST, GAB, Visit = `Sequence Num`, BirthCohort) %>% recode_ist
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   `TPHE IST` = col_character()
    ## )

    ## Joining, by = "SampleID"

    ## Joining, by = "IST"

``` r
ics <- read_delim(file.path(pth, "ics_md.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)%>% 
  mutate(IST = gsub('ICS IST ([0-9])', 'ICS_\\1', `ICS IST`)) %>% left_join(timeline) %>% select(Subject, CGA = cga, IST, GAB,Visit = `Sequence Num`, BirthCohort) %>% recode_ist
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   `ICS IST` = col_character()
    ## )

    ## Joining, by = "SampleID"

    ## Joining, by = "IST"

``` r
# ics <- read_tsv(file.path(pth, "ics_basic.txt")) %>% recode_ist()
# tphe <- read_tsv(file.path(pth, "tphe_basic.txt")) %>% recode_ist()


ggplot(ics, aes(y=CGA, x=oIST, color=GAB )) + scale_color_gradient2(midpoint=37, low="red", mid="blue", high="darkblue", space ="Lab" ) + geom_quasirandom() + coord_flip()
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

These were reordered by average PMA in the version of this figure in the
paper.

``` r
summary(lm(CGA ~ IST, data = ics))
```

    ## 
    ## Call:
    ## lm(formula = CGA ~ IST, data = ics)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -53.997  -7.202  -0.354   2.285  58.074 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   34.022      2.133  15.948  < 2e-16 ***
    ## ISTICS_2       2.684      2.901   0.925   0.3554    
    ## ISTICS_3       5.264      3.001   1.754   0.0801 .  
    ## ISTICS_4       5.332      2.901   1.838   0.0668 .  
    ## ISTICS_5      12.821      3.090   4.150 4.08e-05 ***
    ## ISTICS_6      14.015      2.767   5.066 6.25e-07 ***
    ## ISTICS_7      28.148      3.525   7.984 1.55e-14 ***
    ## ISTICS_8      60.115      2.705  22.221  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 14.31 on 396 degrees of freedom
    ## Multiple R-squared:  0.6883, Adjusted R-squared:  0.6828 
    ## F-statistic: 124.9 on 7 and 396 DF,  p-value: < 2.2e-16

``` r
ggplot(tphe, aes(y=CGA, x=oIST, color=GAB )) + scale_color_gradient2(midpoint=37, low="red", mid="blue", high="darkblue", space ="Lab" ) + geom_quasirandom() + coord_flip()
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
summary(lm(CGA ~ IST, data = tphe))
```

    ## 
    ## Call:
    ## lm(formula = CGA ~ IST, data = tphe)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -42.376  -3.280   0.068   1.880  49.829 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   33.411      1.150  29.047  < 2e-16 ***
    ## ISTTPHE_2      4.872      1.611   3.024 0.002654 ** 
    ## ISTTPHE_3      6.075      1.701   3.571 0.000398 ***
    ## ISTTPHE_4      6.096      1.928   3.162 0.001686 ** 
    ## ISTTPHE_5     13.184      1.861   7.085 6.16e-12 ***
    ## ISTTPHE_6     45.535      2.247  20.268  < 2e-16 ***
    ## ISTTPHE_7     61.853      1.611  38.389  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 10.03 on 407 degrees of freedom
    ## Multiple R-squared:  0.847,  Adjusted R-squared:  0.8447 
    ## F-statistic: 375.5 on 6 and 407 DF,  p-value: < 2.2e-16

## Co-occurance of t cell measures

``` r
ics_tphe = bind_rows(ICS = ics, TPHE = tphe, .id = 'assay')
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector, coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector, coercing into character vector

``` r
assays_avail = ics_tphe %>% group_by(Subject, Visit, BirthCohort) %>% summarize(assays = str_c(assay, collapse = '_')) %>% mutate(assays = factor(assays, levels = c('ICS', 'TPHE', 'ICS_TPHE')))

counts_by_subj = assays_avail %>% group_by(Subject, assays, BirthCohort) %>% summarize(n = n()) %>% ungroup() %>% arrange(Subject, desc(n))
# Take modal scenario if a subject had different assays available at different time points (uncommon)
counts_by_subj = counts_by_subj[!duplicated(counts_by_subj$Subject),]
assays_by_term = with(counts_by_subj, table(n, BirthCohort, assays))
ftab = ftable(assays_by_term, row.vars = c('n', 'BirthCohort'))
ftab
```

    ##               assays ICS TPHE ICS_TPHE
    ## n BirthCohort                         
    ## 1 Full-term            0    2        5
    ##   Pre-term             1    2       10
    ## 2 Full-term            5    0       36
    ##   Pre-term             2    4       47
    ## 3 Full-term            0    2       37
    ##   Pre-term             1    1       30

``` r
write.ftable(ftab, file.path(refined, 'assay_consort_alternative.txt'))
```

Number of subjects with 1, 2 or 3 samples of the various assays,
stratified by Term.

``` r
ics_tphe = ics_tphe %>% mutate(Subjectf = fct_reorder(factor(Subject), GAB))
traj_plot = ggplot(ics_tphe, aes(y = Subjectf, x = CGA, fill = oIST_num)) + 
  geom_point(pch = 22) + scale_fill_distiller('IST', palette = 'GnBu') + facet_wrap(~assay) + 
  theme_minimal() + scale_y_discrete(breaks = NULL) + ylab("Subjects") + xlab('PMA') + geom_text(aes(label = oIST_num), size = 1.5) + theme(legend.position = 'bottom')
trajs = ics_tphe %>% group_by(assay) %>% do(plot = {
 out = traj_plot %+% .
  print(out)
  out
})
```

![](02_ist_cst_results_files/figure-gfm/t_ist_traj-1.png)<!-- -->![](02_ist_cst_results_files/figure-gfm/t_ist_traj-2.png)<!-- -->

``` r
#Figure 4 CST Composition Heatmap

library("readr")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("viridis")
library("ggplot2")
library("igraph")
library("slam")
library("scales")

CLUSTER_COLUMN <- "Renamed_CST"

md.rec <- read_delim(file.path(pth, "rec_basic.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   DOL = col_double(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_double(),
    ##   CST = col_character(),
    ##   Renamed_CST = col_character(),
    ##   PreviousCST = col_character(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character()
    ## )

``` r
rns <- md.rec[[1]]
md.rec <- md.rec[ , 2:ncol(md.rec)]
rownames(md.rec) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
md.nas <- read_delim(file.path(pth, "nas_basic.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   DOL = col_double(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_double(),
    ##   CST = col_character(),
    ##   Renamed_CST = col_character(),
    ##   PreviousCST = col_character(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character()
    ## )

``` r
rns <- md.nas[[1]]
md.nas <- md.nas[ , 2:ncol(md.nas)]
rownames(md.nas) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
genera.rec <- read_delim(file.path(refined, "REC_top_taxa.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Taxon = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
rns <- genera.rec[[1]]
genera.rec <- genera.rec[ , 2:ncol(genera.rec)]
rownames(genera.rec) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
genera.nas <- read_delim(file.path(refined, "NAS_top_taxa.txt"), "\t", escape_double = FALSE, trim_ws = TRUE, guess_max = 3600)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Taxon = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
rns <- genera.nas[[1]]
genera.nas <- genera.nas[ , 2:ncol(genera.nas)]
rownames(genera.nas) <- rns
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
md.rec <- data.frame(md.rec)
md.nas <- data.frame(md.nas)

genera.rec <- data.matrix(genera.rec)
genera.nas <- data.matrix(genera.nas)

genera.rec <- genera.rec[ , (colnames(genera.rec) %in% rownames(md.rec))]
genera.nas <- genera.nas[ , (colnames(genera.nas) %in% rownames(md.nas))]

anno.rec <- md.rec[colnames(genera.rec), CLUSTER_COLUMN, drop=FALSE]
anno.nas <- md.nas[colnames(genera.nas), CLUSTER_COLUMN, drop=FALSE]

cga.rec <- md.rec[colnames(genera.rec), c(CLUSTER_COLUMN, "CGA"), drop=FALSE]
cga.nas <- md.nas[colnames(genera.nas), c(CLUSTER_COLUMN, "CGA"), drop=FALSE]

anno.rec[[CLUSTER_COLUMN]] <- fct_reorder(cga.rec[[CLUSTER_COLUMN]], cga.rec$CGA, mean)
anno.nas[[CLUSTER_COLUMN]] <- fct_reorder(cga.nas[[CLUSTER_COLUMN]], cga.nas$CGA, mean)

anno.rec <- anno.rec[order(anno.rec[[CLUSTER_COLUMN]]), , drop=FALSE]
anno.nas <- anno.nas[order(anno.nas[[CLUSTER_COLUMN]]), , drop=FALSE]

top25.rec <- head(names(rev(sort(rowSums(genera.rec)))), 25)
top25.nas <- head(names(rev(sort(rowSums(genera.nas)))), 25)

colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100)

mat.rec <- genera.rec[top25.rec, rownames(anno.rec), drop=FALSE]
mat.nas <- genera.nas[top25.nas, rownames(anno.nas), drop=FALSE]

mat.rec <- t(apply(mat.rec, 1L, scales::rescale))
mat.nas <- t(apply(mat.nas, 1L, scales::rescale))

pheatmap(mat = mat.rec, color = colors, annotation_col = anno.rec, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = cumsum(unname(table(anno.rec[[CLUSTER_COLUMN]]))))
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
pheatmap(mat = mat.nas, color = colors, annotation_col = anno.nas, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = cumsum(unname(table(anno.nas[[CLUSTER_COLUMN]]))))
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
#Figure 4 CST Occurence Over PMA

library(readr)
library(ggplot2)
library(tidyverse)
library(ggbeeswarm)


recode_cst = function(tab){
  s = tab %>% group_by(CST) %>% summarize(mcga = mean(CGA)) %>% arrange(mcga)
  s = s %>% mutate(oCST = fct_reorder(factor(str_c('oCST_', seq_len(nrow(.)))), mcga), oCST_num = as.numeric(oCST))
  left_join(tab, s)
}

rec <- read_tsv(file.path(pth, "rec_basic.txt")) %>% recode_cst
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   DOL = col_double(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_double(),
    ##   CST = col_character(),
    ##   Renamed_CST = col_character(),
    ##   PreviousCST = col_character(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character()
    ## )

    ## Joining, by = "CST"

``` r
nas <- read_tsv(file.path(pth, "nas_basic.txt")) %>% recode_cst
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   DOL = col_double(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_double(),
    ##   CST = col_character(),
    ##   Renamed_CST = col_character(),
    ##   PreviousCST = col_character(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character()
    ## )
    ## Joining, by = "CST"

``` r
ggplot(rec, aes(y=CGA, x=Renamed_CST, color=gaBirth)) + scale_color_gradient2(midpoint=37, low="red", mid="blue", high="darkblue", space ="Lab" ) + geom_quasirandom() + coord_flip()
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggplot(nas, aes(y=CGA, x=Renamed_CST, color=gaBirth)) + scale_color_gradient2(midpoint=37, low="red", mid="blue", high="darkblue", space ="Lab" ) + geom_quasirandom() + coord_flip()
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
communities = bind_rows(
  list(ics = ics, tphe = tphe), .id = 'community'
) %>% rename(Renamed_CST = IST) %>%
  bind_rows(bind_rows(list(nas = nas, rec = rec), .id = 'community')) %>% 
  mutate(Subject = factor(Subject), Renamed_CST = factor(Renamed_CST))
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector, coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector, coercing into character vector

``` r
hospital_humilk = read_csv('data/milk_hospital.csv') %>% rename(perinatal_milk = `Any Human Milk Perinatal`)
```

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   `Any Human Milk Perinatal` = col_logical()
    ## )

``` r
nabx_polish = read_csv('data/antibiotic_exposure.csv')
```

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   discharge = col_logical(),
    ##   group = col_character(),
    ##   `Number of systemic antibiotic` = col_double()
    ## )

``` r
covariates = read_csv("data/subject_covariates.csv")
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
timeline = read_csv("data/subject_timeline.csv") %>% filter()
```

    ## Parsed with column specification:
    ## cols(
    ##   `Sequence Num` = col_double(),
    ##   DOL = col_double(),
    ##   cga = col_double(),
    ##   Subject = col_character()
    ## )

``` r
nabx_time = nabx_polish %>% select(-group) %>% rename(n_antibiotics = `Number of systemic antibiotic`) %>% spread(discharge, n_antibiotics) %>% rename(n_antibiotics_discharge = 'TRUE', n_antibiotics_pre = 'FALSE')

discharge_humilk = read_csv( 'intermediates/milk_subject.csv')
```

    ## Parsed with column specification:
    ## cols(
    ##   Subject = col_character(),
    ##   months_surveyed = col_double(),
    ##   milk_months = col_double(),
    ##   any_milk_quarter = col_double()
    ## )

``` r
covariates =  purrr::reduce(list(covariates, hospital_humilk, nabx_time, discharge_humilk), left_join)
```

    ## Joining, by = "Subject"

    ## Joining, by = "Subject"
    ## Joining, by = "Subject"

``` r
covariates$n_antibiotics_pre = ifelse(covariates$preterm_weeks <=0, NA, covariates$n_antibiotics_pre)


has_cst = communities %>%
  group_by(Subject, Renamed_CST, .drop = FALSE) %>%
  summarize(has_cst = n()>0) %>% ungroup()

was_sampled = communities %>% group_by(community, Subject, .drop = FALSE) %>% 
  summarize(was_sampled = n() > 0) %>% 
  filter(was_sampled) %>% 
  ungroup()  %>%
  inner_join(communities %>% group_by(community, Renamed_CST) %>% summarize())
```

    ## Joining, by = "community"

``` r
nrow(has_cst)
```

    ## [1] 7585

``` r
has_cst = has_cst %>% semi_join(was_sampled)
```

    ## Joining, by = c("Subject", "Renamed_CST")

``` r
nrow(has_cst)
```

    ## [1] 6444

``` r
has_cst = has_cst %>% left_join(covariates, by = 'Subject') %>% 
  mutate(mode_delivery = fct_recode(mode_delivery, vaginal = c('Vaginal Vertex'), other = 'Vaginal Breech', other= 'Caesarean Section'),
         n_antibiotics_pre = ifelse(preterm_weeks<=0, 0, n_antibiotics_pre))
```

    ## Warning: Column `Subject` joining factor and character vector, coercing into character vector

``` r
cst_assoc = has_cst %>% group_by(Renamed_CST) %>% do({
    data = .
    full = glm(has_cst ~ preterm_weeks + mode_delivery + scale(birth_wt_gms) + scale(auc14) + preg_antibiotics + milk_months + perinatal_milk + n_antibiotics_pre + n_antibiotics_discharge, family = 'binomial', data = data)
    drop_term = update(full, . ~ . - preterm_weeks)
    full_aov = anova(drop_term, full, test = 'Chisq')
    term_only = update(full, . ~ preterm_weeks)
    drop_termonly = update(full, . ~ 1)
    term_only_aov = anova(term_only, drop_termonly, test = 'Chisq')
    tibble(full = list(full), drop_term = list(full_aov), term_only = list(term_only), drop_termonly = list(term_only_aov))
})
```

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: glm.fit: algorithm did not converge

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: glm.fit: algorithm did not converge

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

``` r
cst_assoc = pivot_longer(cst_assoc, -Renamed_CST, values_to = 'model') %>% 
  rowwise() %>% mutate(result = list(suppressWarnings(broom::tidy(model))))
cst_coef = unnest(cst_assoc %>% select(-model), cols = c(result))
cst_coef = cst_coef %>% group_by(is_intercept = term == "(Intercept)", name, is_term = term == 'preterm_weeks') %>% 
  mutate(fdr = p.adjust(p.value, method = 'fdr')) %>% ungroup()
```

## Associations between term and EVER CST

``` r
term_adj_results = filter(cst_coef,name == 'drop_term') %>% filter(!is.na(p.value)) %>% 
  select(Renamed_CST, p.value, fdr) %>% arrange(p.value)

term_noadj_results = filter(cst_coef,name == 'drop_termonly') %>% filter(!is.na(p.value)) %>% 
  select(Renamed_CST, p.value, fdr) 
knitr::kable(left_join(term_adj_results, term_noadj_results, by = "Renamed_CST", suffix = c('_controlled', '_simple')))
```

| Renamed\_CST | p.value\_controlled | fdr\_controlled | p.value\_simple | fdr\_simple |
|:-------------|--------------------:|----------------:|----------------:|------------:|
| REC\_11      |           0.0005169 |       0.0211912 |       0.3841019 |   0.5080058 |
| REC\_13      |           0.0055127 |       0.1130098 |       0.5854525 |   0.6858157 |
| NAS\_1       |           0.0097843 |       0.1337185 |       0.0000000 |   0.0000000 |
| ICS\_3       |           0.0146095 |       0.1497474 |       0.0000000 |   0.0000000 |
| NAS\_8       |           0.0207843 |       0.1704313 |       0.0002782 |   0.0009504 |
| NAS\_2       |           0.0360914 |       0.1830087 |       0.0000000 |   0.0000000 |
| REC\_6       |           0.0364280 |       0.1830087 |       0.0151703 |   0.0327359 |
| TPHE\_2      |           0.0369244 |       0.1830087 |       0.0000003 |   0.0000013 |
| REC\_3       |           0.0465117 |       0.1830087 |       0.0000000 |   0.0000001 |
| TPHE\_5      |           0.0495731 |       0.1830087 |       0.1978148 |   0.3119387 |
| ICS\_2       |           0.0514535 |       0.1830087 |       0.0386065 |   0.0753745 |
| ICS\_1       |           0.0535635 |       0.1830087 |       0.0038746 |   0.0117319 |
| NAS\_11      |           0.0790296 |       0.2492472 |       0.0605285 |   0.1078986 |
| ICS\_7       |           0.0985262 |       0.2885410 |       0.4931584 |   0.6127119 |
| REC\_5       |           0.1185858 |       0.3186147 |       0.6805234 |   0.7154221 |
| REC\_1       |           0.1243374 |       0.3186147 |       0.0106492 |   0.0256834 |
| NAS\_9       |           0.1494822 |       0.3605159 |       0.3166683 |   0.4477034 |
| TPHE\_4      |           0.2114238 |       0.4815765 |       0.0000029 |   0.0000133 |
| ICS\_8       |           0.2720756 |       0.5618162 |       0.2505320 |   0.3668504 |
| REC\_7       |           0.2740567 |       0.5618162 |       0.6126618 |   0.6868382 |
| ICS\_5       |           0.3255686 |       0.6356339 |       0.0080025 |   0.0205065 |
| NAS\_12      |           0.3516331 |       0.6468893 |       0.0001106 |   0.0004535 |
| TPHE\_6      |           0.3628891 |       0.6468893 |       0.4576067 |   0.5863086 |
| REC\_4       |           0.3911241 |       0.6469094 |       0.0433714 |   0.0808286 |
| NAS\_4       |           0.4101315 |       0.6469094 |       0.3667133 |   0.5011748 |
| NAS\_13      |           0.4102352 |       0.6469094 |       0.5233694 |   0.6311219 |
| ICS\_4       |           0.4442992 |       0.6746765 |       0.0000002 |   0.0000011 |
| ICS\_6       |           0.4891254 |       0.7162194 |       0.7765904 |   0.7765904 |
| REC\_10      |           0.5070235 |       0.7168263 |       0.0121089 |   0.0275814 |
| NAS\_10      |           0.5600719 |       0.7554856 |       0.7155748 |   0.7334642 |
| NAS\_6       |           0.5712208 |       0.7554856 |       0.1855594 |   0.3043174 |
| TPHE\_7      |           0.6214058 |       0.7961762 |       0.6198296 |   0.6868382 |
| NAS\_7       |           0.6627537 |       0.8234213 |       0.0049588 |   0.0135540 |
| TPHE\_1      |           0.7486164 |       0.9006041 |       0.0000000 |   0.0000000 |
| REC\_12      |           0.7688084 |       0.9006041 |       0.2327137 |   0.3533800 |
| REC\_8       |           0.8457613 |       0.9582787 |       0.0002631 |   0.0009504 |
| NAS\_5       |           0.8695595 |       0.9582787 |       0.6744430 |   0.7154221 |
| NAS\_3       |           0.9075440 |       0.9582787 |       0.1483642 |   0.2534555 |
| REC\_2       |           0.9115334 |       0.9582787 |       0.0000000 |   0.0000000 |
| REC\_9       |           0.9506361 |       0.9744020 |       0.0040060 |   0.0117319 |
| TPHE\_3      |           0.9905694 |       0.9905694 |       0.0314239 |   0.0644191 |

Controlling for preterm\_weeks + mode\_delivery + scale(birth\_wt\_gms)
+ scale(auc14) + preg\_antibiotics + milk\_months + perinatal\_milk +
n\_antibiotics\_pre + n\_antibiotics\_discharge

## Other associations between EVER CST (top 20)

``` r
other_assoc = filter(cst_coef,name == 'full', !is_intercept) %>% 
  select(Renamed_CST, term, estimate, std.error, fdr, p.value) %>% arrange(p.value) 
other_assoc %>% slice(1:20) %>% knitr::kable()
```

| Renamed\_CST | term                      |   estimate | std.error |       fdr |   p.value |
|:-------------|:--------------------------|-----------:|----------:|----------:|----------:|
| NAS\_4       | preg\_antibioticsYes      | -1.8602139 | 0.4804656 | 0.0286206 | 0.0001081 |
| TPHE\_2      | preg\_antibioticsYes      | -2.0715261 | 0.5519200 | 0.0286206 | 0.0001745 |
| REC\_11      | preterm\_weeks            |  0.4262421 | 0.1310133 | 0.0467492 | 0.0011402 |
| NAS\_4       | scale(auc14)              | -1.2719729 | 0.4043432 | 0.1811077 | 0.0016565 |
| ICS\_7       | n\_antibiotics\_discharge |  0.5313692 | 0.1736543 | 0.1815413 | 0.0022139 |
| REC\_6       | mode\_deliveryvaginal     |  1.1991511 | 0.4114048 | 0.2032382 | 0.0035595 |
| ICS\_2       | scale(birth\_wt\_gms)     | -2.2312997 | 0.7901970 | 0.2032382 | 0.0047469 |
| NAS\_9       | preg\_antibioticsYes      | -1.1992710 | 0.4309669 | 0.2032382 | 0.0053901 |
| TPHE\_5      | preg\_antibioticsYes      |  1.5007428 | 0.5428849 | 0.2032382 | 0.0057030 |
| REC\_1       | perinatal\_milkTRUE       | -1.9493134 | 0.7131852 | 0.2032382 | 0.0062713 |
| REC\_13      | scale(birth\_wt\_gms)     |  1.6446209 | 0.6106350 | 0.2032382 | 0.0070749 |
| REC\_11      | perinatal\_milkTRUE       | -1.6919789 | 0.6313174 | 0.2032382 | 0.0073605 |
| REC\_5       | mode\_deliveryvaginal     | -1.2373067 | 0.4622541 | 0.2032382 | 0.0074355 |
| REC\_13      | preterm\_weeks            |  0.3271683 | 0.1233737 | 0.1481942 | 0.0080052 |
| NAS\_1       | preterm\_weeks            |  0.5238912 | 0.2056328 | 0.1481942 | 0.0108435 |
| TPHE\_1      | mode\_deliveryvaginal     | -1.1112531 | 0.4421790 | 0.2828143 | 0.0119665 |
| REC\_11      | preg\_antibioticsYes      | -1.1818975 | 0.4708658 | 0.2828143 | 0.0120713 |
| REC\_1       | mode\_deliveryvaginal     |  1.1770991 | 0.4754450 | 0.2907043 | 0.0132944 |
| NAS\_8       | n\_antibiotics\_pre       |  0.0744229 | 0.0312557 | 0.3294952 | 0.0172612 |
| TPHE\_4      | mode\_deliveryvaginal     |  1.3895171 | 0.5919405 | 0.3294952 | 0.0189053 |

``` r
write_csv(other_assoc, 'intermediates/cst_assoc.csv')
```

Top 20 associations listed above, others [are
here](intermediates/cst_assoc.csv).

# Trajectories

``` r
rec_nas = bind_rows(REC = rec, NAS = nas, .id = 'assay') %>% rename(GAB = gaBirth)
rec_nas = rec_nas %>% mutate(Subjectf = fct_reorder(factor(Subject), GAB))
traj_plot = ggplot(rec_nas, aes(y = Subjectf, x = CGA, fill = oCST_num)) + 
  geom_point(pch = 22) + scale_fill_distiller('CST', palette = 'GnBu') + facet_wrap(~assay) + 
  theme_minimal() + scale_y_discrete(breaks = NULL) + ylab("Subjects") + xlab('PMA') + geom_text(aes(label = oCST_num), size = 1.5) + 
  theme(legend.position = 'bottom')
trajs2 = rec_nas %>% group_by(assay) %>% do(plot = {
  out = traj_plot %+% .
  print(out)
  out
})
```

![](02_ist_cst_results_files/figure-gfm/cst_traj-1.png)<!-- -->![](02_ist_cst_results_files/figure-gfm/cst_traj-2.png)<!-- -->

``` r
cowplot::plot_grid(plotlist = trajs$plot, ncol = 2)
```

![](02_ist_cst_results_files/figure-gfm/combined_traj-1.png)<!-- -->

``` r
cowplot::plot_grid(plotlist = trajs2$plot, ncol = 2)
```

![](02_ist_cst_results_files/figure-gfm/combined_traj-2.png)<!-- -->

\#Supplementary Figure 3 PCoA Plots

Runs in qiime.

``` sh


qiime diversity core-metrics-phylogenetic --i-phylogeny nas_rooted_tree.qza --i-table nas_table.qza --p-sampling-depth 1200 --m-metadata-file nas_basic.txt  --output-dir nas_cda --p-n-jobs 24

qiime diversity core-metrics-phylogenetic --i-phylogeny rec_rooted_tree.qza --i-table rec_table.qza --p-sampling-depth 2250 --m-metadata-file rec_basic.txt  --output-dir rec_cda --p-n-jobs 24
```

``` r
#Supplementary Figure 3 Axis Density by Term

library(readr)
library(ggplot2)

rec <- read_tsv(file.path(pth, "rec_basic.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   DOL = col_double(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_double(),
    ##   CST = col_character(),
    ##   Renamed_CST = col_character(),
    ##   PreviousCST = col_character(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character()
    ## )

``` r
rec$BirthCohort <- ifelse(rec$gaBirth >= 37, "Full term", "Preterm")

pc1 <- read_tsv(file.path('intermediates', "rec_pc1.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   PC1 = col_double()
    ## )

``` r
df <- merge(rec, pc1, by = "SampleID")

ggplot(data = df, aes(x = PC1, group = BirthCohort, fill = BirthCohort)) + geom_density(adjust=1.5, position="fill")
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
nas <- read_tsv(file.path(pth, "nas_basic.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   DOL = col_double(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_double(),
    ##   CST = col_character(),
    ##   Renamed_CST = col_character(),
    ##   PreviousCST = col_character(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character()
    ## )

``` r
nas$BirthCohort <- ifelse(nas$gaBirth >= 37, "Full term", "Preterm")

pc1 <- read_tsv(file.path('intermediates', "nas_pc1.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   PC1 = col_double()
    ## )

``` r
df <- merge(nas, pc1, by = "SampleID")

ggplot(data = df, aes(x = PC1, group = BirthCohort, fill = BirthCohort)) + geom_density(adjust=1.5, position="fill")
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->
