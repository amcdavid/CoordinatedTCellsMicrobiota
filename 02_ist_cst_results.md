IST/CST Results
================

``` r
#Figure 3 IST Composition Heatmaps

library("readr")
library("tidyverse")
```

    ## ── Attaching packages ───────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.0.0     ✔ purrr   0.2.5
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
    ## ✔ tidyr   0.8.1     ✔ stringr 1.3.1
    ## ✔ ggplot2 3.0.0     ✔ forcats 0.3.0

    ## Warning: package 'tibble' was built under R version 3.5.2

    ## Warning: package 'dplyr' was built under R version 3.5.2

    ## ── Conflicts ──────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library("RColorBrewer")
library("pheatmap")
```

    ## Warning: package 'pheatmap' was built under R version 3.5.2

``` r
library("viridis")
```

    ## Loading required package: viridisLite

``` r
library("ggplot2")
library("igraph")
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     compose, simplify

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     crossing

    ## The following object is masked from 'package:tibble':
    ## 
    ##     as_data_frame

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library("slam")
```

    ## Warning: package 'slam' was built under R version 3.5.2

``` r
library("scales")
```

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:viridis':
    ## 
    ##     viridis_pal

    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard

    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

``` r
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

ics <- read_tsv(file.path(pth, "ics_basic.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   Subject = col_character(),
    ##   Visit = col_character(),
    ##   IST = col_character(),
    ##   Renamed_IST = col_character(),
    ##   Sex = col_character(),
    ##   Cohort = col_character(),
    ##   MOD = col_character(),
    ##   GAB = col_double(),
    ##   BirthCohort = col_character(),
    ##   DOL = col_integer(),
    ##   CGA = col_double()
    ## )

``` r
tphe <- read_tsv(file.path(pth, "tphe_basic.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   Subject = col_character(),
    ##   Visit = col_character(),
    ##   IST = col_character(),
    ##   Renamed_IST = col_character(),
    ##   Sex = col_character(),
    ##   Cohort = col_character(),
    ##   MOD = col_character(),
    ##   GAB = col_double(),
    ##   BirthCohort = col_character(),
    ##   DOL = col_integer(),
    ##   CGA = col_double()
    ## )

``` r
ggplot(ics, aes(y=CGA, x=Renamed_IST, color=GAB )) + scale_color_gradient2(midpoint=37, low="red", mid="blue", high="darkblue", space ="Lab" ) + geom_quasirandom() + coord_flip()
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

These were reordered by average PMA in the version of this figure in the
paper.

``` r
ggplot(tphe, aes(y=CGA, x=Renamed_IST, color=GAB )) + scale_color_gradient2(midpoint=37, low="red", mid="blue", high="darkblue", space ="Lab" ) + geom_quasirandom() + coord_flip()
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

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
    ##   DOL = col_integer(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_integer(),
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
    ##   DOL = col_integer(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_integer(),
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
    ##   Taxon = col_character(),
    ##   B160QC3F.03 = col_integer(),
    ##   E160RC6F.03 = col_integer(),
    ##   H160RC7D.03 = col_integer(),
    ##   G160R60V.03 = col_integer(),
    ##   E160W873.03 = col_integer(),
    ##   K160WBN8.03 = col_integer()
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
    ##   Taxon = col_character(),
    ##   G160QNB9.04 = col_integer(),
    ##   A160PNVS.04 = col_integer(),
    ##   B160NVWP.04 = col_integer(),
    ##   C160PHZT.04 = col_integer(),
    ##   C160PQRM.04 = col_integer(),
    ##   J160R3VV.04 = col_integer(),
    ##   J160RC5Z.04 = col_integer()
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

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
pheatmap(mat = mat.nas, color = colors, annotation_col = anno.nas, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = cumsum(unname(table(anno.nas[[CLUSTER_COLUMN]]))))
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
#Figure 4 CST Occurence Over PMA

library(readr)
library(ggplot2)
library(tidyverse)
library(ggbeeswarm)

rec <- read_tsv(file.path(pth, "rec_basic.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   DOL = col_integer(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_integer(),
    ##   CST = col_character(),
    ##   Renamed_CST = col_character(),
    ##   PreviousCST = col_character(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character()
    ## )

``` r
nas <- read_tsv(file.path(pth, "nas_basic.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   DOL = col_integer(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_integer(),
    ##   CST = col_character(),
    ##   Renamed_CST = col_character(),
    ##   PreviousCST = col_character(),
    ##   PostInitialDischarge = col_character(),
    ##   IllnessVisit = col_character()
    ## )

``` r
ggplot(rec, aes(y=CGA, x=Renamed_CST, color=gaBirth)) + scale_color_gradient2(midpoint=37, low="red", mid="blue", high="darkblue", space ="Lab" ) + geom_quasirandom() + coord_flip()
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggplot(nas, aes(y=CGA, x=Renamed_CST, color=gaBirth)) + scale_color_gradient2(midpoint=37, low="red", mid="blue", high="darkblue", space ="Lab" ) + geom_quasirandom() + coord_flip()
```

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

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
    ##   DOL = col_integer(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_integer(),
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

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
nas <- read_tsv(file.path(pth, "nas_basic.txt"))
```

    ## Parsed with column specification:
    ## cols(
    ##   SampleID = col_character(),
    ##   ID = col_character(),
    ##   Subject = col_character(),
    ##   DOL = col_integer(),
    ##   MOD = col_character(),
    ##   Sex = col_character(),
    ##   gaBirth = col_double(),
    ##   CGA = col_double(),
    ##   Reads = col_integer(),
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

![](02_ist_cst_results_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

\#Supplementary Figure 5 gCST 3 Time to Occurrence Based on Tphe IST at
Discharge

This works exactly the same as the survival analysis script in
`02_network_modeling`, except you can specify that it run for only a
subset of CSTs and immune variables, and for all of those that are
significant a plot is generated. Optionally, you can save the fitted
models and filtered input data for the chosen CST-immune subset.

See 02\_network\_modeling “time to occurrence”/survival analysis script
for additional comments and description.

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
library(survminer)
```

    ## Warning: package 'survminer' was built under R version 3.5.2

    ## Loading required package: ggpubr

    ## Warning: package 'ggpubr' was built under R version 3.5.2

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     set_names

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

``` r
library(dplyr)
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:slam':
    ## 
    ##     rollup

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## The following object is masked from 'package:purrr':
    ## 
    ##     transpose

``` r
library(readr)

type <- "alt_surv"

sites = c("REC")

#Included for fig making
sigCSTs <- c("REC_3")

#Included for fig making
sigVars <- c("TPHE_IST_Disch")

#Included for fig making
fits.alt_surv <- list()

#Included for fig making
table.alt_surv <- list()

for (site in sites) {
  
  mapping <- read_tsv(file.path(refined, sprintf("%s_Surv_Mapping.txt", site)))
  
  raw_table <- read_tsv(file.path(refined, sprintf("%s_Surv_Input.txt", site)))
  
  mapping_length <- ncol(mapping)
  
  table_length <- ncol(raw_table)
  
  for (cst_i in 3:table_length) {
    
    tmp.time <- Sys.time()
    
    cst <- colnames(raw_table[,cst_i])
    
    #Included for fig making 
    if (!(cst %in% sigCSTs)) { next ; }
    
    #Included for fig making
    fits.alt_surv[[cst]] <- list()
    
    #Included for fig making
    table.alt_surv[[cst]] <- list()
    
    cst_in <- data.frame(raw_table[, 1:2], raw_table[, cst_i])
    
    cst_in <- cst_in[complete.cases(cst_in[,]),]
    
    var_names <- list()
    voi_names <- list()
    var_coeffs <- list()
    var_exp_coeffs <- list()
    variable_pvals <- list()
    k <- 1
    
    for (var_i in 4:mapping_length) {
      
      var <- colnames(mapping[,var_i])
      
      #Included for fig making  
      if (!(var %in% sigVars)) { next ; }
      
      mapping_in <- data.frame(mapping[, 1:3], mapping[, var_i])
      
      fac = FALSE
      
      working_table <- merge(mapping_in, cst_in, by = "Subject")
      
      working_table <- working_table[complete.cases(working_table[,]), ]
      
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
      
      if ((nrow(working_table) < 20) || ((fac) && (nrow(table(working_table[[var]])) < 2))) {
        
        next
        
      }
      
      working_table$cst_obs <- working_table[[cst]]
      
      if (fac) {
        
        working_table$voi <- factor(working_table[[var]])
        
      } else {
        
        working_table$voi <- c(scale(working_table[[var]], center = TRUE, scale = TRUE))
        
      }
      
      working_table$GAB <- ((working_table$gaBirth)/37) - 1
      
      cst_surv <- Surv(working_table$PrevSampleDOL, working_table$cst_obs, type = "interval2")
      
      #Included for fig making
      fits.alt_surv[[cst]][[var]] <- ic_par(cst_surv ~ MOD + GAB + voi, data = working_table, model = "aft", dist = "loglogistic")
      
      #Included for fig making
      table.alt_surv[[cst]][[var]] <- working_table
      
      #Included for fig making
      next
      
      #Included for fig making
      if(!(dir.exists(fig_path <- file.path(figure, "cst_occurence_figs", type))))
      dir.create(fig_path)
      if (fac) {
        
        fac_lvls <- levels(as.data.frame(table(working_table[[var]]))[[1]])
        df.new <- data.frame(voi = fac_lvls, GAB = c(-.2, -.2), MOD = c("Caesarean_Section", "Caesarean_Section"))
        rownames(df.new) <- fac_lvls
      } else {
        df.new <- data.frame(voi = c(1, -1), GAB = c(-.2, -.2), MOD = c("Caesarean_Section", "Caesarean_Section"))
        rownames(df.new) <- c('High', 'Low')
      }
       pdf(file.path(fig_path, sprintf("%s_%s_plot.pdf", cst, var)))
        plot(fits.alt_surv[[cst]][[var]], df.new, xlab = "Day of Life")
        dev.off()
        
      
      #Included for fig making
      next
      
      fit <- ic_par(cst_surv ~ MOD + GAB + voi, data = working_table, model = "aft", dist = "loglogistic")
      
      fit_summary <- summary(fit)
      
      for (p in 5:nrow(fit_summary$summaryParameters)) {
        
        var_names[[k]] <- var
        voi_names[[k]] <- rownames(fit_summary$summaryParameters)[p]
        var_coeffs[[k]] <- fit_summary$summaryParameters[p, 1]
        var_exp_coeffs[[k]] <- fit_summary$summaryParameters[p, 2]
        variable_pvals[[k]] <- fit_summary$summaryParameters[p, 5]
        k <- k + 1
        
      }
      
      if(!(dir.exists(place <- file.path(refined, type))))
      dir.create(place)
      
      cat("Current Variable: ", file = file.path(place, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat(var, file = file.path(place, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat('\n', file = file.path(place, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      capture.output(summary(fit), file = file.path(place, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
      cat('\n###\n###\n###\n\n', file = file.path(place, sprintf("%s_results/%s_fit_log.txt", type, cst)), append = TRUE)
        
    }
    
    #Included for fig making
    next
    
    adj_variable_pvals <- p.adjust(variable_pvals, method = "fdr")
    results <- cbind(var_names, voi_names, var_coeffs, var_exp_coeffs, variable_pvals, adj_variable_pvals)
    
    write.csv(results, file.path(place, sprintf("%s_results/%s_model_pvals.csv", type, cst)))
    
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

``` r
#Included for fig making
#save(fits.alt_surv, file = "selected_fitted_alt_survival_models.rda")
#save(table.alt_surv, file = "selected_fitted_alt_surv_dfs.rda")

#Print total runtime
cat("\nTotal runtime:")
```

    ## 
    ## Total runtime:

``` r
cat(difftime(Sys.time(), start.time, units = "hours"))
```

    ## 0.0003581128

``` r
cat("\n\n\n")
```
