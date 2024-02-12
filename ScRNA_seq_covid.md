Single cell RNA-Seq analysis
================
Kousalya Devi Murugesan
2024-02-12

### Outline of workflow followed in the analysis

1.  Data download

2.  Create count matrix and create Seurat objects

3.  Merged Seurat Objects

4.  Quality control and filtering

5.  Normalize and Scale data

6.  Find neighbors and clusters

7.  Identification of conserved markers

8.  Assigning cell-types based on both manual and SingleR annotation

9.  Perform pseudobulk differential expression analysis using DESeq2

10. Identifying blood transfusion modules over-represented for
    significantly differentialy expressed genes

11. Visualization of expression of some marker genes

### Data download

The CITE-Seq data was downloaded from Gene Expression Omnibus (GEO)
under accession number GSE155673

### Load libraries

``` r
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 'SeuratObject' was built under R 4.3.1 but the current version is
    ## 4.3.2; it is recomended that you reinstall 'SeuratObject' as the ABI
    ## for R may have changed

    ## 'SeuratObject' was built with package 'Matrix' 1.6.3 but the current
    ## version is 1.6.5; it is recomended that you reinstall 'SeuratObject' as
    ## the ABI for 'Matrix' may have changed

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following object is masked from 'package:base':
    ## 
    ##     intersect

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

### Set working directory

``` r
setwd("~/Documents/Postdoc/PSA_lab/Rawdata/")
```

### Get data location

``` r
dirs <- list.dirs(path = '~/Documents/Postdoc/PSA_lab/Rawdata', recursive = F, full.names = F)
```

### Create the count matrix and Seurat object for each sample

``` r
for (i in dirs){
  # create variable name
  name <- gsub('GSE155673_','',i)
  
  # Create count matrix
  cts <- ReadMtx(mtx = paste0(i, '/matrix.mtx'), 
                features = paste0(i, '/features.tsv'), 
                cells = paste0(i, '/barcodes.tsv'))
  
  # Create seurat object
  assign(name, CreateSeuratObject(counts = cts)) 
  }
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

### Add sample characteristics to each Seurat object

``` r
library(readxl)
covid_19 <- read_xlsx("../covid-19_scrna_sample_characteristics.xlsx")

S_obj = c(cov01, cov02, cov03, cov04, cov07, cov08, cov09, cov10, cov11, cov12, cov17, cov18)

cov01 <- AddMetaData(cov01, metadata = c(covid_19[1,]), col.name = colnames(covid_19))
cov02 <- AddMetaData(cov02, metadata = c(covid_19[2,]), col.name = colnames(covid_19))
cov03 <- AddMetaData(cov03, metadata = c(covid_19[3,]), col.name = colnames(covid_19))
cov04 <- AddMetaData(cov04, metadata = c(covid_19[4,]), col.name = colnames(covid_19))
cov07 <- AddMetaData(cov07, metadata = c(covid_19[5,]), col.name = colnames(covid_19))
cov08 <- AddMetaData(cov08, metadata = c(covid_19[6,]), col.name = colnames(covid_19))
cov09 <- AddMetaData(cov09, metadata = c(covid_19[7,]), col.name = colnames(covid_19))
cov10 <- AddMetaData(cov10, metadata = c(covid_19[8,]), col.name = colnames(covid_19))
cov11 <- AddMetaData(cov11, metadata = c(covid_19[9,]), col.name = colnames(covid_19))
cov12 <- AddMetaData(cov12, metadata = c(covid_19[10,]), col.name = colnames(covid_19))
cov17 <- AddMetaData(cov17, metadata = c(covid_19[11,]), col.name = colnames(covid_19))
cov18 <- AddMetaData(cov18, metadata = c(covid_19[12,]), col.name = colnames(covid_19))
```

### Merge Seurat objects

``` r
merged_seurat <- merge(cov01, y = c(cov02, cov03, cov04, cov07, cov08, cov09, cov10, cov11, cov12, cov17, cov18), add.cell.ids = ls()[1:12], project = "Covid_19")
```

### Add information about the individual and the condition to the metadata

``` r
## Create a sample column
merged_seurat$cells <- rownames(merged_seurat@meta.data)

## Split the sample column 
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'cells', into = c('sample_id', 'barcode'), sep = '_')
```

### Standard pre-processing workflow

#### Calculate percentage reads mapped to mitochondrial genome

``` r
#Ccalculate percent mitochondrial reads
merged_seurat[['percent.mt']] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
```

#### Determining the complexity of cells

The complexity of cells were determined by calculating number of genes
per UMI for each cell. This was calculated to get an idea of possible
doublets.

``` r
# Add number of genes per UMI for each cell to metadata

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
```

## Understanding the raw data through vizualization

#### QC and selecting cells for further analysis

#### Quality control plots

##### Number of cell counts per sample

``` r
# Visualize the number of cell counts per sample

ggplot(merged_seurat@meta.data, aes(x = Sample_id)) + 
  geom_bar() + ggtitle("Ncells") + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) + 
  geom_text(stat = 'count', aes(label=..count..), vjust = -1)
```

    ## Warning: The dot-dot notation (`..count..`) was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `after_stat(count)` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

#### Visualize the distribution of mitochondrial gene expression detected per cell

``` r
ggplot(merged_seurat@meta.data, aes(x=percent.mt, fill=Sample_id)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    geom_vline(xintercept = 10)
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI

``` r
ggplot(merged_seurat@meta.data, aes(x=log10GenesPerUMI, color = Sample_id, fill=Sample_id)) +
    geom_density(alpha = 0.2) + 
    geom_vline(xintercept = 0.8) +
    ggtitle("Distribution of genes detected per UMI") + 
    theme(plot.title = element_text(hjust=0.5, face="bold"))
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

#### Visualize the distribution of number of UMIs per cell

``` r
## Density plot
ggplot(merged_seurat@meta.data, (aes(x=nCount_RNA, fill = Disease_status))) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("log10 Cell density") +
    geom_vline(xintercept = 1000) +
    ggtitle("The distribution of number of UMIs per cell before filtering") + 
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    facet_wrap(~Sample_id)
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
## Violin plot
ggplot(merged_seurat@meta.data, (aes(x=Sample_id, y= nCount_RNA, fill = Disease_status))) + 
  theme_classic() +
  scale_y_log10() +
  ylab("log10nCount_RNA") +
    geom_violin() + 
    geom_boxplot(width=0.1) +
    stat_summary(fun = "mean",
               geom = "point",
               aes(color = "Mean")) +
  stat_summary(fun = "median",
               geom = "point",
               aes(color = "Median")) +
  scale_colour_manual(values = c("red", "blue"),name = "") +
    geom_hline(yintercept = 500) +
    ggtitle("The distribution of number of UMIs per cell before filtering") + 
    theme(plot.title = element_text(hjust=0.5, face="bold"))
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

#### Visualize the distribution of genes detected per cell

``` r
#Violin plot
ggplot(merged_seurat@meta.data, (aes(x=Sample_id, y= nFeature_RNA, fill = Disease_status))) + 
  theme_classic() +
  scale_y_log10() +
  ylab("log10nFeature_RNA") +
    geom_violin() + 
    geom_boxplot(width=0.1) +
    stat_summary(fun = "mean",
               geom = "point",
               aes(color = "Mean")) +
  stat_summary(fun = "median",
               geom = "point",
               aes(color = "Median")) +
  scale_colour_manual(values = c("red", "blue"),name = "") +
    geom_hline(yintercept = 300) +
    ggtitle("The distribution of number of genes per cell before filtering") + 
    theme(plot.title = element_text(hjust=0.5, face="bold"))
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# Boxplot
ggplot(merged_seurat@meta.data, aes(x=Sample_id, y=nFeature_RNA), fill = Sample_id) +        geom_boxplot() + 
  ggtitle("Distribution of genes in each sample") + 
  theme(plot.title = element_text(hjust=0.5, face="bold"))
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

#### Visualize the correlation between genes detected and number of UMIs

``` r
ggplot(merged_seurat@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color = percent.mt)) + geom_point() + 
scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 50) +  
scale_x_log10() +
scale_y_log10() +
geom_vline(xintercept = 1000) +
geom_hline(yintercept = 300) +
facet_wrap(~Sample_id)
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### Filtering

- Remove empty droplets with very few no. of genes

- Remove doublets with very high no. of genes

- Remove low quality or dead cells with high mitochondrial contamination

  ### Thresholds used for filtering

  - Immune cells are highly complex and have low mitochondrial ratio

  - Percentage of Mitochondria \< 10 (To remove dead and damaged cells)

  - log10 Genes per UMI \> 0.8 (To remove lysed cells and possible
    doublets(less complex cells)

  - No. of UMI per cell \> 500

  - No. of genes per cell \> 300

``` r
merged_seurat_filtered<- subset(merged_seurat, subset = percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA > 300 & log10GenesPerUMI >0.8)
```

### Re-assess QC metrics

``` r
# Visualize the distribution of number of UMIs per cell
# ggplot(merged_seurat_filtered@meta.data, (aes(x=nCount_RNA, fill = Disease_status))) + 
#     geom_density(alpha = 0.2) + 
#     theme_classic() +
#     ggtitle("The distribution of number of UMIs per cell after filtering") + 
#     theme(plot.title = element_text(hjust=0.5, face="bold")) +
#     scale_x_log10() + 
#     ylab("log10 Cell density") +
#     geom_vline(xintercept = 1000) + geom_vline(intercept = 500) +
#     facet_wrap(~Sample_id)

## Violin plot
ggplot(merged_seurat_filtered@meta.data, (aes(x=Sample_id, y= nCount_RNA, fill = Disease_status))) + 
  theme_classic() +
  scale_y_log10() +
  ylab("log10nCount_RNA") +
    geom_violin() + 
    geom_boxplot(width=0.1) +
    stat_summary(fun = "mean",
               geom = "point",
               aes(color = "Mean")) +
  stat_summary(fun = "median",
               geom = "point",
               aes(color = "Median")) +
  scale_colour_manual(values = c("red", "blue"),name = "") +
    geom_hline(yintercept = 500) +
    ggtitle("The distribution of number of UMIs per cell after filtering") + 
    theme(plot.title = element_text(hjust=0.5, face="bold"))
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

#### Visualize the distribution of genes detected per cell via histogram

``` r
# Density plot
# ggplot(merged_seurat_filtered@meta.data, aes(x=nFeature_RNA, fill = Disease_status)) + 
#   geom_density(alpha=0.2) + 
#   scale_x_log10() + 
#   geom_vline(xintercept = 300) +
#   ggtitle("Distribution of genes after filtering") + theme(plot.title = element_text(hjust=0.5,       face="bold")) +
#   facet_wrap(~Sample_id)

# Violin plot
ggplot(merged_seurat_filtered@meta.data, (aes(x=Sample_id, y= nFeature_RNA, fill = Disease_status))) + theme_classic() + scale_y_log10() + ylab("log10nFeature_RNA") +
    geom_violin() + 
    geom_boxplot(width=0.1) +
    stat_summary(fun = "mean",
               geom = "point",
               aes(color = "Mean")) +
  stat_summary(fun = "median",
               geom = "point",
               aes(color = "Median")) +
  scale_colour_manual(values = c("red", "blue"),name = "") +
    geom_hline(yintercept = 300) +
    ggtitle("The distribution of number of genes per cell after filtering") + 
    theme(plot.title = element_text(hjust=0.5, face="bold"))
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

#### Visualize the correlation between genes detected and number of UMIs

``` r
# Colored by percentage of mitochondria
ggplot(merged_seurat_filtered@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color = percent.mt)) + geom_point() + 
scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 50) +  
scale_x_log10() +
scale_y_log10() +
geom_vline(xintercept = 1000) +
geom_hline(yintercept = 300) +
ggtitle("Correlation between number of UMIs and genes detected after filtering") + 
theme(plot.title = element_text(hjust=0.5, face="bold")) +
facet_wrap(~Sample_id)
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
# Color by log10GenesPerUMI (i.e. complexity of cells)
ggplot(merged_seurat_filtered@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color = log10GenesPerUMI)) + geom_point() + 
scale_color_gradient() +  
scale_x_log10() +
scale_y_log10() +
geom_vline(xintercept = 1000) +
geom_hline(yintercept = 300) +
facet_wrap(~Sample_id)
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
rm(merged_seurat_filtered)
```

### Normalizing the data

``` r
# m_seurat_filtered_no10 <- NormalizeData(m_seurat_filtered_no10, normalization.method = "LogNormalize", scale.factor = 10000)
```

### Highly variable features (features selection)

``` r
# m_seurat_filtered_no10 <- FindVariableFeatures(m_seurat_filtered_no10)
```

### Scaling the data

``` r
# For scaling all genes
# all.genes <- rownames(merged_seurat_filtered)
# merged_seurat_filtered <-  ScaleData(merged_seurat_filtered, features = all.genes)
# 
# # For scaling only the highly variable features identified in previous step
# m_seurat_filtered_no10 <- ScaleData(m_seurat_filtered_no10)
```

### Run PCA

``` r
merged_seurat_filtered <- readRDS("../output/m_seurat_filtered.rds")
DimPlot(merged_seurat_filtered, reduction = "pca", group.by = "Disease_status")
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
rm(merged_seurat_filtered)
```

The sample Cov_10 clustered separately. Since the sample also had a
different distribution compared to other samples and the patient had
severe covid and deceased, this sample was removed from further
analysis.

### PCA after removing the Cov_10 sample

``` r
m_seurat_filtered_no10 <- readRDS("../output/m_seurat_filtered_no10_clus.rds")
DimPlot(m_seurat_filtered_no10, reduction = "pca", group.by = "Disease_status")
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

### Elbow plot to identify PCs for downstream analysis

``` r
ElbowPlot(m_seurat_filtered_no10, ndims = 40)
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

Dimensions upto 25 explained most of the variance. The dimesion 25 was
used in the paper and the same was used.

## Find neighbours and clusters of unintegrated data

In the PCA plot the samples from both diseased and healthy individuals
over-lapped well. Hence batch correction and integration of data was not
performed. The paper also mentioned the absence of batch effects.

Resolution 0.5 was used to cluster the cells by UMAP

``` r
# m_seurat_filtered_no10 <- FindNeighbors(m_seurat_filtered_no10, dims = 1:25, reduction = "pca")
# 
# m_seurat_filtered_no10 <- FindClusters(m_seurat_filtered_no10, resolution = 0.5, cluster.name = "unintegrated_clusters")
```

### RunUMAP to visualize the clusters

``` r
# m_seurat_filtered_no10 <- RunUMAP(m_seurat_filtered_no10, dims = 1:25, reduction = "pca", reduction.name = "umap.unintegrated")

# DimPlot(m_seurat_filtered_no10, reduction = "umap.unintegrated", label = TRUE, label.size = 6)
```

### Exploration of quality control metrics

``` r
## Distribution of cells per cluster in each sample
## Extract cluster and sample information from seurat object to determine number of cells per cluster per sample
n_cells <- FetchData(m_seurat_filtered_no10, 
            vars = c("seurat_clusters", "Disease_status")) %>% dplyr::count(seurat_clusters, Disease_status) %>% tidyr::spread(seurat_clusters, n)

View(n_cells)

#Visualize cells per cluster using UMAP
DimPlot(m_seurat_filtered_no10, label = F, split.by = "Disease_status", group.by = "Disease_status")
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

### Segregation of clusters by various sources of un-interesting variation

``` r
metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI")

FeaturePlot(m_seurat_filtered_no10, reduction = "umap.unintegrated", features = metrics, pt.size = 0.4, min.cutoff = 'q10', label = TRUE)
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

### Identification of conserved markers in all condition

``` r
# DefaultAssay(m_seurat_filtered_no10) <-"RNA"
# 
# #Test it out in one cluster "0"
# m_seurat_filtered_no10 <- JoinLayers(m_seurat_filtered_no10)
# cluster0_cons_markers <- FindConservedMarkers(m_seurat_filtered_no10, 
#                          ident.1 = 0, 
#                          grouping.var = "Disease_status", 
#                          only.pos = FALSE,
#                          logfc.threshold = 0.25)
# 
# # Create function to get conserved markers for any given cluster
# get_conserved <- function(cluster){
#   FindConservedMarkers(m_seurat_filtered_no10,
#                        ident.1 = cluster,
#                        grouping.var = "Disease_status",
#                        only.pos = TRUE) %>%
#     rownames_to_column(var = "gene") %>%
#     left_join(y = unique(annotations[, c("gene_name", "description")]),
#                by = c("gene" = "gene_name")) %>%
#     cbind(cluster_id = cluster, .)
# }
# 
# # Iterate function across desired clusters
# conserved_markers <- map_dfr(c(0:24), get_conserved)
```

#### Evaluating marker genes and manual annotation of clusters

- Extract top 10 markers per cluster based on average FC

- Extract top 10 markers (with avg.FC. \> 2) based on the metric
  (pct.1 - pct.2) that differentiates the cluster from other clusters

The significant markers thus obtained were used to manually annotate
cell-type of the clusters using “The Human Protein Atlas” database
(<https://www.proteinatlas.org/>).

``` r
# # Extract top 10 markers per cluster
# top10_conserved_markers_pct <- conserved_markers %>% 
#   mutate(avg_fc = (covid_19_avg_log2FC + Healthy_avg_log2FC) /2) %>% 
#   mutate(covid_19_pct = (covid_19_pct.1-covid_19_pct.2)) %>%
#   mutate(Healthy_pct = (Healthy_pct.1-Healthy_pct.2)) %>%
#   group_by(cluster_id) %>% 
#   top_n(n = 10, 
#         wt = avg_fc)
# 
# # Add avg_fc, difference between pct.2-pct.1
# conserved_markers_full_unf <- conserved_markers %>% 
#   mutate(avg_fc = (covid_19_avg_log2FC + Healthy_avg_log2FC) /2) %>% 
#   mutate(covid_19_pct = (covid_19_pct.1-covid_19_pct.2)) %>%
#   mutate(Healthy_pct = (Healthy_pct.1-Healthy_pct.2)) %>%
#   group_by(cluster_id)
# 
# # Filter rows with avg_fc greater than or equal to 2
# conserved_markers_full_fil_avgfc <- dplyr::filter(conserved_markers_full_unf, `avg_fc` >= 2)
# 
# #order the table by covid_19 pct in desc
# conserved_markers_full_fil_avgfc <- conserved_markers_full_fil_avgfc %>% group_by(cluster_id) %>% arrange(desc(covid_19_pct), .by_group = TRUE) 
# 
# # Extract top 10 genes with highest pct and log2fc >= 2 
# conserved_markers_top10_covid_pct <- conserved_markers_full_fil_avgfc %>% group_by(cluster_id) %>% top_n(n=10, wt = covid_19_pct)
# 
# # Visualize top 10 markers per cluster
# View(top10_conserved_markers_pct)
# 
# #Order the rows by difference in pct
# 
# top10_conserved_markers_pct <- top10_conserved_markers_pct %>% group_by(cluster_id) %>% arrange(covid_19_pct, .by_group = TRUE)
# 
# #Save_files
# write.csv(top10_conserved_markers_pct, 
#           file = "../output/top10_conserved_markers_pct.csv", 
#           quote = FALSE, 
#           row.names = FALSE)
# 
# write.csv(conserved_markers, 
#           file = "../output/conserved_markers_full_list.csv", 
#           quote = FALSE, 
#           row.names = FALSE)
# 
# write.csv(conserved_markers_top10_covid_pct, 
#           file = "../output/conserved_markers_top10_pct.csv", 
#           quote = FALSE, 
#           row.names = FALSE)
# 
# write.csv(conserved_markers_full_fil_avgfc, 
#           file = "../output/conserved_markers_log2fc_morethan_2.csv", 
#           quote = FALSE, 
#           row.names = FALSE)
```

### Using singleR to annotate the cells

SingleR was used to annotate cells based on “HumanPrimaryCellAtlasData”
database.

``` r
# # Load reference library
# library(SingleR)
# library(celldex)
# 
# # Load reference data with Ensembl annotations
# ref <- celldex::HumanPrimaryCellAtlasData(ensembl = TRUE)
# View(as.data.frame(colData(ref)))
# assays(ref)$logcounts
# 
# # Extract raw counts from filtered seurat object
# library(SingleCellExperiment)
# m_seurat_filtered_no10_sce <- as.SingleCellExperiment(m_seurat_filtered_no10)
# 
# require(EnsDb.Hsapiens.v110)
# ens <- mapIds(ens110,
#   keys = rownames(m_seurat_filtered_no10_sce),
#   column = 'GENEID',
#   keytype = 'SYMBOL')
# all(rownames(m_seurat_filtered_no10_sce) == names(ens))
# 
# keep <- !is.na(ens)
# ens <- ens[keep]
# m_seurat_filtered_no10_sce <- m_seurat_filtered_no10_sce[keep,]
# rownames(m_seurat_filtered_no10_sce) <- ens
# 
# #Run singleR
# results_singleR <- SingleR(test = m_seurat_filtered_no10_sce, ref = ref, labels = ref$label.main ) 
# 
# # Copy the singleR labels to Seurat object
# 
# m_seurat_filtered_no10$singleR.labels <- results_singleR$labels[match(rownames(m_seurat_filtered_no10@meta.data), rownames(results_singleR))]
# 
# m_seurat_filtered_no10$singleR.finelabels <- results_singleR$labels[match(rownames(m_seurat_filtered_no10@meta.data), rownames(results_singleR))]
# 
# DimPlot(m_seurat_filtered_no10, reduction = 'umap.unintegrated', group.by = 'singleR.finelabels', )
# #
# 
# # summarise the number of cell types in each cluster
# # Identify the cell-type of each cluster based on the maximum cell-type(95%)
# 
# clus_singleR_uniq <- m_seurat_filtered_no10@meta.data %>% group_by(seurat_clusters, singleR.finelabels) %>% summarise(n())
# 
# write_csv(clus_singleR_uniq, path = "../output/singleR_summarised_output.csv")
```

### Rename all identities with name of cell-type

The cell-types were annotated based on both manual and SingleR
annotation

``` r
# current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
# new.cluster.ids <- c("NK",
# "CD4 Memory",
# "CD8",
# "CD4 Memory or Naïve",
# "C_MONO_1",
# "C_MONO_2",
# "B",
# "C_MONO_3",
# "Platelets",
# "MONO_1",
# "GMP/Plasmacytoid_DC",
# "NC_MONO_1",
# "NC_MONO_2",
# "CD4",
# "Myeloid_DC_1",
# "Plasma cell",
# "CMP",
# "CD8",
# "Cell cycle",
# "Myeloid_DC_2",
# "B/T/NK",
# "NK/T","NK/CD8","NC_MONO_3","MONO_2")
# 
# m_seurat_filtered_no10@meta.data$seurat_clusters <- plyr::mapvalues(x = m_seurat_filtered_no10@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
```

#### Visualization of annotated clusters

``` r
DimPlot(object = m_seurat_filtered_no10, 
        reduction = "umap.unintegrated", 
        group.by = "seurat_clusters",
        label = TRUE,
        label.size = 3.5,
        repel = TRUE) + theme_classic() + ggtitle("Annotated clusters") + theme(plot.title = element_text(hjust=0.5, face="bold")) 
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

### Download annotation

``` r
# library(AnnotationHub)
# #Connect to annnotation hub
# ah <- AnnotationHub()
# 
# #Access the Ensemble database for organism
# ahdb <- query(ah, pattern =c("Homo sapiens", "EnsDb"), ignore.case=TRUE)
# 
# # Acquire the latest annotation files
# id <- ahdb %>% mcols() %>% rownames() %>% tail(n=1)
# 
# #Download the appropriate Ensembl database
# edb <- ah[id]
# 
# # View edb and find "AH113665"
# edb
# 
# # Retrieve the record for ensbsb
# ens110 <- ah[["AH113665"]]
# 
# # Extract gene-level information from database
# annotations <- genes(ens110, return.type = "data.frame")
# 
# # select annotations of interest
# annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
```

### Differential expression - pseudo bulk RNA-Seq analysis

Though Seurat has functions to calculate differential expression, the
p-values from the functions are often inflated because it considers the
each cell as independent sample. However, this assumption is not true as
cells that belong to one individual are not independent (as they have
same genetics and exposed to same environment). Hence, the better
approach proposed is to aggregate counts (i.e. pseudo bulk) of cells for
each sample and perform differential expression using softwares that
perform bulk RNA-Seq analysis.

In this analysis sample counts were aggregated using Seurat package and
differential expression analysis was performed using DESeq2 package.
Initially, differential expression was calculated by contrasting
Covid-19 affected vs. Healthy individuals. The significantly
differential expressed genes were then used to perform
over-representation analysis using tmod package.

The blood transcription modules over represented by the significantly
differential expressed genes are presented below.

![](../Plots/BTM_output.png)

The modules such as: Innate antiviral response, Antigen processing and
presentation, Interferon, genes enriched in B,NK and dendritic cells
were over-represented in significantly up-regulated genes. The modules
associated with translation and protein synthesis were over-represented
in down-regulated genes. This corroborated with the modules reported in
the paper.

However, the analysis could not confirm the differential expression of
the marker for severe Covid-19 “S100A12” . It was significantly
up-regulated only in Plasma cell type. The other markers “OSM”,
“TNFSF14” were not significantly differentially expressed between
diseased and healthy condition.

Hence, I re-ran the analysis by contrasting the conditions as : Severe
vs. Healthy; Severe vs. Moderate and Moderate vs. Healthy.

The figures of BTMs identified are given below,

1.  Severe vs. Healthy. The immune related pathways were
    over-represented in Severe covid patients vs. Healthy individuals

![](../Plots/BTM_sev_healthy.png)

2.  Moderate vs. Healthy: Most of immune related modules were
    over-represented in up=regulated genes like the Severe vs. Healthy
    comparison

    ![](../Plots/BTM_mod_vs_h.png)

3.  Severe vs. Moderate: Immune related pathways related to interferon
    response like type-I INF response, interferon and antiviral IFN
    signature were down-regulated. This was in confirmation with studies
    (Hadjadj et al., 2020, smith et.al 2022) that reported impaired type
    I interferon response in severe/hospitalized Covid patients.

    ![](../Plots/BTM_sev_vs_Mod.png)

### Visualizing marker genes

### Plot gene markers severe covid

``` r
FeaturePlot(m_seurat_filtered_no10, reduction = "umap.unintegrated", 
                                    features = c("S100A12", "OSM", "TNFSF14"), 
                                    label = 'seurat_annotations', label.size = 0.5,            
                                    order = TRUE,            
                                    min.cutoff = 'q10', 
                                    repel = FALSE, 
                                    split.by = 'Disease_severity') 
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

``` r
VlnPlot(object = m_seurat_filtered_no10, features = c("S100A12", "OSM", "TNFSF14", "AGER", "IL6", "HLA-DRA"), split.by = 'Disease_severity', group.by = 'Sample_id')
```

    ## The default behaviour of split.by has changed.
    ## Separate violin plots are now plotted side-by-side.
    ## To restore the old behaviour of a single split violin,
    ## set split.plot = TRUE.
    ##       
    ## This message will be shown once per session.

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-36-2.png)<!-- -->

### Plotting ISGs

``` r
FeaturePlot(m_seurat_filtered_no10, reduction = "umap.unintegrated", 
                                    features =c("IFITM3", "ISG15","IFI27"), 
                                    label = T, label.size = 0.5,            
                                    order = TRUE,            
                                    min.cutoff = 'q10', repel = FALSE, 
                                    split.by = 'Disease_severity')
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
VlnPlot(object = m_seurat_filtered_no10, features = c("IFITM3", "ISG15","IFI27"), group.by = 'Disease_severity')
```

![](ScRNA_seq_covid_files/figure-gfm/unnamed-chunk-37-2.png)<!-- -->
