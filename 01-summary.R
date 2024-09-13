# Reading the data
library(broadSeq)
library(ggplot2)

se <- readRDS(system.file("extdata","rat_vole_mouseSE_salmon.rds", package = "broadSeq"))
SummarizedExperiment::assayNames(se)

#Sample Metadata
as.data.frame(colData(se)) %>% dplyr::count(stage,species) %>% tidyr::spread(stage,n)

#Filtering out low expression genes
assays(se)[["counts"]][,5] %>% ggpubr::ggdensity(y = "count")+
  ggplot2::geom_vline(xintercept = 10)+ggplot2::scale_x_log10()

keep <- (assays(se)[["counts"]] >= 3) %>% rowSums() >= 5
# smallest Group Size is 5
table(keep)

#Normalization
#CPM
se <- broadSeq::normalizeEdgerCPM(se ,method = "none",cpm.log = TRUE )
## The normalized values are added with the assay name "logCPM"
SummarizedExperiment::assayNames(se)

#TMM
se <- broadSeq::normalizeEdgerCPM(se , method = "TMM", cpm.log = FALSE )
## The normalized values are added with the assay name "TMM"
SummarizedExperiment::assayNames(se)

#access
assays(se)[["counts"]][1:5,1:5]

assays(se)[["TMM"]][1:5,1:5]

assays(se)[["logCPM"]][1:5,1:5]

# Transformation
# VST
se <- broadSeq::transformDESeq2(se,method = "vst"  )

# Normalized counts transformation
se <- broadSeq::transformDESeq2(se, method = "normTransform"  )

#rlog
se <- broadSeq::transformDESeq2(se, method = "rlog")
SummarizedExperiment::assayNames(se)

# Comparision
p <- broadSeq::sampleAssay_plot(se[, se$species=="Mouse" ], 
                                assayName = "counts", fill = "stage", 
                                yscale = "log2")+ rremove("x.text")

p1 <- broadSeq::sampleAssay_plot(se[, se$species=="Mouse"], 
                                 assayName = "vst", fill = "stage")+ rremove("x.text")

p2 <- broadSeq::sampleAssay_plot(se[, se$species=="Mouse"], 
                                 assayName = "TMM", fill = "stage", 
                                 yscale = "log10")+ rremove("x.text")

p3 <- broadSeq::sampleAssay_plot(se[, se$species=="Mouse"], 
                                 assayName = "logCPM", fill = "stage")+ rremove("x.text")


ggarrange(p,p1,p2,p3, common.legend = TRUE, labels = c("A","B","C"))


if (requireNamespace("vsn", quietly = TRUE)) {
  library("vsn")
  x <- meanSdPlot(
    log2(assays(se[, se$species == "Rat"])[["counts"]]+1),
    plot = FALSE)
  print(x$gg +ggtitle(label = "log2(n+1) "))
  
  x <- meanSdPlot(
    assays(se[, se$species == "Rat"])[["vst"]],
    plot = FALSE)
  
  print(x$gg +ggtitle(label = "Vst"))
  
  x <- meanSdPlot(
    assays(se[, se$species == "Rat"])[["logCPM"]],
    plot = FALSE)
  print(x$gg + ggtitle(label = "logCPM"))
}    

# 5. Visualization of gene Expression
## Multiple assay of a single gene
broadSeq::assay_plot(se, feature = c("Shh"), 
                     assayNames = c("counts","logCPM","vst","TMM"),
                     x = "stage", fill="species", add="dotplot", palette = "npg")

## Expression of multiple genes from a single assay
broadSeq::genes_plot(se, 
                     features = c("Shh","Edar"), 
                     facet.by = "symbol",
                     x = "stage", assayName = "vst", fill="species", palette = "jco")

# Pre-defined or custom color palette based on journals
jco <- broadSeq::genes_plot(se[,se$species == "Mouse"], 
                            features = c("Shh"), facet.by = "symbol", assayName =  "logCPM",
                            x = "stage",  fill="stage", add="dotplot", xlab = "", 
                            title = "Journal of Clinical Oncology", palette = "jco") 

npg <- broadSeq::genes_plot(se[,se$species == "Mouse"], 
                            features = c("Shh"), facet.by = "symbol",assayName =  "logCPM",
                            x = "stage", fill="stage", add="dotplot", xlab = "",
                            title = "Nature Publishing Group", palette = "npg") 

aaas <- broadSeq::genes_plot(se[,se$species == "Mouse"], 
                             features = c("Shh"), facet.by = "symbol", assayName = "logCPM",
                             x = "stage", fill="stage", add="dotplot", xlab = "",
                             title = "Science", palette = "aaas")

nejm <- broadSeq::genes_plot(se[,se$species == "Mouse"], 
                             features = c("Shh"), facet.by = "symbol", assayName = "logCPM",
                             x = "stage", fill="stage", add="dotplot",  xlab = "",
                             title = "New England Journal of Medicine",palette = "nejm")


ggarrange(jco+ggpubr::rotate_x_text(), npg+ggpubr::rotate_x_text(),
          aaas+ggpubr::rotate_x_text(),nejm+ggpubr::rotate_x_text(),
          nrow = 1, common.legend = TRUE,legend = "none",
          labels = c("A","B","C","D")) %>% 
  annotate_figure( top = text_grob("Color palette")) 

# 6. QC with Clustering
# MDS Plot
broadSeq::plot_MDS(se, scaledAssay = "vst", ntop=500, 
                   color = "species", shape = "stage", 
                   ellipse=TRUE, legend = "bottom")
head(rowData(se))

# Hierarchical clustering and Heatmap
p_vst <- broadSeq::plotHeatmapCluster(
  se,
  scaledAssay = "vst",
  annotation_col = c("species", "stage"),
  annotation_row = c("Class","gene_biotype"),
  ntop = 30, show_geneAs = "symbol",
  cluster_cols = TRUE, cluster_rows = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  main = "Top 30 variable gene vst"
)

# PCA Plot 
#prcompTidy
computedPCA_logCPM <- broadSeq::prcompTidy(se, scaledAssay = "logCPM", ntop = 500)

computedPCA_vst <- broadSeq::prcompTidy(se, scaledAssay = "vst", ntop = 500)

#logCPM
plotAnyPC(computedPCA = computedPCA_logCPM,
          x = 1, y = 2, color = "species", shape = "stage",
          legend = "bottom")

#VST
pca_vst <- plotAnyPC(computedPCA = computedPCA_vst,
                     x = 2, y = 3,  color = "species", shape = "stage", 
                     legend = "bottom") 

#Other PCs
computedPCA_vst$eigen_values %>%
  dplyr::filter(var_exp >= 2) %>%
  ggbarplot(x="PC",y="var_exp", label = TRUE, label.pos = "out")

pca_vst_2_3 <-plotAnyPC(computedPCA = computedPCA_vst,
                        x = 2, y = 3,  
                        color = "species", shape = "stage", legend = "bottom")

#Gene Loading
computedPCA_vst %>% broadSeq::getFeatureLoadRanking(keep = c("symbol","Class")) %>% head()

#  Top 5 genes of PC2
computedPCA_vst$loadings %>% top_n(5,abs(PC2)  ) %>% dplyr::select(gene,PC2)

pca_vst_loading <- computedPCA_vst %>% 
  broadSeq::getFeatureLoadRanking(keep = c("symbol","Class"), topN = 50, pcs=1:10) %>% 
  dplyr::count(Class, PC) %>%
  ggbarplot(
    x = "PC", y = "n", fill = "Class",
    legend = "bottom", palette = c("red","blue","orange","purple","white","grey")
  ) 

#Biplot
# By default it plots top 2 genes from each PC axis
pca_vst_bi <- broadSeq::biplotAnyPC(computedPCA = computedPCA_vst, 
                                    x = 1, y = 2, genesLabel = "symbol", 
                                    color = "species", shape = "stage", 
                                    legend = "bottom")

ggarrange(
  ggarrange(pca_vst_bi+ggtitle(label =  ""),
            pca_vst_2_3+ggtitle(label =  ""), common.legend = TRUE),
  pca_vst_loading, nrow = 2)

# User defined genes : Now plotting top 5 genes from PC3
# Top 5 genes of PC3
biplotAnyPC(computedPCA = computedPCA_vst,x = 2, y = 3, 
            color = "species", shape = "stage",
            genes= computedPCA_vst$loadings %>% 
              top_n(5,abs(PC3)) %>% pull(gene),
            genesLabel = "symbol")

## Plot progression gene "Shh" 
biplotAnyPC(computedPCA = computedPCA_vst,x = 2, y = 3, 
            color = "species", shape = "stage",
            genes=c("Shh"),
            genesLabel = "symbol")

# 7. Compare Differential expression
#Data
se <- readRDS(system.file("extdata","rat_vole_mouseSE_salmon.rds", package = "broadSeq"))

# To reduce the run time, subset of the data used here
se <- se[,colData(se)$species == "Mouse"]

#Gene information
head(rownames(se))

head(rowData(se))

head(colData(se))

table(colData(se)$stage)

# Differential Expression
result_Noiseq <- 
  use_NOIseq(se = se, 
             colData_id = "stage", control = "Bud", treatment = "Cap",
             rank = TRUE, 
             r = 10) # r is an argument of NOISeq::noiseqbio
head(result_Noiseq)

pg <- broadSeq::genes_plot(se, x = "stage", assayName =  "counts",  
                           features = result_Noiseq %>% dplyr::filter(rank <5) %>% rownames(),
                           fill="stage", facet.by = "symbol",
                           palette="jco", add = "dotplot")+rotate_x_text()

pg_sc <- ggscatter(result_Noiseq, x="Bud_mean", y="Cap_mean",color = "prob")+ 
  scale_x_log10()+scale_y_log10()

pg+pg_sc


# limma 
?use_limma_trend(se, colData_id, control, treatment, rank = FALSE, ...)
?use_limma_voom(se, colData_id, control, treatment, rank = FALSE, ...)

# edgeR 
?use_edgeR_exact(se, colData_id, control, treatment, rank = FALSE, ...)
?use_edgeR_GLM(se, colData_id, control, treatment, rank = FALSE, ...)

# deseq2
?use_deseq2(se, colData_id, control, treatment, rank = FALSE, ...)

# DELocal
?use_DELocal(se, colData_id, control, treatment, rank = FALSE, ...)

# noiseq
?use_NOIseq(se, colData_id, control, treatment, rank = FALSE, ...) 

# EBSeq 
?use_EBSeq(se, colData_id, control, treatment, rank = FALSE, ...)

# samseq
?use_SAMseq(se, colData_id, control, treatment, rank = FALSE, ...)

# Compare DE results
# First define a named list of functions
funs <- list(limma_trend = use_limma_trend, limma_voom = use_limma_voom,
             edgeR_exact = use_edgeR_exact, edgeR_glm = use_edgeR_GLM,
             deseq2 = use_deseq2, 
             DELocal = use_DELocal, noiseq = use_NOIseq, 
             EBSeq = use_EBSeq) 


multi_result <- broadSeq::use_multDE(
  se = se, 
  deFun_list = funs, return.df = TRUE,  
  colData_id = "stage", control = "Bud", treatment = "Cap", 
  rank = TRUE)

head(multi_result)

# nrow(multi_result) == nrow(se)
colnames(multi_result)

#Similarity of methods
clusters <- multi_result %>% dplyr::select(ends_with("rank")) %>% t() %>% dist() %>% hclust()
plot(clusters,main =  "distance: Euclidean")

#Plots
#Volcano
multi_result %>% broadSeq::volcanoPlot(
  pValName = "deseq2_padj",
  lFCName = "deseq2_log2FoldChange",
  labelName = "symbol",
  palette = "lancet" ,
  selectedLabel =
    multi_result %>% dplyr::arrange(deseq2_padj) %>% pull(symbol) %>% head()
)

multi_result %>% broadSeq::volcanoPlot(
  pValName = "deseq2_padj",
  lFCName = "deseq2_log2FoldChange",
  labelName = "symbol",
  palette = c("purple","orange","grey"),
  selectedLabel = list(criteria = "(`x` > 5 | `x` < -2) & (`y` > 10)")
) +xlim(-7.5,7.5)  

sessionInfo()
  