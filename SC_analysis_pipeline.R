# Analysis pipeline
source("Functions.R")
tenx_paths <- c("x", "y", "z")
tenx_files <- file.path("~/10X outputs", tenx_paths, "GRCh38")

#read in files
library(SoupX)
dataDirs <- tenx_files
ed_sx <- soupX_emptydrops(dataDirs, FDR_cut = 0.05)
sce <- ed_sx$SCE
sce <- do.call(cbind, sce)

#correct counts for mature kidney samples in the following way
scl <- ed_sx$Channels
scl = inferNonExpressedGenes(scl)
candidate_soup_genes <- lapply(names(scl$channels), function(x){
  return(add_gene_names(scl$channels[[x]]$nonExpressedGenes)[c(1:15), ])
})
ranked_soup_counts <- lapply(names(scl$channels), function(x){
  soup_profile <- scl$channels[[x]]$soupProfile
  soup_profile <- add_gene_names(soup_profile[order(soup_profile$cnts, decreasing = TRUE)[1:100], ])
  return(soup_profile)
}) 
Reduce(intersect, lapply(ranked_soup_counts, function(a){a$Gene}))
HB_genes <- c("HBB","HBD","HBG1","HBG2", "HBE1","HBZ","HBM","HBA2", "HBA1","HBQ1")
HB_genes <- rowData(sce)[match(HB_genes, rowData(sce)$Symbol), 1]
gene_use_list <- list(HB=HB_genes)
for(i in names(scl$channels)){
  scl = calculateContaminationFraction(scl,i, gene_use_list)
}
for(i in names(scl$channels)){
  scl = interpolateCellContamination(scl,i,useGlobal=FALSE)
}

scl = adjustCounts(scl, verbose = TRUE)

colnames(sce) <- colnames(scl$atoc)
counts(sce) <- scl$atoc
assay(sce, "original_counts") <- scl$toc
colData(sce)$rhos <- unlist(lapply(unique(names(scl$channels)) ,function(x){return(scl$channels[[x]]$rhos)}))
short_samp <- unlist(lapply(strsplit(colData(sce)$Sample,"/"), function(x){return(x[4])}))
colData(sce)$Short_Sample <- factor(short_samp, 
                                    levels = unique(short_samp) )
colnames(sce) <- paste0(sce$Short_Sample, "_", colData(sce)$Barcode)

sce$Experiment <- droplevels(sce$Experiment)

#### Quality control metrics

sce <- calculateQCMetrics(sce, feature_controls = list(
  "Mito" = grep("^MT-", rowData(sce)$Symbol),
  "HSP" = grep("^HSP", rowData(sce)$Symbol),
  "Ribo" = c(grep("^RPS", rowData(sce)$Symbol), 
             grep("^RPL", rowData(sce)$Symbol)),
  "HB" = grep("^HB", rowData(sce)$Symbol)),
  percent_top = c(10, 50, 100, 200, 500))

#initial cutoffs -  paper methods for compartment specific QC cutoffs used.
feature.drop <- sce$total_features_by_counts < 200
mit.drop <- sce$pct_counts_Mito > 70
master.drop <-  feature.drop | mit.drop 
sce <- sce[, master.drop == FALSE]

#### Normalisation
sce <- seurat_normalise(sce)
#####Feature selection
#x_low and x_high selected based on assessment of distribution of mean expression values.
x_low = 0.001 #x_low tune
x_high = 8 #x_high tune
hvg <- seurat_hvg(sce, x_low = x_low, x_high = x_high, y_cut = 1, num_bin = 20)
metadata(sce)$hvg <- hvg$HVG
#define genes that drive lots of technical variance.
hkGeneREGEX='^(DNAJ[ABC]|EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+',rowData(sce)$Symbol,value=TRUE), #Poorly characterised
                            grep('MALAT1',rowData(sce)$Symbol,value=TRUE), #Contamination or highly expressed poorly characterised
                            grep('^MT-',rowData(sce)$Symbol,value=TRUE), #Mitochondria
                            grep("XIST", rowData(sce)$Symbol, value = TRUE), #gender
                            grep(hkGeneREGEX,rowData(sce)$Symbol,value=TRUE) #Housekeeping genes
))
coreExcludeEnsemble <- rowData(sce)[match(coreExcludeGenes, rowData(sce)$Symbol), 1]

metadata(sce)$hvg <- metadata(sce)$hvg[!metadata(sce)$hvg %in% coreExcludeEnsemble]

#do MNN batch correction in PCA space.
set.seed(100)
original <- pblapply(unique(sce$Experiment), function(e){exprs(sce)[metadata(sce)$hvg, sce$Experiment == e]})
mnn <- do.call(fastMNN, c(original, pc.input = FALSE, auto.order=  TRUE)) #auto-order = TRUE
reducedDim(sce, "PCA") <- mnn$corrected

#or for the mature kidney, we have calculated PCA and dropped low PCs
set.seed(100)
sce <- runPCA(sce, ncomponents = 50, method = "irlba", exprs_values = "logcounts", feature_set = metadata(sce)$hvg, 
              scale_features = TRUE)
plot(attr(reducedDim(sce, "PCA"), "percentVar"), pch = 15, col = "dodgerblue", xlab = "PC", ylab = "percentVar")
n.PC = 10 #10 for initial embedding
abline(v = n.PC, col = "orange")
reducedDim(sce, "PCA") <- reducedDim(sce, "PCA")[, 1:n.PC]

#UMAP calculation
umap.config <- umap.defaults
umap.config$n_neighbors = 50 #50 for initial mature embedding 
umap.out <- umap(reducedDim(sce, "PCA"), config = umap.config, method = "umap-learn")
metadata(sce)$umap.out <- umap.out
reducedDim(sce, "UMAP")  <- umap.out$layout

#graph clustering based on the UMAP KNN graph
metadata(sce)$graph <- get_umap_graph(metadata(sce)$umap.out)
clusters <- cluster_SLM(metadata(sce)$graph) #resolution tune
sce$cluster <- as.factor(clusters)

##### Celltype similarity assessment

#calculation of joint PCA space between two samples
original <- pblapply(unique(sce$Tissue), function(t)
{sce_sub <- sce[, sce$Tissue == t]
exprs_in <- exprs(sce_sub)[metadata(sce_sub)$hvg,]
return(exprs_in)})
names(original) <- unique(sce$Tissue)
set.seed(100)
mbpca <- multiBatchPCA(original$tissue_1, original$tissue_2, approximate = TRUE)
names(mbpca@listData) <- names(original)
reducedDim(sce, "PCA") <- do.call(rbind, mbpca)

#similarity assessment using random forests  tool
#train is a reference SCE object, test is a query SCE object
features <- union(metadata(train)$hvg, metadata(test)$hvg)
rf_model <- ranger::ranger(celltype~., 
                           data = data.frame(t(as.matrix(exprs(train)[features, ])), 
                                             "celltype" = train$celltype),
                           probability = TRUE)
predictions <- predict(rf_model, data.frame(t(as.matrix(exprs(test)[features, ]))))
predictions <- data.frame(predictions$predictions, row.names = colnames(test))

#similarity assessment using knn regression tool
train = mbpca$tissue_1
train_celltype <- sce[, sce$tissue %in% "tissue_1"]$celltype
test = mbpca$tissue_2
test_celltype <-  sce[, sce$tissue %in% "tissue_2"]$celltype
k=20

#knn regression...
KNN_regression_output_list <- pblapply(unique(train_celltype), function(c){
  y = train_celltype == c
  y = as.numeric(y)
  KNN_classifier <- FNN::knn.reg(train =as.matrix(train), test =  as.matrix(test), 
                                 k=k, y = y)
  return(KNN_classifier$pred)
})
KNN_regression_output <- data.frame(do.call(cbind, KNN_regression_output_list))
colnames(KNN_regression_output) <- unique(train_celltype)

#aggregate by mean
enrich_data <- aggregate(KNN_regression_output, by =list(factor(test_celltype)), mean)
enrich_data <- data.frame(row.names = enrich_data[, 1], enrich_data[, -1])
pheatmap(enrich_data, color = viridis::magma(20), border_color = NA, cluster_rows = FALSE, cluster_cols = FALSE)

#### Geneset scoring
#as an example we use the M1/M2 score on a combined SCE of mature & fetal myeloiod found in Fig3, but this applies for all applications of this function
M1_M2 <- readRDS("Data/Genesets/M1 M2 markers.RDS")
M1_M2_scores <- seurat_test_genesets(sce, genelist = genesets)
M1_M2_scores <- scale(M1_M2_scores[sce$celltype %in% c("Adult_MNP1", "Adult_MNP2", "Adult_MNP4", "Fetus_Macrophage 1", "Fetus_Macrophage 2", "Fetus_Monocyte"), ])
means <- dplyr::group_by(data.frame(M1_M2_scores), factor(grps)) %>% dplyr::summarise_all(mean)
mean_scores <- t(data.frame(means[, 2:ncol(means)], row.names = means$`factor(grps)`))
pheatmap(mean_scores, color =viridis::cividis(10), cluster_cols = FALSE, cluster_rows = FALSE)

#### Marker genes
ave <- Matrix::rowMeans(exprs(sce))
markers <- tfidf_all_markers(sce[names(ave)[ave > 0.01], ], sce$cluster)
markers_use <- unique(unlist(lapply(markers, function(m){rownames(m)[1:5]})))
heat_in <- as.matrix(exprs(sce)[markers_use, order(sce$cluster)])
rownames(heat_in) <- rowData(sce)[rownames(heat_in), "Symbol"]
ggplot_heatmap(heat_in, scale_cap = 5, list("low" =  "blue",
                                            "mid" = "grey90",
                                            "high" = "red"))

#### GO term enrichment
ave <- calculateAverage(sce, exprs_values = "logcounts", use_size_factors = FALSE)
genes.use <-  names(ave)[ave >0.01]

#do this test comparing each cell type vs all other cell types
sce_wilcox_test <- pblapply(unique(sce$celltype), function(x){
  idx <- sce$celltype
  idx[!idx %in% x] <- "other"
  wilcox_test_out <- wilcox_test_sc(sce, idx, group1 = x, group2= "other", genes.use = genes.use)
  return(wilcox_test_out)
})
names(sce_wilcox_test) <- unique(sce$celltype)

#do the enrichments with clusterprofiler
library(clusterProfiler)
enrichments <- pblapply(names(sce_wilcox_test), function(x){
  test_set <- sce_wilcox_test[[x]]
  gene_vector <- test_set$gene[test_set$fdr.p.val < 0.05 & test_set$LFC > 0.5]
  universe <- test_set$gene
  enr <- clusterProfiler::enrichGO(gene_vector, 'org.Hs.eg.db', universe = universe, ont = "BP", keyType = "SYMBOL")
  enr <- simplify(enr)
  return(enr)
})
enrichments_results <- lapply(enrichments, function(x){x@result})
names(enrichments_results) <- names(sce_wilcox_test)

#plot the enrichment heatmap
enrichr_output <- lapply(enrichments_results, function(x){
  x <- x[1:10, ]
  return(x)})

terms <- unique(unlist(pblapply(enrichr_output, function(e){return(e$Description)})))
#term_p_value_map
enrichment_p_map <- matrix(NA, nrow = length(terms), ncol = length(enrichr_output))
rownames(enrichment_p_map) <- terms
colnames(enrichment_p_map) <- 1:length(enrichr_output)
for(i in seq_along(enrichr_output)){
  for(j in terms){
    en <- enrichr_output[[i]]
    idx <- en$Description %in% j
    if(all(!idx)){p = NA}else{
      p <- en[en$Description %in% j, "qvalue"]
      p <- -log(p)}
    enrichment_p_map[j, i] <- p
  }
}
colnames(enrichment_p_map) <- names(enrichr_output)
#we will use magma here...
pheatmap(enrichment_p_map, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE,
         color = viridis::plasma(5), labels_row = terms, fontsize_row = 7)

##### Slingshot pseudotime
#We first look for features that vary with stage
require(gam)
t <- as.numeric(sce$stages)
#calculate average expression
ave <- calculateAverage(sce, exprs_values = "logcounts", use_size_factors = FALSE)
features_use <- names(ave)[ave > 0.01]
hkGeneREGEX='^(DNAJ[ABC]|EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN|SLC[0-9])'
coreExcludeGenes = unique(c(grep('\\.[0-9]+',rowData(sce)$Symbol,value=TRUE), #Poorly characterised
                            grep('MALAT1',rowData(sce)$Symbol,value=TRUE), #Contamination or highly expressed poorly characterised
                            grep('^HB[BGMQDAZE][12]?',rowData(sce)$Symbol,value=TRUE), #Contamination from Hb
                            grep('^MT-',rowData(sce)$Symbol,value=TRUE), #Mitochondria
                            grep("XIST", rowData(sce)$Symbol, value = TRUE), #gender
                            grep(hkGeneREGEX,rowData(sce)$Symbol,value=TRUE) #Housekeeping genes
))
coreExcludeEnsemble <- rowData(sce)[match(coreExcludeGenes, rowData(sce)$Symbol), 1]
features_use <- features_use[!features_use %in% coreExcludeEnsemble]
exprs_mat <- exprs(sce)[features_use, ]
Y <-exprs_mat

# do the general additive model test
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
qvals = p.adjust(gam.pval,method='BH')
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:200]
heatdata <- exprs_mat[rownames(exprs_mat) %in% topgenes, 
                      order(t, na.last = NA)]
dim_red_in <- scale(t(heatdata))

#now we do a diffusion map
library(destiny)
set.seed(100)
dm <- DiffusionMap(dim_red_in, k=20)
diff_map <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
factor_plot(diff_map, sce$stages)
reducedDim(sce, "Diff_map") <- diff_map
#find a pseudotime with slingshot
library(slingshot)
sce <- slingshot(sce,clusterLabels = "stages", end.clus = "5", reducedDim = 'Diff_map')
plot(reducedDims(sce)$Diff_map, pch=16)
lines(SlingshotDataSet(sce), lwd=2)
curves <- slot(SlingshotDataSet(sce), "curves")


