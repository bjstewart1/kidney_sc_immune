## Bulk pseudodepth approach
pseudobulk <- pblapply(unique(sce$Short_Sample), function(x){
  message(x)
  rs <- Matrix::rowSums(counts(sce)[, sce$Short_Sample == x])
  return(rs)
})
pseudobulk <- do.call(cbind, pseudobulk)
colnames(pseudobulk) <- unique(sce$Short_Sample)
#then calculate TPM
TPM <- calculateTPM(pseudobulk)
pseudobulk <- log2(TPM + 1)

#read in the TCGA data
#now get the bulk TCGA data
TCGA <- readRDS(file.path(dirname(getwd()), "Kidney RNAseq datasets", "Data/TCGA kidney normal_processed.RDS"))
#subset this to samples which have a spatial label
TCGA <- TCGA[, colData(TCGA)$cluster %in% c("Cortex enriched", "Corticomedulla enriched", "Medulla enriched")]

require(gam)
t <-  colData(TCGA)$pt
# for time, only look at the most variable genes
Y <- assay(TCGA, "exprs")[metadata(TCGA)$hvg, ]
#select genes that have high p values testing against pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:1000]

#RF predictions
set.seed(100)
common_genes <- intersect(rownames(pseudobulk), topgenes)
cluster <- factor(colData(TCGA)$cluster)
data_train <- data.frame(t(assay(TCGA, "exprs")[common_genes, ]), cluster = cluster)

#set up test data
data_test <- data.frame(t(pseudobulk[common_genes, ]))
cs <- colSums(data_test) == 0 
cleandat <- sva::ComBat(t(as.matrix(data_test[ !cs])), droplevels(annot$Project))
data_test <- t(cleandat)

#random forests model training
rf_model <- ranger::ranger(formula = cluster~., data =data_train[, !cs] , probability = TRUE)
predictions <- predict(rf_model, data =data_test)
pt <- get_pseudotime(pc$x, 1)
pseudospace_score <- pt$Pseudotime
names(pseudospace_score) <- colnames(pseudobulk)
plot_pseudotime(pt)

#categorical predictions
common_genes <- intersect(rownames(pseudobulk), topgenes)
cluster <- factor(colData(TCGA)$cluster)
data_train <- data.frame(t(assay(TCGA, "exprs")[common_genes, ]), cluster = cluster)
data_test <- cleandat
rf_model <- ranger::ranger(formula = cluster~., data =data_train[, !cs] , probability = FALSE)
predictions <- predict(rf_model, data = t(data_test))
pseudospace_category <- predictions$predictions
names(pseudospace_category) <- colnames(pseudobulk)

#### Bulk GSVA scores
eset <- readRDS( file.path(dirname(getwd()), "Kidney RNAseq datasets", "Data/Lindgren processed data.RDS"))
pt <- readRDS(file.path(dirname(getwd()), "Kidney RNAseq datasets", "Data/Lindgren pt.RDS"))
pt <- pt$Pseudotime

library(GSVA)
gsva_in <- Y <- exprs(eset)
rownames(gsva_in) <- fData(eset)$Symbol
nephron_markers <- readRDS("Data/nephron_markers.RDS")
markers_use <- pblapply(nephron_markers, function(m){m[1:50, "Gene"]})

gsva_out <- as.matrix(gsva(gsva_in, markers_use, method = "gsva", parallel.sz=1))
gsva_plot <- pblapply(c(rownames(gsva_out)), function(g){
  ggplot(data.frame("Pseudodepth" = pt,
                    "Score" = gsva_out[g, ]),
         aes(x = Pseudodepth, y = Score)) + geom_point(size = 1, color = "grey") + geom_smooth(color = "black", fill = "grey")  + ggtitle(g)
})
plot_grid(plotlist = gsva_plot, ncol = 4)

#Plot chemokines in pseudospace
chemokines <- c("CXCL12", "CXCL14", "CXCL1", "IL8", "CXCL5", "CXCL6", "CXCL17", "CX3CL1", "CCL20", "CCL2", "CCL15")
rownames(gsva_in)[rownames(gsva_in) %in% chemokines]
chemo_plot <- pblapply(chemokines, function(g){
  ggplot(data.frame("Pseudodepth" = pt,
                    "Expression" = scale(gsva_in[g, ])),
         aes(x = Pseudodepth, y = Expression)) + geom_point(size = 1, color = "grey") + geom_smooth(color = "black", fill = "grey")  + ggtitle(g)
})
plot_grid(plotlist = chemo_plot, ncol = 4)






