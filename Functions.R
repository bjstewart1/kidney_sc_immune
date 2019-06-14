#load some key libraries
library(scran)
library(scater)
library(cowplot)
library(reticulate)
use_python("/usr/bin/python")
library(umap)
library(pheatmap)
library(DropletUtils)
library(igraph)
library(SoupX)
library(pbapply)

#plot out a feature on a ggplot - blue grey red color scheme that we like
#' @param layout - a layout to use, be it PCA, tSNE, UMAP, graph layout
#' @param gene - the gene symbol to use (character)
#' @param sce - the sce object to reference
#' @param xlab - the x axis label
#' @param ylab - the y axis label
feature_plot <- function(layout, gene, sce, xlab = "Dim1", ylab="Dim2", size = 0.05, scale = FALSE){
  idx <- which(rowData(sce)$Symbol == gene)
  set.seed(100)
  #this prevents biased overplotting
  scramble <- sample(1:nrow(layout), nrow(layout))
  if(scale == TRUE){
    exprs_in <- scale(exprs(sce)[idx, scramble ])
  }else{(exprs_in <- exprs(sce)[idx, scramble ])}
  ggplot(data.frame("x" = layout[scramble, 1], "y" = layout[scramble, 2], "Expression" = exprs_in),
         aes(x = x, y=y, col = Expression)) + geom_point(pch=19, cex=size) + xlab(xlab) + ylab(ylab)+
    scale_color_gradientn(colours = c("blue", "grey", "red")) + ggtitle(gene) + theme_classic()
}

#feature violin plots
#' @param sce - sce object to use
#' @param gene - gene to use
#' @param group - groupings (eg clusters, celltypes)
feature_violin_plot <- function(sce, gene, group){
  idx <- which(rowData(sce)$Symbol == gene)
  dat <- data.frame("x" = as.factor(group), "y" = exprs(sce)[idx, ])
  ggplot(dat, aes(x= x, y=y, fill = x)) + geom_violin(scale = "width") + 
    #ggbeeswarm::geom_quasirandom(size= 0.01)
    ylab("Expression") +
    xlab("") + theme_minimal() + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(gene)
}

#plot by factor
#' @param layout - a layout to use, be it PCA, tSNE, UMAP, graph layout
#' @param factor - a factor to color by
#' @param xlab - the x axis label
#' @param ylab - the y axis label
factor_plot <- function(layout, factor, xlab = "Dim1", ylab ="Dim2", size = 0.05){
  set.seed(100)
  #this prevents biased overplotting
  scramble <- sample(1:nrow(layout), nrow(layout))
  ggplot(data.frame("x" = layout[scramble, 1], "y" = layout[scramble, 2], "Factor" = as.factor(factor)[scramble]), 
         aes(x = x, y=y, col = Factor)) + geom_point(pch=19, cex=size) + xlab(xlab) +ylab(ylab)+ theme_classic() + 
    guides(colour = guide_legend(override.aes = list(size=5)))
}

#annotated factor plot
#' @param layout  - a layout to use, be it PCA, tSNE, UMAP, graph layout
#' @param groups - a factor to group by
average.coords <- function(layout, groups){
  groups <- as.factor(groups)
  dat <- data.frame("X" = layout[, 1], "Y" = layout[, 2])
  centroids <- do.call(rbind, lapply(levels(groups), function(x){
    cm <- colMeans(dat[groups == x, ])
    return(cm)
  }))
  rownames(centroids) <- levels(groups)
  return(centroids)
}

#' @param layout  - a layout to use, be it PCA, tSNE, UMAP, graph layout
#' @param groups - a factor to group by
#' @param xlab - x axis label
#' @param ylab - y axis label
annotated_factor_plot <- function(layout, factor, xlab = "Dim1", ylab = "Dim2", size = 0.05){
  av.coords  <- average.coords(layout, groups = factor)
  set.seed(100)
  #this prevents biased overplotting
  scramble <- sample(1:nrow(layout), nrow(layout))
  ggplot(data.frame("x" = layout[scramble, 1], "y" = layout[scramble, 2], "Factor" = as.factor(factor)[scramble]), 
         aes(x = x, y=y, col = Factor)) + geom_point(pch=19, cex=size) + xlab(xlab) +ylab(ylab)+ theme_classic() + 
    theme(legend.position="none") + annotate("text", x = av.coords[, 1], y=av.coords[, 2],label= rownames(av.coords))
}


#convenience function to wrap seurat's normalisation appraoch...
#' @param sce - sce object to normalise
seurat_normalise <- function(sce){
  require(Seurat)
  object <- CreateSeuratObject(raw.data = counts(sce))
  object <- Seurat::NormalizeData(object)
  exprs(sce) <- object@data
  return(sce)}

#' Convenience wrapper for seurat's excellent find variable genes function
#' @param sce - an SCE object which we are using
#' @param x_low - lower bound for average count expression
#' @param x_high - upper bound for average count expression
#' @param y_cut - cutoff on binned z scores
#' @param num_bin - number of bins to use.
seurat_hvg <- function(sce, x_low = 0.05, x_high = 4, y_cut = 0.5, num_bin = 20, s_plot = FALSE){
  require(Seurat)
  object <- CreateSeuratObject(raw.data= counts(sce))
  object@data <- exprs(sce)
  object <- FindVariableGenes(object, x.low.cutoff = x_low, x.high.cutoff = x_high, y.cutoff = y_cut, 
                              num.bin = num_bin, do.plot = s_plot, set.var.genes = TRUE)
  HVG_info = add_gene_names(object@hvg.info)
  HVG_info$is_hvg <- rownames(HVG_info) %in% object@var.genes
  return(list("HVG" = object@var.genes, "HVG_info" = HVG_info))
}


#SLM clustering algorithm (Blondel et al.). Requires modularity optimiser in /Code in the WD
#' @param graph An igraph object to cluster.
#' @param modularity Modularity to pass to the clustering algorithm.
#' @param resolution Resolution parameter.
#' @param algorithm Algorithm to use (1 Louvain, 3 SLM)
#' @param n.start n.start for algorithm
#' @param n.iter n.iter for algorithm
#' @random.seed Random seed for algorithm
#' @print.output Print output?
#' @return A vector of cluster identities
cluster_SLM <- function(graph, modularity = 1, resolution = 1, algorithm = 3, 
                        n.start = 100, n.iter = 10, random.seed = 0, print.output = 0){
  require(igraph)
  options(scipen = 999)
  #establish where the java file is
  modularity.dir <- file.path(normalizePath(dirname(getwd())), "Resources", "Modularity optimiser")
  ModularityJarFile <- file.path(modularity.dir, "ModularityOptimizer.jar")
  #make an edgelist
  message("getting edges and edgeweights")
  e.l <- get.data.frame(graph)
  #make it zero index friendly
  e.l[, 1] <- as.numeric(e.l[, 1]) - rep(1, nrow(e.l))
  e.l[, 2] <- as.numeric(e.l[, 2]) - rep(1, nrow(e.l))
  temp.file.location <- modularity.dir
  unique_ID <- sample(10000:99999, 1)
  edge_file <- paste(temp.file.location, "/edge_", unique_ID, 
                     ".txt", sep = "")
  output_file <- paste(temp.file.location, "/output_", unique_ID, 
                       ".txt", sep = "")
  write.table(e.l, file = edge_file, sep = "\t", row.names = FALSE, 
              col.names = FALSE)
  command <- paste("java -d64 -jar -Xmx200000m", shQuote(ModularityJarFile), 
                   shQuote(edge_file), shQuote(output_file), modularity, 
                   resolution, algorithm, n.start, n.iter, random.seed, 
                   print.output, sep = " ")
  message("Clustering algorithm at work...")
  system(command, wait = TRUE, intern = FALSE)
  ident.use <- read.table(file = output_file, header = FALSE, 
                          sep = "\t")[, 1]
  file.remove(edge_file)
  file.remove(output_file)
  options(scipen = 0)
  return(ident.use + 1)
}


#' Add gene names to a data frame from an sce
#' @param x - data frame to add gene names to
#' @param sce - SCE object
add_gene_names <- function(x){
  rn <- rownames(x)
  gn <- rowData(sce)[match(rn, rownames(rowData(sce))), "Symbol"]
  x$Gene <- gn 
  return(x)
}


#' Calculates the tf-idf for a set of target cells against a background of the cells given in "universe".
#' From Matt Youngs script for kidney cells.
#'
#' @param data The data matrix to use.
#' @param target Columns that are the target.
#' @param universe Columns that we should consider (target must be a subset).
tfidf = function(data,target,universe){
  if(!all(target %in% universe))
    stop('Target must be a subset of universe')
  nObs = Matrix::rowSums(data[,target,drop=FALSE]>0)
  nTot = Matrix::rowSums(data[,universe,drop=FALSE]>0)
  tf = nObs/length(target)
  idf = log(length(universe)/nTot)
  score = tf*idf
  #Calculate p-value for significance based on using a hypergeometric distribution to simulate the results of infinite random sampling
  pvals = phyper(nObs-1,nTot,length(universe)-nTot,length(target),lower.tail=FALSE)
  qvals = p.adjust(pvals,method='BH')
  ntf = (exp(-idf)*length(universe)-tf*length(target))/(length(universe)-length(target))
  return(data.frame(geneFrequency=tf,
                    geneFrequencyOutsideCluster=ntf,
                    geneFrequencyGlobal=exp(-idf),
                    geneExpression=Matrix::rowMeans(data[,target,drop=FALSE]),
                    geneExpressionOutsideCluster = Matrix::rowMeans(data[,universe[!(universe%in%target)],drop=FALSE]),
                    geneExpressionGlobal = Matrix::rowMeans(data),
                    idf=idf,
                    tfidf=score,
                    qval=qvals)[order(score,decreasing=TRUE),])
}

#find all markers with tfid
#' @param sce - an SCE object
#' @param groups - a set of cluster identities length of ncol(sce)
tfidf_all_markers <- function(sce, groups){
  groups <- factor(groups)
  tfid.list <- list()
  pb <- txtProgressBar(min = 1, max = length(levels(groups)), style = 3)
  for(i in levels(groups)){
    tfid <- add_gene_names(tfidf(data = exprs(sce), target = colnames(sce)[groups == i], universe = colnames(sce)))
    tfid.list[[i]] <- tfid
    setTxtProgressBar(pb, match(i, levels(groups)))
  }
  return(tfid.list)    
}

#' Wilcox rank sum test between two groups
#' @param sce - an SCE object
#' @param groups - cluster identities same length as ncol(sce)
#' @param group1 - first group comparitor
#' @param group2 - second group comparitor
#' @param genes.use - genes test (rownames of sce)
wilcox_test_sc <- function(sce, groups, group1, group2, genes.use){
  groups <- factor(groups)
  sce.use <-sce
  sce.use$groups <- groups
  sce.use <- sce.use[genes.use, sce.use$groups %in% c(group1, group2)]
  #remove genes which are not expressed
  unexpressed <- Matrix::rowSums(exprs(sce.use)) == 0 
  sce.use <- sce.use[unexpressed == FALSE, ]
  #calculate lfc between the two groups
  global.av <- Matrix::rowMeans(exprs(sce.use))
  gp1.dat <- Matrix::rowMeans(exprs(sce.use)[, sce.use$groups %in% group1])
  gp2.dat <- Matrix::rowMeans(exprs(sce.use)[, sce.use$groups %in% group2])
  lfc = gp1.dat - gp2.dat
  
  #now apply a wilcox test between the two clusters
  w.out <-  unlist(lapply(rownames(sce.use), function(x){
    wilcox.test(exprs(sce.use)[x, ]~sce.use$groups)$p.value
  }))
  #multiple testing correction
  fdr = p.adjust(w.out, method='BH')
  
  #assemble and return
  gene.symbols <- rowData(sce.use)[match(rownames(sce.use), rowData(sce.use)$ID), 2]
  df <- data.frame("gene" = as.character(gene.symbols), "Global average" = global.av, "Av expr gp1" = gp1.dat, "Av expr gp2" = gp2.dat, "LFC" = lfc, "absLFC" = abs(lfc), "p.val" = w.out, "fdr.p.val" = fdr) 
  df <- df[order(df$fdr.p.val, decreasing = FALSE), ]
  return(df)
}

#find all markers with wilcox
#' @param sce - an SCE object
#' @param groups - a set of cluster identities length of ncol(sce)
#' @param genes.use - the genes to use in the test
wilcox_markers_all_groups <- function(sce, groups, genes.use){
  groups <- factor(groups)
  wilcox.list <- pblapply(levels(groups), function(i){
    gp <- ifelse(groups == i, i, "o")
    wt <- wilcox_test_sc(sce, groups = gp, group1 = i, group2 = "o", genes.use = genes.use )
    return(wt)    
  })
  return(wilcox.list)    
}

####empty drops
#empty drops wrapper
#' @param sce - an SCE object
#' @param FDR_cut - significance cutoff for calling cells
#' @param data_to_sce - add data to the sce object?
#' @return - an SCE with this data added
empty_drops_wrapper <- function(sce, FDR_cut = 0.01){
  require(DropletUtils)
  message("Running emptyDrops")
  set.seed(100)
  e.out <- DropletUtils::emptyDrops(counts(sce))
  is.cell <- e.out$FDR <= FDR_cut
  idx <- is.cell %in% TRUE
  br.out <- barcodeRanks(counts(sce))
  #subset the sce
  sce <- sce[, idx]
  #add this bunch of data to sce
  message("Adding data to SCE")
  colData(sce)$knee <- (br.out$total > br.out$knee)[idx]
  colData(sce)$inflection <-  (br.out$total > br.out$inflection)[idx]
  colData(sce)$rank <- br.out$rank[idx]
  colData(sce)$total <- br.out$total[idx]
  colData(sce)$emptyFDR <- e.out$FDR[idx]
  return(sce)
}

#function to read in an empty drops sce and also create soupX channellist object all in one
#' @param dataDirs path to raw 10X outputs - "blah blah blah...raw_gene_bc_matrices/GRCh38" 
#' @param FDR_cut - the FDR cut to use in empty drops
#' @return channels which can be used in downstream soupX to return a corrected toc 
soupX_emptydrops <- function (dataDirs, FDR_cut = 0.05) {
  require(SoupX)
  channelNames = sprintf("Channel%d", seq_along(dataDirs))
  channels = list()
  sce_list <- list()
  for (i in seq_along(dataDirs)) {
    message(sprintf("Loading data for 10X channel %s from %s", 
                    channelNames[i], dataDirs[i]))
    dataDir = dataDirs[i]
    require(DropletUtils)
    #here we read in a table of drops...
    tryCatch({
      sce <- DropletUtils::read10xCounts(samples=dataDir, col.names = TRUE)
      tod <- counts(sce)
      #here we do empty drops
      require(DropletUtils)
      message("Running emptyDrops")
      set.seed(100)
      e.out <- DropletUtils::emptyDrops(tod)
      is.cell <- e.out$FDR <= FDR_cut
      cellIdxs <- is.cell %in% TRUE
      channels[[channelNames[i]]] = SoupChannel(tod = tod, 
                                                toc = tod[,cellIdxs, drop = FALSE], 
                                                channelName = channelNames[i], 
                                                dataType = "10X")
      sce_list[[i]] <- sce[, cellIdxs]
    }, error = function(e){message("Error raised, skipping")})
  }
  #bind together the channels into a soupchannellist object
  channels = SoupChannelList(channels)
  return(list("Channels" = channels, "SCE" = sce_list))
}

#' function to get a weighted knn graph
#' @param input data - the data on which we do the KNN search
#' @param k - the k parameter
knn_graph <- function(input_data, k=k){
  library(FNN)
  knn.ir <- get.knn(input_data, k=k)
  distances <- knn.ir$nn.dist
  knn.ir <- knn.ir$nn.index
  p=1
  el = matrix(nrow = nrow(knn.ir)*ncol(knn.ir), ncol = 3)
  for(i in c(1:nrow(knn.ir))){
    for(j in 1:ncol(knn.ir)){
      node.i <- i
      node.j <- knn.ir[i, j]
      weight <- distances[i, j]
      el[p, ] <- c(node.i, node.j, weight)
      p=p+1
    }
  }
  require(igraph)
  g <- graph_from_edgelist(el[, c(1:2)], directed = FALSE)
  E(g)$weight <- el[, 3]
  g <- igraph::simplify(g)
  return(g)
}

# get pseudotime - this is essentially a simple rewrite of embeddr code.
#' @param embedding - a 2D embedding through which to fit a principal curve
#' @param clusters - not strictly necessary, but a grouping of the data
get_pseudotime <- function(embedding, clusters){
  X <- as.matrix(embedding)
  pc <- princurve::principal.curve(x = X)
  pst <- pc$lambda
  ## rescale pseudotimes to be in [0, 1]
  pst <- (pst - min(pst))/(max(pst) - min(pst))
  
  ## find orthogonal distances to projections
  d <- sqrt(rowSums(X - pc$s)^2)
  
  proj_dist <- pseudotime <- trajectory_1 <- trajectory_2 <- rep(NA, nrow(X))
  pseudotime <- pst
  trajectory_1 <- pc$s[, 1]
  trajectory_2 <- pc$s[, 2]
  proj_dist <- d
  
  out <- data.frame("Embedding_1" = embedding[, 1],
                    "Embedding_2" = embedding[, 2],
                    "Trajectory_1" = trajectory_1,
                    "Trajectory_2" = trajectory_2,
                    "Pseudotime" = pseudotime,
                    "Clusters" = clusters
  )
  return(out)
}

# convenience function to plot pseudotime
#' @param pseudotime - output of get_pseudotime()
plot_pseudotime <- function(pseudotime, dot_size = 2, lwd = 2){
  ggplot(pseudotime[order(pseudotime$Pseudotime), ]) +
    geom_point(aes(x = Embedding_1, y=Embedding_2, color = Clusters, size = dot_size)) +
    geom_path(aes(x =Trajectory_1, y=Trajectory_2), linetype = 5, size = lwd, color = "grey") + theme_classic()
}


#uses the AddModuleScore function in seurat 
#' @param sce - the sce to use
#' @param genelist - the genelist (ensemble genes)
#' @param return the scores
seurat_test_genesets <- function(sce, genelist = ensemble_genelist){
  sobj <- Seurat::CreateSeuratObject(counts(sce))
  sobj <- Seurat::NormalizeData(sobj)
  sobj <- Seurat::AddModuleScore(sobj, features = genelist)
  scores <- sobj@meta.data[, grep("Cluster", colnames(sobj@meta.data))]
  colnames(scores) <- names(genelist)
  return(scores)
}


#set a color palette for heatmaps of single cell data
heat_colors_sc <- colorRampPalette(c("blue", "dodgerblue","lightgrey", "red", "darkred"))(20)


#define a color scheme for p values 
p_scale <- colorRampPalette(c("purple4", "limegreen"))

#' Plots a series of barplots and connects them  - from JEFworks
#' Modified from https://stackoverflow.com/questions/22560850/barplot-with-connected-series
#' 
#' @param dat NxM matrix with N rows as features and M columns as samples
#' @param color Vector of N colors
#' @param space Space between barplots
#' @param alpha Alpha for area connecting barplots
#' 
#' @examples
#' dat <- matrix(rnorm(100),10,10)
#' dat <- abs(matrix(rnorm(100),10,10)) 
#' connectedBarplot(dat, color=rainbow(nrow(dat)))
#'
connectedBarplot <- function(dat, color=rainbow(nrow(dat)), space=1, alpha=0.5, ...) {  
  b <- barplot(dat, col=color, space = space, ...)                     
  
  for (i in seq_len(ncol(dat) - 1)) {     
    lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## bottom line       
    
    for (j in seq_len(nrow(dat))) {     
      if (j == 1) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))                       
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(0, dat[j,i], dat[j,i+1], 0),               
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j == 2) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),                       
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j > 2) {                    
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], colSums(dat[1:(j-1),])[i+1]),              
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
    }          
  }              
}      


#do a mean mean plot
#' @cluster_1 - the x axis cluster
#' @cluster_2 - the y axis cluster
#' @cluster_universe - all the cluster IDS in the sce
#' @sce the sce
#' @size - point size
#' @nselect - the number of genes with high fold changes to plot

mean_mean_plot <- function(cluster_1, cluster_2, cluster_universe, sce, size = 0.5, nselect = 10){
  require(ggrepel)
  xax <- calcAverage(sce[, cluster_universe == cluster_1], use_size_factors=FALSE, exprs_values = "logcounts")
  yax <- calcAverage(sce[, cluster_universe == cluster_2], use_size_factors=FALSE, exprs_values = "logcounts")
  lfc <- xax - yax
  dat <-  data.frame(xax, yax, lfc, "Symbol" = rowData(sce)$Symbol)
  lfc <- lfc[order(lfc)]
  top_genes <- c(names(head(lfc, nselect)), names(tail(lfc, nselect)))
  dat$label <- NA
  dat[top_genes, "label"] <- as.character(dat[top_genes, "Symbol"])
  dat$color <-2
  dat[top_genes, "color"] <- 1
  dat[dat$lfc == 0, "label"] <- NA
  dat[dat$lfc == 0, "color"] <- 2
  ggplot(dat, aes(x= xax, y=yax, label = label, color = color)) +xlab(paste(cluster_1))+ylab(paste(cluster_2)) + geom_point(size =size) + 
    ggrepel::geom_text_repel() + theme(legend.position="none") + coord_fixed()
}

