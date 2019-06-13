#cellphonedB
source(file.path(dirname(getwd()), "Resources", "Functions.R"))

path <- "Adult nephron and immune 10%"
res_p <- read.table(file.path("Data/cellphone db", path, "pvalues.txt"), header= TRUE, sep="\t")
significant <- read.table(file.path("Data/cellphone db", path, "significant_means.txt"), header= TRUE, sep="\t")
means <- read.table(file.path("Data/cellphone db", path, "means.txt"), header= TRUE, sep="\t")
#limit to secreted signals
significant <- significant[significant$secreted %in% "True", ]
#don't look at integrins
significant <- significant[significant$isIntegrin %in% "False", ]
nephron_celltypes <- c("Podocyte", "PT", "Loop.of.Henle", "Distal.tubule", "Intercalated.cell", "Principal.cell", "Pelvic.epithelium")
imm_cells <- c("Mast.cell",  "MNP1", "MNP2", "MNP3", "MNP4", "Neutrophil", "Plasmacytoid.dendritic.cell")
#generate index
interact_index <- data.frame(do.call(rbind, strsplit(colnames(significant)[11:ncol(significant)], "_")))
interact_index[, 3] <- 1:nrow(interact_index)
colnames(interact_index) <- c("Partner 1", "Partner 2", "Index")
interact_index <- interact_index[interact_index$`Partner 1` %in% nephron_celltypes, ]
interact_index <- interact_index[interact_index$`Partner 2` %in% imm_cells, ]
interact_index$`Partner 1` <- factor(interact_index$`Partner 1`, levels = nephron_celltypes)
interact_index <- interact_index[order(interact_index$`Partner 1`), ]

significant <- significant[, c(1:10, 10+interact_index$Index)]

#now just take chemokines and cytokines
chemo_idx <- unique(c(
  #grep("IL[[:digit:]]", significant$interacting_pair), 
  grep("CX", significant$interacting_pair), grep("CC", significant$interacting_pair), grep("XC", significant$interacting_pair)))
significant <- significant[chemo_idx, ]
#remove rows with no significant pairs
sign_in <- significant[!rowSums(is.na(significant[, 11:ncol(significant)])) == ncol(significant[, 11:ncol(significant)]), ]
pairs <- sign_in$interacting_pair
sign_in <- sign_in[, 11:ncol(significant)]

rownames(means) <- means$interacting_pair
rownames(sign_in) <- pairs
sign_in <- means[rownames(sign_in), colnames(sign_in)]

#scale the rows
sign_in <- t(scale(t(sign_in)))

#get the associated p values
rownames(res_p) <- res_p$interacting_pair
p_vals <- res_p[rownames(sign_in), colnames(sign_in)]
p_vals <- -log10(p_vals + 0.001)

means_melt <- reshape::melt(data.frame(sign_in, "Pair" = rownames(sign_in)))
p_melt <- stack(p_vals)
df_plot <- data.frame(p = p_melt$values, means_melt)
colnames(df_plot) <- c("P_value", "Interaction", "Cell_Cell", "Mean")
df_plot$Interaction <- as.factor(df_plot$Interaction)
df_plot$Cell_Cell <- as.factor(df_plot$Cell_Cell)
ggplot(df_plot, aes(x = Cell_Cell, y = Interaction, color = Mean, size = P_value )) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_color_gradientn(colours = c("dodgerblue", "orange"), 
                                                                                                
