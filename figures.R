####### figures for Estrada & Marshall (2024) #######

library(tidyverse)
library(phytools)
library(ape)
library(ggpie)
library(ggtree)
library(ggimage)

#read in primate data 
primates <- read.csv(file = "data/primates_final.csv", header = TRUE)

######################## PHYLOGENETIC TREES & ASR #################################################################

#file to align 10kTrees and IUCN taxonomies
taxa <- read.csv(file = "data/10k_primate_tree_taxa.csv", header = TRUE)
colnames(taxa) <- c("X10k_tax", "scientific_name", "remove")

#add terrestriality data
taxa <- left_join(taxa, 
                  primates[,c("scientific_name", "Terrestrial.binary", 
                              "Terrestrial.ordinal","ground.use", "IUCN.family", )], 
                  by = "scientific_name")

#10kTrees phylogeny for taxa with ordinal terrestriality data only
ordinal_tree <- read.nexus(file = "data/10kTrees_Primates_ordinal.nex")

#remove underscores from Latin binomials
ordinal_tree$tip.label <- sub("_"," ", ordinal_tree$tip.label)
ordinal_tree$tip.label <- sub("_"," ", ordinal_tree$tip.label)

#make phylo tree dichotomous
ordinal_tree_di <- multi2di(ordinal_tree)

#ordinal terrestriality data frame
terr_ord_data <- data.frame(taxa[taxa$X10k_tax %in% ordinal_tree$tip.label, 5])
rownames(terr_ord_data) <- taxa[taxa$X10k_tax %in% ordinal_tree$tip.label, 1]
terr_ord_data <- setNames(terr_ord_data$taxa.taxa.X10k_tax..in..ordinal_tree.tip.label..5., 
                          rownames(terr_ord_data))

#setting plot colors
cols2 <- c("white", "grey", "black")

#functions to add outlines to phylo node pie charts (from https://github.com/YuLab-SMU/ggtree)
ggpie <- function(data, y, fill, color, alpha=1, outline.color="transparent", outline.size=0) {
  p <- ggplot(data, aes_(x=1, y=y, fill=fill)) +
    geom_bar(stat='identity', alpha=alpha, color=outline.color, size=outline.size, show.legend = F) +
    coord_polar(theta='y') + theme_inset()
  
  if (missingArg(color) || is.null(color) || any(is.na(color))) {
    ## do nothing
  } else {
    p <- p+scale_fill_manual(values=color)
  }
  return(p)
}

nodepie2 <- function(data, cols, color, alpha=1, outline.color="transparent", outline.size=0) {
  if (! "node" %in% colnames(data)) {
    stop("data should have a column 'node'...")
  }
  type <- value <- NULL
  if (missingArg(color)) {
    color <- NA
  }
  ldf <- gather(data, type, value, !! cols) %>% split(., .$node)
  lapply(ldf, function(df) ggpie(df, y=~value, fill=~type, color, alpha, outline.color, outline.size))
}


### phylogenetic tree with ancestral state reconstruction 

#run the ASR
ace.result2 <- ace(terr_ord_data, ordinal_tree_di, model = "ER", type = "discrete")

#plot the tree and add the pie charts on the nodes, plus an axis showing MYA
plotTree(ordinal_tree_di, fsize = 0.25, ftype = "i", lwd = 1.5)
nodelabels(node = 1:ordinal_tree_di$Nnode + Ntip(ordinal_tree_di),
           pie = ace.result2$lik.anc, piecol = cols2, cex = 0.25)
axis(1, line = -1.5, lwd = 0.25, at = (0:73),
     labels = rev(0:73), cex.axis = 0.5, padj = -2.5)
add.simmap.legend(prompt = FALSE, leg = c("arboreal", "semiterrestrial", "terrestrial"), 
                  colors = cols2, fsize = 0.8, y = 50)


### phylogenetic tree showing terrestriality data via heatmap

#create dataframe
terr_ord_data_df <- data.frame(taxa[taxa$X10k_tax %in% ordinal_tree$tip.label, 5])
rownames(terr_ord_data_df) <- taxa[taxa$X10k_tax %in% ordinal_tree$tip.label, 1]

#creating pies 
ancstats <- as.data.frame(ace.result2$lik.anc)
ancstats$node <- 1:ordinal_tree_di$Nnode + Ntip(ordinal_tree_di)
pies2 <- nodepie2(ancstats, cols = 1:3, outline.color = "black")
pies2 <- lapply(pies2, function(g) g + scale_fill_manual(values = c("white", "grey", "black")))

#phylo tree
p2 <- ggtree(ordinal_tree_di) +
  xlim(0, 100)

p2 <- gheatmap(p2, terr_ord_data_df, width = 0.03, low = "white", high = "black", 
               colnames = FALSE, font.size = 1, color = "black") +
  scale_fill_manual(values = c("white", "grey", "black"), name = "", labels = c("arboreal", "semi-terrestrial", "terrestrial")) +
  geom_cladelabel(418, "Galagoidea", offset = 7, barsize = 1, angle = 0, offset.text = 5, hjust = 0.5, fontsize = 3) + 
  geom_cladelabel(373, "Lemuroidea", offset = 7, barsize = 1, angle = 0, offset.text = 5, hjust = 0.5, fontsize = 3) +
  geom_cladelabel(370, "Tarsioidea", offset = 7, barsize = 1, angle = 0, offset.text = 2, hjust = 0.2, fontsize = 3) +
  geom_cladelabel(329, "Ceboidea", offset = 7, barsize = 1, angle = 0, offset.text = 4, hjust = 0.5, fontsize = 3) +
  geom_cladelabel(312, "Hominoidea", offset = 7, barsize = 1, angle = 0, offset.text = 5, hjust = 0.5, fontsize = 3) +
  geom_cladelabel(223, "Cercopithecoidea", offset = 7, barsize = 1, angle = 0, offset.text = 7.25, hjust = 0.5, fontsize = 3) + 
  theme(legend.position = c(0.3,0.9))

#add images from phylopic
phylopic_info <- data.frame(node = c(418, 373, 370, 338, 312, 223),
                            phylopic = c("7fb9bea8-e758-4986-afb2-95a2c3bf983d",
                                         "bac25f49-97a4-4aec-beb6-f542158ebd23",
                                         "f598fb39-facf-43ea-a576-1861304b2fe4",
                                         "aceb287d-84cf-46f1-868c-4797c4ac54a8",
                                         "0174801d-15a6-4668-bfe0-4c421fbe51e8",
                                         "72f2f854-f3cd-4666-887c-35d5c256ab0f"))
p2 %<+%
  phylopic_info + 
  geom_nodelab(aes(image = phylopic), geom = "phylopic", alpha = .5, color = "brown",
               nudge_x = c(-5,-5,-5,-6,-5,-5),
               nudge_y = c(0,0,-1,-3,0,-3))


### ASR with continuous data 

#create df
terr_con_data <- data.frame(taxa[taxa$X10k_tax %in% binary_tree$tip.label, 8])
rownames(terr_con_data) <- taxa[taxa$X10k_tax %in% binary_tree$tip.label, 1]
terr_con_data <- na.omit(terr_con_data)

#trim tree
con_tree_bi <- drop.tip(ordinal_tree_di, ordinal_tree_di$tip.label[-na.omit(match(rownames(terr_con_data), ordinal_tree_di$tip.label))])

#trim df
terr_con_data2 <- terr_con_data[which(rownames(terr_con_data) %in% con_tree_bi$tip.label),, drop = FALSE]

#change format
terr_con_data2 <- setNames(terr_con_data2$taxa.taxa.X10k_tax..in..binary_tree.tip.label..8., rownames(terr_con_data2))

#run ASR
fast_cont <- fastAnc(con_tree_bi, terr_con_data2, vars = TRUE, CI = TRUE)

#plot tree
obj <- contMap(con_tree_bi, terr_con_data2, plot = FALSE)
plot(obj, legend = 0.410939*max(nodeHeights(con_tree_bi)), leg.txt = "percent time on ground", 
     legend.pos = c(0.1,0.5), fsize = c(0.45, 0.8), lwd = 2)
