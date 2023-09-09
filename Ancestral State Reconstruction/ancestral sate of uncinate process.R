# title: "**Ancestral State Reconstruction (ASR) Using Parsimony and Likelihood method**"
# project for: "Deep reptilian evolutionary roots of a major avian respiratory adaptation" (DOI: )"
#               For Baysian method, please refer to Heritage (2021) (DOI:10.1101/2021.01.10.426107)
# author: "Wani2Y"
# last modified: "16/01/2023"

# load the required libraries
library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(treeio)
library(tidytree)
library(paleotree)
library(Rcpp)
library(phangorn)
library(naniar)
library(TreeTools)

# read in the informal consensus cladogram (the phylogentic tree used to map reconstructions)
arch_phylo <- read.nexus("arch_tree.nex")

### perform ASR using preferred character coding
# read in the preferred character coding of the distribution of uncinate processes as identified by preserved specimens and inferred occurence
u_state_oc <- read.csv("unciante_scar_trait_data_oc.csv", header = TRUE)

# setup the phylogenies for analysis
# create a phylogeny without branch length (set to 1)
arch_phylo_br1 <- root(arch_phylo, "Chanaresuchus")
arch_phylo_br1 <- multi2di(arch_phylo_br1)
arch_phylo_br1 <- compute.brlen(arch_phylo_br1, 1)

# read in fossil age information collected from the literature
tip.age <- read.csv("ages.csv", header = TRUE, row.names = 1)

# time calibrate the phylogeny
arch_phylo_ts <- timePaleoPhy(arch_phylo_br1, timeData = tip.age, "equal", vartime = 2)
# export time scaled tree for other analysis if needed
# write.tree(arch_phylo_ts, file = "tree")

## Maximum Likelihood (ML)
# ML method using phylogeny wihtout branch length
ML_recon_br1_oc <- ace(u_state_oc[, 2], arch_phylo_br1, type = "discrete", method = "ML", model = "ER")

# visualizing the results
# transforming the results for visualization
anstate_ML_br1_oc <- as.data.frame(ML_recon_br1_oc$lik.anc)
anstate_ML_br1_oc$node <- 1:arch_phylo_br1$Nnode+Ntip(arch_phylo_br1)
cols <- setNames(palette()[1:length(unique(na.omit(u_state_oc[, 2])))],sort(unique(u_state_oc[, 2])))
pies_ML_br1_oc <- nodepie(anstate_ML_br1_oc, cols=1:3)
pies_ML_br1_oc <- lapply(pies_ML_br1_oc, function(g) g+scale_fill_manual(values = cols))

# plotting the results
p_ML_br1_oc <- ggtree(arch_phylo_br1, branch.length = "none") +
            geom_tiplab(size = 1.8) +
            geom_treescale(fontsize = 0.3, linesize = 0.25) +
            geom_inset(pies_ML_br1_oc, width = .015, height = .015, x = "node")
ggsave("ML_oc_br1.pdf", width = 80, height = 480, units = "cm", limitsize = FALSE)
write.csv(as.matrix(anstate_ML_br1_oc), file = "ML_oc_br1.csv", row.names = FALSE)

# ML method using time claibrated phylogeny
ML_recon_ts_oc <- ace(u_state_oc[, 2], arch_phylo_ts, type = "discrete", method = "ML", model = "ER")

# visualizing the results
# transforming the results for visualization
anstate_ML_ts_oc <- as.data.frame(ML_recon_ts_oc$lik.anc)
anstate_ML_ts_oc$node <- 1:arch_phylo_ts$Nnode+Ntip(arch_phylo_ts)
pies_ML_ts_oc <- nodepie(anstate_ML_ts_oc, cols=1:3)
pies_ML_ts_oc <- lapply(pies_ML_ts_oc, function(g) g+scale_fill_manual(values = cols))

# plotting the results
P_ML_ts_oc <- ggtree(arch_phylo_ts) +
           geom_tiplab(size = 1.8) +
           geom_treescale(fontsize = 0.3, linesize = 0.25) +
           geom_inset(pies_ML_ts_oc, width = .015, height = .015, x = "node")
ggsave("ML_oc_ts.pdf", width = 80, height = 480, units = "cm", limitsize = FALSE)
write.csv(as.matrix(anstate_ML_ts_oc), file = "ML_oc_ts.csv", row.names = FALSE)

## Maximal Parsimony (MP)
# transform distribution of uncinate process into phyDat format
u_mp_oc <- as.phyDat(matrix(u_state_oc$X, dimnames = list(u_state_oc$label)), type="USER", levels = c("0", "1", "2"), ambiguity=NA)

# MP using phylogeny witout branch length
MP_br1_oc <- ancestral.pars(arch_phylo_br1, u_mp_oc, type = "ACCTRAN", return = "prob")

# visualizing the results
# transforming the results from phangorn for visualization
anstate_mp_br1_oc <- data.frame(matrix(unlist(MP_br1_oc[1026:2049]), nrow = length(MP_br1_oc[1026:2049]), byrow = TRUE))
colnames(anstate_mp_br1_oc) <- c("0", "1", "2")

# plotting the results
anstate_mp_br1_oc <- as.data.frame(anstate_mp_br1_oc)
anstate_mp_br1_oc$node <- 1:arch_phylo_br1$Nnode+Ntip(arch_phylo_br1)
pies_MP_br1_oc <- nodepie(anstate_mp_br1_oc, cols=1:3)
pies_MP_br1_oc <- lapply(pies_MP_br1_oc, function(g) g+scale_fill_manual(values = cols))
P_MP_br1_oc <- ggtree(arch_phylo_br1, branch.length = "none") +
            geom_tiplab(size = 1.8) +
            geom_treescale(fontsize = 0.3, linesize = 0.25) +
            geom_inset(pies_MP_br1_oc, width = .015, height = .015, x = "node")
ggsave("MP_oc_br1.pdf", width = 80, height = 480, units = "cm", limitsize = FALSE)
write.csv(as.matrix(anstate_mp_br1_oc), file = "MP_oc_br1.csv", row.names = FALSE)

# MP using time calibrated phylogeny
MP_ts_oc <- ancestral.pars(arch_phylo_ts, u_mp_oc, type = "ACCTRAN", return = "prob")

# visualizing the results
# transforming the results from phangorn for visualization
anstate_mp_ts_oc <- data.frame(matrix(unlist(MP_ts_oc[1026:2049]), nrow = length(MP_ts_oc[1026:2049]), byrow = TRUE))
colnames(anstate_mp_ts_oc) <- c("0", "1", "2")
anstate_mp_ts_oc <- as.data.frame(anstate_mp_ts_oc)
anstate_mp_ts_oc$node <- 1:arch_phylo_ts$Nnode+Ntip(arch_phylo_ts)
pies_MP_ts_oc <- nodepie(anstate_mp_ts_oc, cols=1:3)
pies_MP_ts_oc <- lapply(pies_MP_ts_oc, function(g) g+scale_fill_manual(values = cols))

# plotting the results
P_MP_ts_oc <- ggtree(arch_phylo_ts) +
            geom_tiplab(size = 1.8) +
            geom_treescale(fontsize = 0.3, linesize = 0.25) +
            geom_inset(pies_MP_ts_oc, width = .015, height = .015, x = "node")
ggsave("MP_oc_ts.pdf", width = 100, height = 500, units = "cm", limitsize = FALSE)
write.csv(as.matrix(anstate_mp_ts_oc), file = "MP_oc_ts.csv", row.names = FALSE)

### perform ASR using the alternate coding
# read in the preferred character coding of the distribution of uncinate processes as identified by preserved specimens and inferred occurence
u_state_ac <- read.csv("unciante_scar_trait_data_ac.csv", header = TRUE)

## Maximum Likelihood (ML)
# ML method using phylogeny wihtout branch length
ML_recon_br1_ac <- ace(u_state_ac[, 2], arch_phylo_br1, type = "discrete", method = "ML", model = "ER")

# visualizing the results
# transforming the results from phangorn for visualization
anstate_ML_br1_ac <- as.data.frame(ML_recon_br1_ac$lik.anc)
anstate_ML_br1_ac$node <- 1:arch_phylo_br1$Nnode+Ntip(arch_phylo_br1)
pies_ML_br1_ac <- nodepie(anstate_ML_br1_ac, cols=1:3)
pies_ML_br1_ac <- lapply(pies_ML_br1_ac, function(g) g+scale_fill_manual(values = cols))

# plotting the results
p_ML_br1_ac <- ggtree(arch_phylo_br1, branch.length = "none") +
            geom_tiplab(size = 1.8) +
            geom_treescale(fontsize = 0.3, linesize = 0.25) +
            geom_inset(pies_ML_br1_ac, width = .015, height = .015, x = "node")
ggsave("ML_ac_br1.pdf", width = 80, height = 480, units = "cm", limitsize = FALSE)
write.csv(as.matrix(anstate_ML_br1_ac), file = "ML_ac_br1.csv", row.names = FALSE)

# ML method using time claibrated phylogeny
ML_recon_ts_ac <- ace(u_state_ac[, 2], arch_phylo_ts, type = "discrete", method = "ML", model = "ER")

# visualizing the results
# transforming the results from phangorn for visualization
anstate_ML_ts_ac <- as.data.frame(ML_recon_ts_ac$lik.anc)
anstate_ML_ts_ac$node <- 1:arch_phylo_ts$Nnode+Ntip(arch_phylo_ts)
pies_ML_ts_ac <- nodepie(anstate_ML_ts_ac, cols=1:3)
pies_ML_ts_ac <- lapply(pies_ML_ts_ac, function(g) g+scale_fill_manual(values = cols))

# plotting the results
P_ML_ts_ac <- ggtree(arch_phylo_ts) +
            geom_tiplab(size = 1.8) +
            geom_treescale(fontsize = 0.3, linesize = 0.25) +
            geom_inset(pies_ML_ts_ac, width = .015, height = .015, x = "node")
ggsave("ML_ac_ts.pdf", width = 80, height = 480, units = "cm", limitsize = FALSE)
write.csv(as.matrix(anstate_ML_ts_ac), file = "ML_ac_ts.csv", row.names = FALSE)

## Maximal Parsimony (MP)
# transform distribution of uncinate process into phyDat format
u_mp_ac <- as.phyDat(matrix(u_state_ac$X, dimnames = list(u_state_ac$label)), type="USER", levels = c("0", "1", "2"), ambiguity=NA)

# MP using phylogeny witout branch length
MP_br1_ac <- ancestral.pars(arch_phylo_br1, u_mp_ac, type = "ACCTRAN", return = "prob")

# visualizing the results
# transforming the results from phangorn for visualization
anstate_mp_br1_ac <- data.frame(matrix(unlist(MP_br1_ac[1026:2049]), nrow = length(MP_br1_ac[1026:2049]), byrow = TRUE))
colnames(anstate_mp_br1_ac) <- c("0", "1", "2")
anstate_mp_br1_ac <- as.data.frame(anstate_mp_br1_ac)
anstate_mp_br1_ac$node <- 1:arch_phylo_br1$Nnode+Ntip(arch_phylo_br1)
pies_MP_br1_ac <- nodepie(anstate_mp_br1_ac, cols=1:3)
pies_MP_br1_ac <- lapply(pies_MP_br1_ac, function(g) g+scale_fill_manual(values = cols))

# plotting the results
P_MP_br1_ac <- ggtree(arch_phylo_br1, branch.length = "none") +
            geom_tiplab(size = 1.8) +
            geom_treescale(fontsize = 0.3, linesize = 0.25) +
            geom_inset(pies_MP_br1_ac, width = .015, height = .015, x = "node")
ggsave("MP_ac_br1.pdf", width = 80, height = 480, units = "cm", limitsize = FALSE)
write.csv(as.matrix(anstate_mp_br1_ac), file = "MP_ac_br1.csv", row.names = FALSE)

# MP using time calibrated phylogeny
MP_ts_ac <- ancestral.pars(arch_phylo_ts, u_mp_ac, type = "ACCTRAN", return = "prob")

# visualizing the results
# transforming the results from phangorn for visualization
anstate_mp_ts_ac <- data.frame(matrix(unlist(MP_ts_ac[1026:2049]), nrow = length(MP_ts_ac[1026:2049]), byrow = TRUE))
colnames(anstate_mp_ts_ac) <- c("0", "1", "2")
anstate_mp_ts_ac <- as.data.frame(anstate_mp_ts_ac)
anstate_mp_ts_ac$node <- 1:arch_phylo_ts$Nnode+Ntip(arch_phylo_ts)
pies_MP_ts_ac <- nodepie(anstate_mp_ts_ac, cols=1:3)
pies_MP_ts_ac <- lapply(pies_MP_ts_ac, function(g) g+scale_fill_manual(values = cols))

# plotting the results
P_MP_ts_ac <- ggtree(arch_phylo_ts) +
           geom_tiplab(size = 1.8) +
           geom_treescale(fontsize = 0.3, linesize = 0.25) +
           geom_inset(pies_MP_ts_ac, width = .015, height = .015, x = "node")
ggsave("MP_ac_ts.pdf", width = 100, height = 500, units = "cm", limitsize = FALSE)
write.csv(as.matrix(anstate_mp_ts_ac), file = "MP_ac_ts.csv", row.names = FALSE)