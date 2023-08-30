

phyll.dat <- read.csv("/home/m/Dropbox/Projects/Phyllostomid call evolution/data simulation Nazareth/phyllostomidae1.csv", stringsAsFactors = F)


head(phyll.dat)


names(phyll.dat)[2:3] <- c("Especie", "Gremio")

library(ape)

murcis <- read.nexus("~/MEGA/Phylogenies & data/murcis1.nex")
phy_rabosky <- read.tree("/home/m/Dropbox/Projects/Phyllostomid call evolution/Datos Rabosky/treePL_ML/chiroptera.rooted.dated.tre")

library("readxl")

mm_list <- read_xlsx("/home/m/Dropbox/Projects/Phyllostomid call evolution/global-mammal-checklist-at-may-2018.xlsx", range = "C4:H5740")

str(mm_list)


phyllos_list <- mm_list[mm_list$Family == "PHYLLOSTOMIDAE", ]
phyllos_list <- phyllos_list[!is.na(phyllos_list$Species), ]
phyllos_list$Species.name <- paste(phyllos_list$Genus, phyllos_list$Species, sep = "_")


plot.phylo(phy_rabosky)

plot.phylo(murcis)

rab_gen <- as.vector(sapply(phy_rabosky$tip.label, function(x) strsplit(x, "_")[[1]][1], simplify = T))

prun_rab_tree <- drop.tip(phy_rabosky, tip = which(!rab_gen %in% phyllos_list$Genus))

any(duplicated(prun_rab_tree$tip.label))

plot.phylo(prun_rab_tree, cex = 0.6,no.margin = T, type = "fan")

Ntip(prun_rab_tree)


