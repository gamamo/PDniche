# Run phylo niche -------------------------------------------------------------------

library(tidyverse)
library(picante)
library(treeio)
library(tidytree)
library(ggrepel)

source("phylo_niche.R")

# PART 1: Metanetworks ----------------------------------------------------

# Load data
data <- read_csv("DataForAnalysis_v26ago26_encodingFix_utf_NoDupli.csv")
data <- data |> select(-interaction) |> mutate(interaction = int) |> select (-int)

length(unique(data$plant_id))
length(unique(data$Scientific))
length(unique(data$database))

## Run the phyloniche for individual networks --------------------------------

# Transform to data wide
wider <- data |> 
  dplyr::select(Scientific, plant_id,interaction) |> 
  dplyr::rename(birds = Scientific, plants = plant_id) |> 
 # dplyr::mutate(values = 1) |> 
  dplyr::filter(!is.na(plants)) |> 
  dplyr::mutate(birds = str_replace_all(birds, " ", "_")) |> 
  dplyr::mutate(plants = str_replace_all(plants, " ", "_")) |> 
  dplyr::distinct(plants, birds, .keep_all = TRUE) |> 
  tidyr::pivot_wider(names_from = birds, values_from  = interaction,values_fill = 0) 


# Transform to a community matrix
wider <- wider |> 
  column_to_rownames("plants") |>   # This makes 'plants' the row names
  as.matrix()

# Load the phylogenetic tree
t1 <- read.tree("WRLD_phylo.tre")

intree <- intersect(rownames(wider), t1$tip.label)

length(unique(intree))

pt1 <- keep.tip(t1, intree)

# This resolve polytomies
pt1 <- multi2di(pt1) 

# now make sure that only the species in the tree are in the community data
wider <- wider[rownames(wider) %in% pt1$tip.label, ]

# make sure that all colummns and rows have sums > 0
wider <- wider[which(rowSums(wider) > 0) , which(colSums(wider) > 0)]

length(which(colSums(wider) == 0))

# the community data should be as list
nets <- list()
nets[[1]] <- wider

# calculate the PDN
pdn <- phylo_niche(networks = nets, phylotree = pt1)
b <- pdn[[2]]

# what are the species with larger proportion of PDniche?
b <- b |> 
  mutate(prop = (pdint.b*100)/PD) |> 
  arrange(desc(prop)) |> 
  mutate(s = sum(prop))

b |> write_csv("PDniche_meta_v260826.csv")


# PART 2: Individual networks ----------------------------------------------------

# Load data
data <- read_csv("DataForAnalysis_v26ago26_encodingFix_utf_NoDupli.csv")
data <- data |> select(-interaction) |> mutate(interaction = int) |> select (-int)

indCode  <- data |> 
  mutate(netcode= as.numeric(factor(database))) |> 
  select(database, netcode) |> 
  distinct()

ind <- data |> 
  mutate(netcode = as.numeric(factor(database))) |>
  group_split(database)

# Transform to data wide
ind2 <- list()

for (i in 1:length(ind)){
  print(i)

ind[[i]] |> 
  dplyr::select(Scientific, plant_id,interaction) |> 
  dplyr::rename(birds = Scientific, plants = plant_id) |> 
  dplyr::filter(!is.na(plants)) |> 
  #dplyr::mutate(values = 1) |> 
  dplyr::mutate(birds = str_replace_all(birds, " ", "_")) |> 
  dplyr::mutate(plants = str_replace_all(plants, " ", "_")) |> 
  dplyr::distinct(plants, birds, .keep_all = TRUE) |> 
  tidyr::pivot_wider(names_from = birds, values_from  = interaction, values_fill = 0 ) -> temp
  
# Transform to a community matrix
temp2 <- temp |> 
    column_to_rownames("plants") %>%  # This makes 'plants' the row names
    as.matrix() 

ind2[[i]] <- temp2

}

## Run the phyloniche for individual networks ------------------------------------------

# Load the phylogenetic tree
t1 <- read.tree("WRLD_phylo.tre")

PDlist <- list()

for (j in 1:length(ind2)){
  #for (j in 120:125){
  print(j)
  
  wider <- ind2[[j]]
  
  # j is the netcode
  intree <- intersect(rownames(wider), t1$tip.label)
  
  if(length(intree)==0){ # if tree returns no species
    PDlist[[j]] <- 0 
    next
  } 
    pt1 <- keep.tip(t1, intree)
    
    # This resolve polytomies
    pt1 <- multi2di(pt1) 
    
    # now make sure that only the species in the tree are in the community data
    wider <- wider[rownames(wider) %in% pt1$tip.label, ]
    
    # make sure that all colummns and rows have sums > 0
    if(dim(as.matrix(wider))[2] < 2 | dim(as.matrix(wider))[1] <= 2 ){ # x' must be an array of at least two
      PDlist[[j]] <- 0
      next
    } 
      wider <- wider[which(rowSums(wider) > 0) , which(colSums(wider) > 0)]
      
      # the community data should be as list
      nets <- list()
      nets[[1]] <- wider
      
      # calculate the PDN
      pdn <- phylo_niche(networks = nets, phylotree = pt1)
      PDlist[[j]] <- pdn[[2]]
      PDlist[[j]]$netcode <- j
}

# convert list to dataframe

PDind <-  PDlist |> 
  keep(is.data.frame) |>  # delete elements of hte list that are 0s
  map_df(~.x) |> 
  rename(birds=species)
      
PDind$birds <- gsub("_"," ", PDind$birds)
PDind <- PDind |> left_join(indCode, by="netcode")

PDind |> write_csv("PDniche_ind_v260826.csv")
