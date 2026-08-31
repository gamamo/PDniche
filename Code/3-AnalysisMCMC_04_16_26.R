
# Load packages -----------------------------------------------------------

library(MCMCglmm)
library(ggplot2)
library(ape)
library(tidyverse)
library(MuMIn)
library(ggridges)
library(rio)
library(car)
library(tictoc)
library(here)


# 1) Metanetwork models ---------------------------------------------------
if(F){
phy <- read.tree("AllBirdsHackett1.tre")
set.seed(123)
phy_sample <- sample(phy, size = 100)
write.tree(phy_sample, file = "sampled_trees_100.tre")
}

phy <- read.tree("sampled_trees_100.tre")

# make a looping to run MCMC for trees
tic()
list_posteriori = list()
list_intervals  = list()
list_summaries  = list()
list_sol        = list()
list_DIC        = list()
list_summaryGeral = list()
list_random = list()

for (i in 1:length(phy)){
#for (i in 1:10){
  
  print(i)
  # load data metanetwork ####
  data <- read_csv("Metanetwork_dataToAnalysis_v260826.csv")
  data <- data |> select(-interaction) |> mutate(interaction = int) |> select (-int)
  data <- data |> 
    mutate(Mass_log = log10(Mass),
           Range.Size_log = log10(Range.Size),
           pdint.b_log = log10(pdint.b),
           Beak.Width_log = log10(Beak.Width),
           Wing.Length_log = log10(Wing.Length)) |> 
    dplyr::filter(interaction==1) |> 
    #dplyr::filter(!str_detect(ref, regex("from beak to fruit", ignore_case = TRUE))) |> 
    dplyr::select(birds, pdint.b, pdint.b_log, NicheBreadth, NichePosition, Beak.Width,
                  Beak.Length_Culmen, Beak.Depth, Tarsus.Length, Wing.Length,
                  Kipps.Distance, hwi, Tail.Length, Mass, Beak.Width_log,
                  Range.Size, Mhuman_impact, Mass_log, Range.Size_log, family, Realm,
                  Mhuman_impact,Mdhuman_impact,ED,DR,Wing.Length_log) |> 
    dplyr::rename(Family = family) |> 
    distinct() |> 
    drop_na()
  
  dim(data)

  # modelling
  
  # transform some variables
  
  data$NichePosition = abs(data$NichePosition)
  
  # scale variables that will be used in the modelling and exclude NAs
  data[,c("pdint.b","pdint.b_log","NicheBreadth", "NichePosition","Beak.Width",
           "Beak.Length_Culmen", "Beak.Depth", "Tarsus.Length", "Wing.Length",
           "Kipps.Distance","hwi" ,"Tail.Length", "Mass",
           "Range.Size","Mhuman_impact","Mass_log","Range.Size_log",
           "Beak.Width_log","Mdhuman_impact","ED","DR","Wing.Length_log")] <- 
    scale(data[,c("pdint.b","pdint.b_log","NicheBreadth", "NichePosition","Beak.Width",
                   "Beak.Length_Culmen", "Beak.Depth", "Tarsus.Length", "Wing.Length",
                   "Kipps.Distance","hwi" ,"Tail.Length", "Mass",
                   "Range.Size","Mhuman_impact","Mass_log","Range.Size_log",
                   "Beak.Width_log","Mdhuman_impact","ED","DR","Wing.Length_log")])
  
  dataM <- as.data.frame(data)
  
  dataM$Family = as.factor(dataM$Family)
  dataM$Realm = as.factor(dataM$Realm)
 
  # Run MCMCglmm models metanetwork ######
  
  #use only one tree
  tree <-  phy[[i]]
  #tree <- phy
  tree$tip.label  <-  gsub("_"," ",tree$tip.label)
  
  # prepare data #####################################################
  # remove species not in the tree ####
  notintree = setdiff(dataM$birds,tree$tip.label)
  notindata = setdiff(tree$tip.label, dataM$birds)
  toexclude = c(notindata,notintree)
  
  #prune tree
  treeprune = drop.tip(tree,toexclude)
  
  #exclude species notintree from the dataset
  dataM = dataM[-which(dataM$birds %in% notintree),]
  
  # prepare the covariance matrix
  Ainv <- inverseA(treeprune)$Ainv
  
  # model without biogeographic realm
  prior<-list(R=list(V=1, nu=0.002),
              G=list(G1=list(V=1, nu=0.002),
                     G2=list(V=1, nu=0.002)))
  
  MCmodmeta <-  MCMCglmm(pdint.b_log~NicheBreadth+NichePosition+Beak.Width+Wing.Length+
                           Mass_log+ hwi+Range.Size_log+Mdhuman_impact+ED,
                       random= ~birds + Realm,
                       prior = prior,
                       ginverse=list(birds=Ainv),
                       nitt=100000, burnin = 25000 ,thin = 50,
                       data=na.omit(dataM))
  
  
  list_posteriori[[i]] = posterior.mode(MCmodmeta$Sol)
  list_intervals[[i]] = HPDinterval(MCmodmeta$Sol)
  list_summaries[[i]]  = summary(MCmodmeta$Sol)[[2]]
  list_sol[[i]]        = MCmodmeta$Sol
  list_DIC[[i]]        = DIC(MCmodmeta)
  list_summaryGeral[[i]] = summary(MCmodmeta)[[5]]
  list_random[[i]] =    summary(MCmodmeta)[[6]]
  
} # fecha looping
toc()


#save products ----------------------

list_posteriori_df <- map_df(list_posteriori,~.x) |> 
  mutate(model = "Metanetwork") |> 
  write_csv(here("productsR1","resu_mcmc_meta","list_posteriori.csv"))

list_intervals_df <-map_df(list_intervals, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "Metanetwork") |> 
  write_csv(here("productsR1","resu_mcmc_meta","list_intervals.csv"))

list_summaries_df <- map_df(list_summaries, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "Metanetwork") |> 
  write_csv(here("productsR1","resu_mcmc_meta","list_summaries.csv"))

list_sol_df <- map_df(list_sol, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "Metanetwork") |> 
  write_csv(here("productsR1","resu_mcmc_meta","list_sol.csv"))

list_dic_df <- map_df(list_DIC, ~ {
  .x |> 
    as.data.frame() 
}, .id = "iteration") |> 
  mutate(model = "Metanetwork") |> 
  rename(DIC = ".x") |> 
  write_csv(here("productsR1","resu_mcmc_meta","list_dic.csv"))

list_summary_df <- map_df(list_summaryGeral, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "Metanetwork") |> 
  write_csv(here("productsR1","resu_mcmc_meta","list_summaryGeral.csv"))

list_random_df <- map_df(list_random, ~ {
  .x |> 
    as.data.frame()|> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "Metanetwork") |> 
  write_csv(here("productsR1","resu_mcmc_meta","list_random.csv"))


if(F){
setwd("C:/Users/gabri/Dropbox/postdocINECOL/dadosCompilados/workingfolder_GM/scripts_MCMC_saturn_7abr22")
export(list_posteriori,"Posterior.xlsx")
export(list_intervals,"Intervals.xlsx")
export(list_summaries,"Summaries.xlsx")
export(list_sol,"SOL.xlsx")
export(list_DIC,"DIC.xlsx")
export(list_summaryGeral,"SummaryGeneral2.xlsx")

setwd("C:/Users/gabri/Dropbox/postdocINECOL/dadosCompilados/workingfolder_GM/scripts_MCMC_saturn_7abr22")
save.image("workspace_Modglobal.RData")


# make the graphics###
# density graphic


den <- map_df(list_posteriori, ~.x) |> 
  pivot_longer(everything(), names_to = "variable", values_to = "value") |> 
  mutate(variable = as.factor(variable))


ggplot(subset(den,variable!="(Intercept)"), aes(x=value,y=variable))+
  geom_density_ridges(scale = 2,fill = "lightblue", alpha = 0.5)+
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(-1,1))+
  geom_vline(xintercept = 0,linetype="dashed")

# distribution graphic

a= list_sol
blist=list()
for(i in 1:length(a)){
  xx= a[[i]]
  xx=as.data.frame(xx)
  xx$tree = rep(i, nrow(xx))
  xx <- xx |> pivot_longer(cols = -tree, names_to = "variable", values_to = "value")
  xx$variable = as.factor(xx$variable)
  xx$tree = as.factor(xx$tree)
  blist[[i]] = xx
}

c= map_df(blist,~.x)


g = ggplot(subset(c,variable!="(Intercept)"), aes(x=value,color=tree))+
  geom_density()+
  geom_vline(xintercept = 0,linetype="dashed")+
  theme(legend.position="none")+
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_color_grey()+
  facet_wrap(~variable,scales="free")+
  ylab("Posterior density")+
  coord_flip()
g

save.image("testMCMC.RData")
}

# 2) Individual level analysis -----------------------------------------------

phy <- read.tree("sampled_trees_100.tre")

# make a looping to run MCMC for trees
tic()
list_posterioriI = list()
list_intervalsI  = list()
list_summariesI  = list()
list_solI        = list()
list_DICI        = list()
list_summaryGeralI = list()
list_randomI = list()

for (i in 1:length(phy)){
#for (i in 1:10){
  
  print(i)
  # load data ind ####
  data <- read_csv( "Indnetwork_dataToAnalysis_v260826.csv")
  data <- data |> 
    mutate(Mass_log = log10(Mass),
           Range.Size_log = log10(Range.Size),
           pdint.b_log = log10(pdint.b),
           Beak.Width_log = log10(Beak.Width),
           Wing.Length_log = log10(Wing.Length),
           ctrbpdint.b_log = log10(ctrbpdint.b)) |> 
    #dplyr::filter(interaction==1) |> 
    #dplyr::filter(!str_detect(ref, regex("from beak to fruit", ignore_case = TRUE))) |> 
    dplyr::select(birds, pdint.b, pdint.b_log, NicheBreadth, NichePosition, Beak.Width,
                  Beak.Length_Culmen, Beak.Depth, Tarsus.Length, Wing.Length,
                  Kipps.Distance, hwi, Tail.Length, Mass, Beak.Width_log,
                  Range.Size, Mass_log, Range.Size_log, family, Realm,
                  human_impact,ED,DR,Wing.Length_log,ctrbpdint.b_log) |> 
    dplyr::rename(Family = family) |> 
    distinct() |> 
    drop_na()
  
  dim(data)

  # modelling
  
  # transform some variables
  
  data$NichePosition = abs(data$NichePosition)
  
  # scale variables that will be used in the modelling and exclude NAs
  data[,c("pdint.b","pdint.b_log","NicheBreadth", "NichePosition","Beak.Width",
          "Beak.Length_Culmen", "Beak.Depth", "Tarsus.Length", "Wing.Length",
          "Kipps.Distance","hwi" ,"Tail.Length", "Mass",
          "Range.Size","Mass_log","Range.Size_log",
          "Beak.Width_log","human_impact","ED","DR","Wing.Length_log",
          "ctrbpdint.b_log")] <- 
    scale(data[,c("pdint.b","pdint.b_log","NicheBreadth", "NichePosition","Beak.Width",
                  "Beak.Length_Culmen", "Beak.Depth", "Tarsus.Length", "Wing.Length",
                  "Kipps.Distance","hwi" ,"Tail.Length", "Mass",
                  "Range.Size","Mass_log","Range.Size_log",
                  "Beak.Width_log","human_impact","ED","DR","Wing.Length_log",
                  "ctrbpdint.b_log")])
  
  dataM <- as.data.frame(data)
  
  dataM$Family = as.factor(dataM$Family)
  dataM$Realm = as.factor(dataM$Realm)
  
  # Run MCMCglmm models metanetwork ######
  
  #use only one tree
  tree <-  phy[[i]]
  #tree <- phy
  tree$tip.label  <-  gsub("_"," ",tree$tip.label)
  
  # prepare data #####################################################
  # remove species not in the tree ####
  notintree = setdiff(dataM$birds,tree$tip.label)
  notindata = setdiff(tree$tip.label, dataM$birds)
  toexclude = c(notindata,notintree)
  
  #prune tree
  treeprune = drop.tip(tree,toexclude)
  
  #exclude species notintree from the dataset
  dataM = dataM[-which(dataM$birds %in% notintree),]
  
  # prepare the covariance matrix
  Ainv <- inverseA(treeprune)$Ainv
  
  # model without biogeographic realm
  prior<-list(R=list(V=1, nu=0.002),
              G=list(G1=list(V=1, nu=0.002),
                     G2=list(V=1, nu=0.002)))
  
  if(F){
  MCmodmeta <-  MCMCglmm(pdint.b_log~NicheBreadth+NichePosition+Beak.Width+Wing.Length+
                           Mass_log+ hwi+Range.Size_log+human_impact+ED+DR,
                         random= ~birds + Realm,
                         prior = prior,
                         ginverse=list(birds=Ainv),
                         nitt=100000, burnin = 25000 ,thin = 50,
                         data=na.omit(dataM))
  }
  MCmodmeta <-  MCMCglmm(pdint.b_log~NicheBreadth+NichePosition+Beak.Width+Wing.Length+
                           Mass_log+ hwi+Range.Size_log+human_impact+ED,
                         random= ~birds + Realm,
                         prior = prior,
                         ginverse=list(birds=Ainv),
                         nitt=100000, burnin = 25000 ,thin = 50,
                         data=na.omit(dataM))

  
  list_posterioriI[[i]] = posterior.mode(MCmodmeta$Sol)
  list_intervalsI[[i]] = HPDinterval(MCmodmeta$Sol)
  list_summariesI[[i]]  = summary(MCmodmeta$Sol)[[2]]
  list_solI[[i]]        = MCmodmeta$Sol
  list_DICI[[i]]        = DIC(MCmodmeta)
  list_summaryGeralI[[i]] = summary(MCmodmeta)[[5]]
  list_randomI[[i]] = summary(MCmodmeta)[[6]]
  
} # fecha looping
toc()


#save products ----------------------

list_posterioriI_df <- map_df(list_posterioriI,~.x) |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind","list_posteriori.csv"))

list_intervalsI_df <-map_df(list_intervalsI, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind","list_intervals.csv"))

list_summariesI_df <- map_df(list_summariesI, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind","list_summaries.csv"))

list_solI_df <- map_df(list_solI, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind","list_sol.csv"))

list_dicI_df <- map_df(list_DICI, ~ {
  .x |> 
    as.data.frame() 
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  rename(DIC = ".x") |> 
  write_csv(here("productsR1","resu_mcmc_ind","list_dic.csv"))

list_summaryI_df <- map_df(list_summaryGeralI, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind","list_summaryGeral.csv"))

list_randomI_df <- map_df(list_randomI, ~ {
  .x |> 
    as.data.frame()|> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind","list_random.csv"))



if(F){
  setwd("C:/Users/gabri/Dropbox/postdocINECOL/dadosCompilados/workingfolder_GM/scripts_MCMC_saturn_7abr22")
  export(list_posteriori,"Posterior.xlsx")
  export(list_intervals,"Intervals.xlsx")
  export(list_summaries,"Summaries.xlsx")
  export(list_sol,"SOL.xlsx")
  export(list_DIC,"DIC.xlsx")
  export(list_summaryGeral,"SummaryGeneral2.xlsx")
  
  setwd("C:/Users/gabri/Dropbox/postdocINECOL/dadosCompilados/workingfolder_GM/scripts_MCMC_saturn_7abr22")
  save.image("workspace_Modglobal.RData")


# make the graphics###
# density graphic


den <- map_df(list_posteriori, ~.x) |> 
  pivot_longer(everything(), names_to = "variable", values_to = "value") |> 
  mutate(variable = as.factor(variable))


ggplot(subset(den,variable!="(Intercept)"), aes(x=value,y=variable))+
  geom_density_ridges(scale = 2,fill = "lightblue", alpha = 0.5)+
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(-1,1))+
  geom_vline(xintercept = 0,linetype="dashed")
}

# 3) Contributions --------------------------------------------------------------

phy <- read.tree("sampled_trees_100.tre")

# make a looping to run MCMC for trees
tic()
list_posterioriC = list()
list_intervalsC  = list()
list_summariesC  = list()
list_solC        = list()
list_DICC        = list()
list_summaryGeralC = list()
list_randomC = list()

for (i in 1:length(phy)){
  #for (i in 1:5){
  
  print(i)
  # load data ind ####
  data <- read_csv( "Indnetwork_dataToAnalysis_v260826.csv")
  data <- data |> 
    mutate(Mass_log = log10(Mass),
           Range.Size_log = log10(Range.Size),
           pdint.b_log = log10(pdint.b),
           Beak.Width_log = log10(Beak.Width),
           Wing.Length_log = log10(Wing.Length),
           ctrbpdint.b_log = log10(ctrbpdint.b)) |> 
    #dplyr::filter(interaction==1) |> 
    #dplyr::filter(!str_detect(ref, regex("from beak to fruit", ignore_case = TRUE))) |> 
    dplyr::select(birds, pdint.b, pdint.b_log, NicheBreadth, NichePosition, Beak.Width,
                  Beak.Length_Culmen, Beak.Depth, Tarsus.Length, Wing.Length,
                  Kipps.Distance, hwi, Tail.Length, Mass, Beak.Width_log,
                  Range.Size, Mass_log, Range.Size_log, family, Realm,
                  human_impact,ED,DR,Wing.Length_log,ctrbpdint.b_log) |> 
    dplyr::rename(Family = family) |> 
    distinct() |> 
    drop_na()
  
  dim(data)
  
  # modelling
  
  # transform some variables
  
  data$NichePosition = abs(data$NichePosition)
  
  # scale variables that will be used in the modelling and exclude NAs
  data[,c("pdint.b","pdint.b_log","NicheBreadth", "NichePosition","Beak.Width",
          "Beak.Length_Culmen", "Beak.Depth", "Tarsus.Length", "Wing.Length",
          "Kipps.Distance","hwi" ,"Tail.Length", "Mass",
          "Range.Size","Mass_log","Range.Size_log",
          "Beak.Width_log","human_impact","ED","DR","Wing.Length_log",
          "ctrbpdint.b_log")] <- 
    scale(data[,c("pdint.b","pdint.b_log","NicheBreadth", "NichePosition","Beak.Width",
                  "Beak.Length_Culmen", "Beak.Depth", "Tarsus.Length", "Wing.Length",
                  "Kipps.Distance","hwi" ,"Tail.Length", "Mass",
                  "Range.Size","Mass_log","Range.Size_log",
                  "Beak.Width_log","human_impact","ED","DR","Wing.Length_log",
                  "ctrbpdint.b_log")])
  
  dataM <- as.data.frame(data)
  
  dataM$Family = as.factor(dataM$Family)
  dataM$Realm = as.factor(dataM$Realm)
  
  # Run MCMCglmm models metanetwork ######
  
  #use only one tree
  tree <-  phy[[i]]
  #tree <- phy
  tree$tip.label  <-  gsub("_"," ",tree$tip.label)
  
  # prepare data #####################################################
  # remove species not in the tree ####
  notintree = setdiff(dataM$birds,tree$tip.label)
  notindata = setdiff(tree$tip.label, dataM$birds)
  toexclude = c(notindata,notintree)
  
  #prune tree
  treeprune = drop.tip(tree,toexclude)
  
  #exclude species notintree from the dataset
  dataM = dataM[-which(dataM$birds %in% notintree),]
  
  # prepare the covariance matrix
  Ainv <- inverseA(treeprune)$Ainv
  
  # model without biogeographic realm
  prior<-list(R=list(V=1, nu=0.002),
              G=list(G1=list(V=1, nu=0.002),
                     G2=list(V=1, nu=0.002)))
  
  MCmodmeta <-  MCMCglmm(ctrbpdint.b_log~NicheBreadth+NichePosition+Beak.Width+Wing.Length+
                           Mass_log+ hwi+Range.Size_log+human_impact+ED,
                         random= ~birds + Realm,
                         prior = prior,
                         ginverse=list(birds=Ainv),
                         nitt=100000, burnin = 25000 ,thin = 50,
                         data=na.omit(dataM))
  
  
  list_posterioriC[[i]] = posterior.mode(MCmodmeta$Sol)
  list_intervalsC[[i]] = HPDinterval(MCmodmeta$Sol)
  list_summariesC[[i]]  = summary(MCmodmeta$Sol)[[2]]
  list_solC[[i]]        = MCmodmeta$Sol
  list_DICC[[i]]        = DIC(MCmodmeta)
  list_summaryGeralC[[i]] = summary(MCmodmeta)[[5]]
  list_randomC[[i]] = summary(MCmodmeta)[[6]]
  
} # fecha looping
toc()

#save products ----------------------

list_posterioriI_df <- map_df(list_posterioriC,~.x) |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind_cont","list_posteriori.csv"))

list_intervalsI_df <-map_df(list_intervalsC, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind_cont","list_intervals.csv"))

list_summariesI_df <- map_df(list_summariesC, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind_cont","list_summaries.csv"))

list_solI_df <- map_df(list_solC, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind_cont","list_sol.csv"))

list_dicI_df <- map_df(list_DICC, ~ {
  .x |> 
    as.data.frame() 
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  rename(DIC = ".x") |> 
  write_csv(here("productsR1","resu_mcmc_ind_cont","list_dic.csv"))

list_summaryI_df <- map_df(list_summaryGeralC, ~ {
  .x |> 
    as.data.frame() |> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind_cont","list_summaryGeral.csv"))

list_randomI_df <- map_df(list_randomC, ~ {
  .x |> 
    as.data.frame()|> 
    rownames_to_column("variable")
}, .id = "iteration") |> 
  mutate(model = "local") |> 
  write_csv(here("productsR1","resu_mcmc_ind_cont","list_random.csv"))



plot(data$)
