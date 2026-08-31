# Load packages -----------------------------------------------------------

library(patchwork)
library(ggplot2)
library(ggridges)
library(here)
library(tidyverse)
library(flextable)
library(ggExtra)

# analysis del DIC ####

dic <- read_csv(here("productsR1","resu_mcmc_meta","list_dic.csv")) |> 
  mutate(model = "Metanetworks")
min(dic$DIC)
max(dic$DIC)

dicN  <-  read_csv(here("productsR1","resu_mcmc_ind","list_dic.csv"))
min(dicN$DIC)
max(dicN$DIC)

dicC  <-  read_csv(here("productsR1","resu_mcmc_ind_cont","list_dic.csv"))
min(dicC$DIC)
max(dicC$DIC)


# p-values ######
list_summaryGeralM  <-  read_csv(here("productsR1","resu_mcmc_meta","list_summaryGeral.csv")) |> 
  mutate(model = "Metanetworks")
list_summaryGeralI <-  read_csv(here("productsR1","resu_mcmc_ind","list_summaryGeral.csv"))|> 
  mutate(model = "Local networks") 
list_summaryGeralC <-  read_csv(here("productsR1","resu_mcmc_ind_cont","list_summaryGeral.csv"))|> 
  mutate(model = "Local networks contribution") 

ps <- rbind(list_summaryGeralM,list_summaryGeralI,list_summaryGeralC )

#range of p values
rp <-  data.frame(aggregate(ps ,pMCMC~variable+model,range))
rp$pMCMC <- round(rp$pMCMC,3)
rp <- rp |> mutate(sig = ifelse(pMCMC<0.05, "**","n.s")) 

rp |> filter(model == "Metanetworks")
rp |> filter(model == "Local networks")
rp |> filter(model == "Local networks contribution")

# random variables effect sizes ####

list_randomM  <-  read_csv(here("productsR1","resu_mcmc_meta","list_random.csv"))
list_randomI  <-  read_csv(here("productsR1","resu_mcmc_ind","list_random.csv"))
list_randomC  <-  read_csv(here("productsR1","resu_mcmc_ind_cont","list_random.csv"))

aggregate(data=list_randomM, post.mean~variable,range)
aggregate(data=list_randomI, post.mean~variable,range)
aggregate(data=list_randomC, post.mean~variable,range)

# densitiy plots ##########################################

list_posteriori <-  read_csv(here("productsR1","resu_mcmc_meta","list_posteriori.csv")) |> 
  mutate(model = "Metanetworks") |> 
  rename(human_impact =Mdhuman_impact)
list_posterioriI <-  read_csv(here("productsR1","resu_mcmc_ind","list_posteriori.csv")) |> 
  mutate(model = "Local networks") 
list_posterioriC <-  read_csv(here("productsR1","resu_mcmc_ind_cont","list_posteriori.csv")) |> 
  mutate(model = "Local networks contribution") 

den <- rbind(list_posteriori,list_posterioriI,list_posterioriC)

den <-  den |> 
  rename(`Beak width`=Beak.Width,
        HWI= hwi,
        Mass = Mass_log, 
        `Climatic niche breadth`=NicheBreadth,
        `Climatic niche position` = NichePosition,
        `Range Size`= Range.Size_log, 
        'Wing length'=Wing.Length,
        `Human modification index` = human_impact) |> 
  pivot_longer(-model, names_to = "variable", values_to = "value") |> 
  mutate(variable = as.factor(variable)) 


Dmeta <- den |> filter(variable !="(Intercept)") |> 
  filter(model == "Metanetworks") |> 
  ggplot(aes(x=value,y=variable))+
  geom_density_ridges(scale = 1.3,fill = "grey", alpha = 0.5)+
  theme_bw()+
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.3)+
  theme(text = element_text(family = "sans"),
        panel.border = element_rect(fill = NA, colour = "black",size=0.1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size=14,colour = "black"),
        plot.subtitle =element_text(size=12,colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12,color="black"),
        axis.title = element_text(size=14,colour = "black"))+
  scale_x_continuous(limits = c(-0.3,1))+
  geom_vline(xintercept = 0,linetype="dashed",size=0.3)+
  xlab("Effect size")+
  ylab("")+
  scale_y_discrete(limits = rev,expand = expansion(mult = c(.05, .1)))+
  ggtitle("Meta-network - PDniche",subtitle = "DIC range =  3623 - 3654")+
  annotate("text",label="0.161 - 0.417",y=1,x=0.8,size=4)+ #Wl
  annotate("text",label="< 0.001",y=2,x=0.8,size=4)+ #range
  annotate("text",label="< 0.001"    ,y=3,x=0.8,size=4)+ #mass
  annotate("text",label="0.641 - 0.999",y=4,x=0.8,size=4)+ #hwi
  annotate("text",label="< 0.001",y=5,x=0.8,size=4)+ #hii
  annotate("text",label="0.071 - 0.176",y=6,x=0.8,size=4)+ #ED
  #annotate("text",label="0.001 - 0.007",y=7,x=0.8,size=4)+ #DR
  annotate("text",label="< 0.001",y=7,x=0.8,size=4)+        #NP
  annotate("text",label="0.648 - 0.856",y=8,x=0.8,size=4)+ #nb
  annotate("text",label="0.063 - 0.311",y=9,x=0.8,size=4)+ #beak wid
  annotate("text",label="pMCMC ranges",y=9.5,x=0.8,size=4)
Dmeta

Dind <- den |> filter(variable !="(Intercept)") |> 
  filter(model == "Local networks") |> 
ggplot(aes(x=value,y=variable))+
  geom_density_ridges(scale = 1.3,fill = "grey", alpha = 0.5)+
  theme_bw()+
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.3)+
  theme(text = element_text(family = "sans"),
        panel.border = element_rect(fill = NA, colour = "black",size=0.1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size=14,colour = "black"),
        plot.subtitle =element_text(size=12,colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12,color="black"),
        #axis.text.y = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title = element_text(size=14,colour = "black"))+
  scale_x_continuous(limits = c(-0.3,1))+
  geom_vline(xintercept = 0,linetype="dashed",size=0.3)+
  xlab("Effect size")+
  ylab("")+
  scale_y_discrete(limits = rev,expand = expansion(mult = c(.05, .1)))+
  ggtitle("Local network - PDniche",subtitle = "DIC range =  14963 - 15051")+
  annotate("text",label="0.457 - 0.921",y=1,x=0.8,size=4)+ #Wl
  annotate("text",label="< 0.001"      ,y=2,x=0.8,size=4)+ #range
  annotate("text",label="0.001 - 0.013",y=3,x=0.8,size=4)+ #mass
  annotate("text",label="0.001 - 0.008",y=4,x=0.8,size=4)+ #hwi
  annotate("text",label="0.027 - 0.084"      ,y=5,x=0.8,size=4)+ #hii
  annotate("text",label="0.397 - 0.787",y=6,x=0.8,size=4)+ #ED
  #annotate("text",label="0.007 - 0.104",y=7,x=0.8,size=4)+ #DR
  annotate("text",label="< 0.001"       ,y=7,x=0.8,size=4)+ #NP
  annotate("text",label="0.001 - 0.007",y=8,x=0.8,size=4)+ #nb
  annotate("text",label="0.589 - 1.000",y=9,x=0.8,size=4)+ #beak wid
  annotate("text",label="pMCMC ranges", y=9.5,x=0.8,size=4)
Dind

Dc <- den |> filter(variable !="(Intercept)") |> 
  filter(model == "Local networks contribution") |> 
  ggplot(aes(x=value,y=variable))+
  geom_density_ridges(scale = 1.3,fill = "grey", alpha = 0.5)+
  theme_bw()+
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.3)+
  theme(text = element_text(family = "sans"),
        panel.border = element_rect(fill = NA, colour = "black",size=0.1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size=14,colour = "black"),
        plot.subtitle =element_text(size=12,colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12,color="black"),
        axis.text.y = element_blank(),
        #axis.text = element_text(size=12,color="black"),
        axis.title = element_text(size=14,colour = "black"))+
  scale_x_continuous(limits = c(-0.3,1))+
  geom_vline(xintercept = 0,linetype="dashed",size=0.3)+
  xlab("Effect size")+
  ylab("")+
  scale_y_discrete(limits = rev,expand = expansion(mult = c(.05, .1)))+
  ggtitle("Local Network - PDcontrib",subtitle = "DIC range =  14602 - 14659")+
  annotate("text",label="0.149 - 0.365",y=1,x=0.8,size=4)+ #Wl
  annotate("text",label="0.001 - 0.003",y=2,x=0.8,size=4)+ #range
  annotate("text",label="0.001 - 0.027" ,y=3,x=0.8,size=4)+ #mass
  annotate("text",label="0.001 - 0.017",y=4,x=0.8,size=4)+ #hwi
  annotate("text",label="< 0.001"      ,y=5,x=0.8,size=4)+ #hii
  annotate("text",label="0.184 - 0.511",y=6,x=0.8,size=4)+ #ED
  #annotate("text",label="0.001 - 0.008",y=7,x=0.8,size=4)+ #DR
  annotate("text",label="< 0.001"      ,y=7,x=0.8,size=4)+ #NP
  annotate("text",label="0.293 - 0.667",y=8,x=0.8,size=4)+ #nb
  annotate("text",label="0.180 - 0.519",y=9,x=0.8,size=4)+ #beak wid
  annotate("text",label="pMCMC ranges",y=9.5,x=0.8,size=4)
Dc

#### save composite file -------------------------------------
if(F){
Dmeta + Dind+
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("productsR1","density_plots.jpeg"),units = "cm",width = 25,height = 12,dpi = 300)


Dc
ggsave(here("productsR1","density_plots_contribution.jpeg"),units = "cm",width = 15,height = 12,dpi = 300)
}

Dind+ Dc+Dmeta + 
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("productsR1","density_plots_ALL_R1.jpeg"),units = "cm",
       width = 33,height = 12,dpi = 600)

# Regression plots --------------------------------------------------------

### load data metanetwork ####
data <- read_csv("Metanetwork_dataToAnalysis_v260826.csv")
data <- data |> 
  mutate(Mass_log = log10(Mass),
         Range.Size_log = log10(Range.Size),
         pdint.b_log = log10(pdint.b),
         Beak.Width_log = log10(Beak.Width),
         Wing.Length_log = log10(Wing.Length)) |> 
  dplyr::filter(int==1) |> 
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

# plot human impact

coefsHI = as.data.frame(matrix(NA, ncol=2, nrow = 100))
colnames(coefsHI) = c("a","b")
estHI=list()
temp <- list_summaryGeralM |> filter(variable == "Mdhuman_impact" | variable == "(Intercept)" )

plot(dataM$Mdhuman_impact,dataM$pdint.b_log,xlab="Human modification index",ylab="PDNiche", pch=19,
     cex.lab=1,cex.axis=1,cex=0.5)

for (i in unique(temp$iteration)){
  coefsHI[i,] <-  temp[temp$iteration==i,c(2,3)] |> 
    pivot_wider(names_from = variable,values_from = post.mean) |> 
    mutate(across(everything(), as.numeric))
  estHI[[i]] <-  coefsHI[i,1] + coefsHI[i,2]*dataM$Mdhuman_impact
  lines(dataM$Mdhuman_impact,estHI[[i]],  col="darkgray")


# Create a dataframe to store all predictions
all_predictions <- data.frame()

for (i in unique(temp$iteration)) {
  coefsHI[i,] <- temp[temp$iteration == i, c(2,3)] |> 
    pivot_wider(names_from = variable, values_from = post.mean) |> 
    mutate(across(everything(), as.numeric))
  
  estHI[[i]] <- coefsHI[i,1] + coefsHI[i,2] * dataM$Mdhuman_impact
  
  # Store predictions for each iteration
  temp_pred <- data.frame(
    Mdhuman_impact = dataM$Mdhuman_impact,
    prediction = estHI[[i]],
    iteration = i
  )
  all_predictions <- rbind(all_predictions, temp_pred)
}

# Create the ggplot
rpm <- ggplot() +
  geom_point(data = dataM, 
             aes(x = Mdhuman_impact, y = pdint.b_log), alpha=0.2) +
  geom_line(data = all_predictions,
            aes(x = Mdhuman_impact, y = prediction, group = iteration),
            color = "darkgray", alpha = 0.3,linewidth=0.1) +
  labs(x = "Human modification index",
       y = "PDNiche") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        panel.border = element_rect(fill = NA, colour = "black",size=0.1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=14,color="black"),
         axis.title = element_text(size=14,colour = "black"))+
  ggtitle("Meta-network - PDniche")
rpm

rpm <- ggMarginal(rpm, type = "density",fill="darkgray")

### load local networks ################################

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

# transform some variables

data$NichePosition = abs(data$NichePosition)

# scale variables that will be used in the modelling and exclude NAs
data[,c("pdint.b","pdint.b_log","NicheBreadth", "NichePosition","Beak.Width",
        "Beak.Length_Culmen", "Beak.Depth", "Tarsus.Length", "Wing.Length",
        "Kipps.Distance","hwi" ,"Tail.Length", "Mass",
        "Range.Size","Mass_log","Range.Size_log","ctrbpdint.b_log",
        "Beak.Width_log","human_impact","ED","DR","Wing.Length_log")] <- 
  scale(data[,c("pdint.b","pdint.b_log","NicheBreadth", "NichePosition","Beak.Width",
                "Beak.Length_Culmen", "Beak.Depth", "Tarsus.Length", "Wing.Length",
                "Kipps.Distance","hwi" ,"Tail.Length", "Mass",
                "Range.Size","Mass_log","Range.Size_log","ctrbpdint.b_log",
                "Beak.Width_log","human_impact","ED","DR","Wing.Length_log")])

dataM <- as.data.frame(data)

dataM$Family = as.factor(dataM$Family)
dataM$Realm = as.factor(dataM$Realm)

# plot human impact

coefsHI = as.data.frame(matrix(NA, ncol=2, nrow = 100))
colnames(coefsHI) = c("a","b")
estHI=list()
temp <- list_summaryGeralI |> filter(variable == "human_impact" | variable == "(Intercept)" )

# Create a dataframe to store all predictions
all_predictions <- data.frame()

for (i in unique(temp$iteration)) {
  coefsHI[i,] <- temp[temp$iteration == i, c(2,3)] |> 
    pivot_wider(names_from = variable, values_from = post.mean) |> 
    mutate(across(everything(), as.numeric))
  
  estHI[[i]] <- coefsHI[i,1] + coefsHI[i,2] * dataM$human_impact
  
  # Store predictions for each iteration
  temp_pred <- data.frame(
    human_impact = dataM$human_impact,
    prediction = estHI[[i]],
    iteration = i
  )
  all_predictions <- rbind(all_predictions, temp_pred)
}

# Create the ggplot
rpl <- ggplot() +
  geom_point(data = dataM, 
             aes(x = human_impact, y = pdint.b_log),alpha=0.2) +
  geom_line(data = all_predictions,
            aes(x = human_impact, y = prediction, group = iteration),
            color = "darkgray", alpha = 0.3,linewidth=0.1) +
  labs(x = "Human modification index",
       y = "PDNiche") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        panel.border = element_rect(fill = NA, colour = "black",size=0.1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=14,colour = "black"))+
  ggtitle("Loca network - PDniche")
rpl

rpl <- ggMarginal(rpl, type = "density",fill="darkgray")

# plot human impact for birds contributios

coefsHI = as.data.frame(matrix(NA, ncol=2, nrow = 100))
colnames(coefsHI) = c("a","b")
estHI=list()
temp <- list_summaryGeralC |> filter(variable == "human_impact" | variable == "(Intercept)" )

# Create a dataframe to store all predictions
all_predictions <- data.frame()

for (i in unique(temp$iteration)) {
  coefsHI[i,] <- temp[temp$iteration == i, c(2,3)] |> 
    pivot_wider(names_from = variable, values_from = post.mean) |> 
    mutate(across(everything(), as.numeric))
  
  estHI[[i]] <- coefsHI[i,1] + coefsHI[i,2] * dataM$human_impact
  
  # Store predictions for each iteration
  temp_pred <- data.frame(
    human_impact = dataM$human_impact,
    prediction = estHI[[i]],
    iteration = i
  )
  all_predictions <- rbind(all_predictions, temp_pred)
}

# Create the ggplot
rplc <- ggplot() +
  geom_point(data = dataM, 
             aes(x = human_impact, y = ctrbpdint.b_log),
             alpha=0.2) +
  geom_line(data = all_predictions,
            aes(x = human_impact, y = prediction, group = iteration),
            color = "darkgray", alpha = 0.3,linewidth=0.1) +
  labs(x = "Human modification index",
       y = "PDcont") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        panel.border = element_rect(fill = NA, colour = "black",size=0.1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=14,colour = "black"))+
  ggtitle("Loca network - PDcontrib")
rplc

rplc <- ggMarginal(rplc, type = "density",fill="darkgray")
rplc

#### save composite graphic ---------------------------------------

rpm + rpl +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("productsR1","regression_plots.jpeg"),units = "cm",width = 28,height = 12,dpi = 300)

wrap_elements(rpl) + wrap_elements(rplc)+ wrap_elements(rpm) +  
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("productsR1","regression_plots_contribution_R1.png"),units = "cm",
       width = 38,height = 12,dpi = 600)


Dc + rplc +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("productsR1","contribution_plots.jpeg"),units = "cm",width = 28,height = 12,dpi = 300)

rplc
ggsave(here("productsR1","contribution_Regressionplots.jpeg"),units = "cm",width =15,height = 10,dpi = 300)

# Get 95% values #############
Mmodel95 <-  cbind(
  aggregate(data=list_summaryGeralM,
            post.mean~variable,mean),
  aggregate(data=list_summaryGeralM,
            `l-95% CI`~variable,mean),
  aggregate(data=list_summaryGeralM,
            `u-95% CI`~variable,mean)
)
Mmodel95 = Mmodel95[,c(1,2,4,6)]
Mmodel95[,c(2,3,4)] = round(Mmodel95[,c(2,3,4)],2)
Mmodel95$range = paste(Mmodel95$post.mean,"(",Mmodel95$`l-95% CI`,"-",Mmodel95$`u-95% CI`,")",sep=" ")


Imodel95 = cbind(
  aggregate(data=list_summaryGeralI,
            post.mean~variable,mean),
  aggregate(data=list_summaryGeralI,
            `l-95% CI`~variable,mean),
  aggregate(data=list_summaryGeralI,
            `u-95% CI`~variable,mean)
)
Imodel95 = Imodel95[,c(1,2,4,6)]
Imodel95[,c(2,3,4)] = round(Imodel95[,c(2,3,4)],2)
Imodel95$range = paste(Imodel95$post.mean,"(",Imodel95$`l-95% CI`,"-",Imodel95$`u-95% CI`,")",sep=" ")

Cmodel95 = cbind(
  aggregate(data=list_summaryGeralC,
            post.mean~variable,mean),
  aggregate(data=list_summaryGeralC,
            `l-95% CI`~variable,mean),
  aggregate(data=list_summaryGeralC,
            `u-95% CI`~variable,mean)
)
Cmodel95 = Cmodel95[,c(1,2,4,6)]
Cmodel95[,c(2,3,4)] = round(Cmodel95[,c(2,3,4)],2)
Cmodel95$range = paste(Cmodel95$post.mean,"(",Cmodel95$`l-95% CI`,"-",Cmodel95$`u-95% CI`,")",sep=" ")


## Export as table ---------------------------
Mmodel <- Mmodel95[,c(1,5)]
Mmodel$model <- "Meta-network - PDNiche"

Imodel <- Imodel95[,c(1,5)]
Imodel$model <- "Local network - PDNiche"

Cmodel <- Cmodel95[,c(1,5)]
Cmodel$model <- "Local network - PDCont"

models <- cbind(Mmodel, Imodel, Cmodel)

ft <- flextable(models) %>%
  theme_vanilla() %>%
  autofit() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all")

# Save to Word
save_as_docx(ft,path = here("productsR1","model_output_table.docx"))

