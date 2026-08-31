# PART 1: load datasets ---------------------------------------------------
##load packages -------------------------------------------------------------

library(tidyverse)
library(terra)
library(lme4)
library(MuMIn)
library(DHARMa)
library(merTools)
library(car)
library(sjPlot)
library(patchwork)
library(here)
library(clipr)
library(lmerTest)
library(car)
library(ggridges)

# Human impacts --------------------------------------------------------------------

#impacts per network
hu <- read_csv("Indnetwork_dataToAnalysis_v260826.csv")
hist(hu$human_impact)

hu |> 
summarise(below_50 = mean(human_impact < 0.5) * 100,
  above_50 = mean(human_impact > 0.5) * 100)
  
h1 <- hu |> 
  drop_na() |> 
  #ggplot(aes(x=human_impact,y=fct_reorder(Realm, human_impact, .fun = "mean",.desc = TRUE),fill=Realm))+
  ggplot(aes(x=human_impact,y=fct_rev(Realm),fill=Realm))+
  geom_density_ridges2(aes(fill=Realm),scale=2,color="gray")+
  theme_bw()+
  theme(text = element_text(size=16),
        panel.grid  = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())+
  labs(x="Human modification index")+
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2))+
  scale_fill_viridis_d()
h1

h2 <- hu |> 
  drop_na() |> 
  ggplot(aes(x=pdint.b_log,y=fct_rev(Realm),fill=Realm))+
  geom_density_ridges2(aes(fill=Realm),scale=2,color="gray")+
  theme_bw()+
  theme(text = element_text(size=16),
        panel.grid  = element_blank(),
        axis.title.y=element_blank(),
        #axis.ticks = element_blank(),
        #axis.text.y = element_blank()
        )+
  labs(x="log10(PDniche)")+
 # scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2))+
  scale_fill_viridis_d()
h2

h2 + h1+
  plot_annotation(tag_levels = "a", tag_suffix = ")")+
  plot_layout(guides = 'collect')
ggsave(here("productsR1","pdniche_impacts_realms.jpeg"),units="cm", width = 35,height = 15)

#impacts per taxonomy
ta <-  read_csv("Metanetwork_dataToAnalysis_v260826.csv")

ta |> 
  drop_na() |> 
  group_by(family) |> 
  summarise(mean = mean(human_impact), sd=sd(human_impact)) |> 
  arrange(desc(mean))

t1 <- ta |> 
  drop_na() |> 
  dplyr::filter(int==1) |> 
  #ggplot(aes(x=human_impact,y=fct_reorder(family, human_impact, .fun = "mean",.desc = TRUE)))+
  ggplot(aes(x=human_impact,y=fct_rev(family)))+
  geom_density_ridges2(scale=2,color="white",fill="gray25")+
  theme_bw()+
  theme(text = element_text(size=12),
        panel.grid  = element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.y = element_blank()
        )+
  labs(x="Human modification index")+
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2))+
  scale_fill_viridis_d()
t1

ta |> 
  drop_na() |> 
  mutate(pdint.b_log = log10(pdint.b)) |> 
  group_by(family) |> 
  summarise(mean = mean(pdint.b_log), sd=sd(pdint.b_log)) |> 
  arrange(desc(mean))

t2 <- ta |> 
  drop_na() |> 
  dplyr::filter(int==1) |> 
  
  mutate(pdint.b_log = log10(pdint.b)) |> 
  ggplot(aes(x=pdint.b_log,y=fct_rev(family)))+
  geom_density_ridges2(scale=2,color="white",fill="gray25")+
  theme_bw()+
  theme(text = element_text(size=12),
        panel.grid  = element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        #axis.text.y = element_blank()
        )+
  labs(x="log(PDniche)")+
  #scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2))+
  scale_fill_viridis_d()
t2

t2 + t1+
  plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave(here("productsR1","impacts_pdniceh_families.jpeg"),units="cm", width = 35,height = 30)


# correlation between metrics ---------------------------------------------

y1 <- nets_human_PD |> 
  ggplot(aes(pdint.b,ctrbpdint.b))+
  geom_point()+
  xlab("PDniche")+
  ylab("Contribution to plant PD")+
  theme_bw()+
  theme(text = element_text(size=12),
        panel.grid  = element_blank(),
        #axis.title.y=element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        #axis.text.y = element_blank()
  )+
  labs(title="Metanetwork")
y1

y2 <- inddata |> 
  ggplot(aes(pdint.b,ctrbpdint.b))+
  geom_point()+
  xlab("PDniche")+
  ylab("PDcont")+
  theme_bw()+
  theme(text = element_text(size=12),
        panel.grid  = element_blank(),
        #axis.title.y=element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        #axis.text.y = element_blank()
        )+
        labs(title="Local networks")
  
y2

cor.test(inddata$pdint.b,inddata$ctrbpdint.b)

y1 + y2+
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("productsR1", "correlations.jpeg"),units="cm", height = 10, width = 20)
