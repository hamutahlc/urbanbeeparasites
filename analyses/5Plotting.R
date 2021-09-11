
setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())

library(ggplot2)
library(viridis)
library(brms)
library(bayesplot)
library(tidybayes)
library(tidyverse)
library(gridExtra)
library(grid)

load("../data/specimens-complete.Rdata")
load("~/Dropbox/urbanbeeparasites/data/CommunityHealthresults.Rdata")

## **************************************************************
## what this do
## **************************************************************

stanplot(fit.parpathRich.apis,
         type = "areas",
         prob = 0.95)

stanplot(fit.parpathRich.bombus,
         type = "areas",
         prob = 0.95)


## **************************************************************
## parasite/path richness in bombus: bee diversity
## **************************************************************

p1.parasite  <- fit.parpathRich.bombus %>%
  spread_draws(b_BeeDiversity_Intercept,
               b_ParPathRichness_Intercept,
               b_ParPathRichness_BeeDiversity) %>%
  mutate(BeeDiversity = 
           list(seq(range(par.and.path$BeeDiversity)[1],
                    range(par.and.path$BeeDiversity)[2] + 0.55,
                    0.01))) %>%
  unnest(BeeDiversity) %>%
  mutate(pred = exp(b_BeeDiversity_Intercept +
                      b_ParPathRichness_Intercept +
                      b_ParPathRichness_BeeDiversity*BeeDiversity)/
           (1+exp(b_BeeDiversity_Intercept +
                    b_ParPathRichness_Intercept +  
                    b_ParPathRichness_BeeDiversity*BeeDiversity))) %>%
  group_by(BeeDiversity) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = BeeDiversity, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="blue4") +
  ylab("Parasite and Pathogen Richness") +
  xlab("Bee Community Diversity") +
  ylim(0,1)  +
  xlim(range(par.and.path$BeeDiversity))  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_point(data=par.and.path[par.and.path$Weights == 1,],
             aes(y=ParPathRichness, x=BeeDiversity))

ggsave("figures/parpath_apis_beeDiversity.pdf",
       height=4, width=5)
  
            


## **************************************************************
## parasite/path richness in bombus: abund perennials
## **************************************************************


#cgane to binomial steps

p2.parasite  <- fit.parpathRich.bombus %>%
  spread_draws(b_ParPathRichness_Intercept,
               b_ParPathRichness_AbundWoodyFlowers) %>%
  mutate(AbundWoodyFlowers = 
           list(seq(range(par.and.path$AbundWoodyFlowers)[1],
                    range(par.and.path$AbundWoodyFlowers)[2],
                    0.01))) %>%
  unnest(AbundWoodyFlowers) %>%
  mutate(pred = exp(b_ParPathRichness_Intercept +
                      b_ParPathRichness_AbundWoodyFlowers*AbundWoodyFlowers)) %>%
  group_by(AbundWoodyFlowers) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = AbundWoodyFlowers, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="blue4") +
  ylab("Parasite and Pathogen Richness") +
  xlab("Abundance Flowering Perennials") +
  ylim(0,1)  +
  xlim(range(par.and.path$AbundWoodyFlowers))  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_point(data=par.and.path[par.and.path$Weights == 1,],
             aes(y=ParPathRichness, x=AbundWoodyFlowers))


ggsave("figures/parpath_apis_AbundWoodyFlowers.pdf",
       height=4, width=5)




## **************************************************************
## parasite/path richness in bombus: garden size
## **************************************************************


## **************************************************************
## parasite/path richness in bombus: apis parasite rates
## **************************************************************


## **************************************************************
## panels combine
## **************************************************************

parasite.apis.all <- grid.arrange(p1.parasite, p2.parasite, p3.parasite, p4.parasite, ncol=2)

ggsave(parasite.apis.all, file="figures/parpath_apis.pdf",
       height=4, width=10)

