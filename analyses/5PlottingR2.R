
setwd("~/Drobbox/urbanbeeparasites")
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
library(scales)
source("src/makeMultiLevelData.R")
source("src/misc.R")


#****************extra functions*********
#*
#*
pdf.f <- function(f, file, ...) {
cat(sprintf('Writing %s\n', file))
pdf(file, ...)
on.exit(dev.off())
f()
}

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3],
              alpha=alpha))
}


standardize <- function(x)
  (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

unstandardize <- function(x, orig){
  (x*sd(orig, na.rm=TRUE)) + mean(orig, na.rm=TRUE)
}


standardize.axis <- function(x, orig)
  (x-mean(orig, na.rm=TRUE))/sd(orig, na.rm=TRUE)
## **************************************************************

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
## plotting, unscaling labs
## ***********************************************************************

par.and.path<- par.and.path[!is.na(par.and.path$ParPathRichness),]

labs.BeeAbund <- (pretty(c(0, par.and.path$BeeAbund), n=5))
axis.BeeAbund <-  standardize.axis(labs.BeeAbund, par.and.path$BeeAbund)

labs.BeeDiversity <- (pretty(c(0, par.and.path$BeeDiversity), n=5))
axis.BeeDiversity <-  standardize.axis(labs.BeeDiversity, par.and.path$BeeDiversity)

labs.natural1000m <- (pretty(c(0, par.and.path$natural1000m), n=5))
axis.natural1000m <-  standardize.axis(labs.natural1000m, par.and.path$natural1000m)

labs.AbundWoodyFlowers <- (pretty(c(0, par.and.path$AbundWoodyFlowers), n=5))
axis.AbundWoodyFlowers <-  standardize.axis(labs.AbundWoodyFlowers, par.and.path$AbundWoodyFlowers)

labs.Size <- (pretty(c(0, par.and.path$Size), n=5))
axis.Size <-  standardize.axis(labs.Size, par.and.path$Size)

labs.apis.parpath.rate <- (pretty(c(0, par.and.path$apis.parpath.rate), n=5))
axis.apis.parpath.rate <-  standardize.axis(labs.apis.parpath.rate, par.and.path$apis.parpath.rate)

labs.ParPathRichness <- (pretty(c(0, par.and.path$ParPathRichness), n=4))
axis.ParPathRichness <-  standardize.axis(labs.ParPathRichness, par.and.path$ParPathRichness)

#**********************************************************
# make sure variables in the data are scaled
all.indiv.vars <- c("BeeDiversity", "BeeAbund",
                    "natural1000m",
                    "AbundWoodyFlowers",
                    "AbundAnnualFlowers",
                    "Size","ParPathRichness")
xvar.par.path.apis <- c(all.indiv.vars,
                        "bombus.parpath.rate")
xvar.par.path.bombus <- c(all.indiv.vars,
                          "apis.parpath.rate")

# scale x variables

par.and.path[, c(all.indiv.vars, "apis.parpath.rate", "bombus.parpath.rate" )] <-
  apply(par.and.path[, c(all.indiv.vars, "apis.parpath.rate", "bombus.parpath.rate" )], 2, scale)
## ***********************************************************************

## **************************************************************
## PARASITE AND PATHOGEN RICHNESS FIGURES 

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
            pred_low_90 = quantile(pred, prob = 0.05),
            pred_high_90 = quantile(pred, prob = 0.95),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = BeeDiversity, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="blue4") +
  geom_ribbon(aes(ymin=pred_low_90, ymax = pred_high_90), alpha = 0.2,
              fill="blue4") +
  geom_ribbon(aes(ymin=pred_low_85, ymax = pred_high_85), alpha = 0.2,
              fill="blue4") +
  ylab("Parasite & Pathogen Richness \nin Bumble Bees") +
  xlab("Bee Community Diversity") + 
   scale_x_continuous(
     breaks = axis.BeeDiversity,
     labels =  labs.BeeDiversity) +
  #   scale_y_continuous(
  #     breaks = axis.ParPathRichness,
  #     labels =  labs.ParPathRichness) +
  # coord_cartesian(ylim = range(axis.ParPathRichness)) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=14)) +
  theme_classic()  +  
  y(0,1)


  
ggsave("figuresR2/bothaxes/parpath_bombus_beeDiversity.pdf",
       height=4, width=5)
  


## **************************************************************
## parasite/path richness in bombus: abund perennials
## **************************************************************


p2.parasite  <- fit.parpathRich.bombus %>%
  spread_draws(b_ParPathRichness_Intercept,
               b_ParPathRichness_AbundWoodyFlowers) %>%
  mutate(AbundWoodyFlowers = 
           list(seq(range(par.and.path$AbundWoodyFlowers)[1],
                    range(par.and.path$AbundWoodyFlowers)[2],
                    0.01))) %>%
  unnest(AbundWoodyFlowers) %>%
  mutate(pred = exp(b_ParPathRichness_Intercept +
                      b_ParPathRichness_AbundWoodyFlowers*AbundWoodyFlowers)/
           (1+exp(b_ParPathRichness_Intercept +
                    b_ParPathRichness_AbundWoodyFlowers*AbundWoodyFlowers))) %>%
  group_by(AbundWoodyFlowers) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_90 = quantile(pred, prob = 0.05),
            pred_high_90 = quantile(pred, prob = 0.95),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = AbundWoodyFlowers, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="blue4") +
  geom_ribbon(aes(ymin=pred_low_90, ymax = pred_high_90), alpha = 0.2,
              fill="blue4") +
  geom_ribbon(aes(ymin=pred_low_85, ymax = pred_high_85), alpha = 0.2,
              fill="blue4") +
  ylab(" ") +
  xlab("Abundance Flowering Perennials") +
  scale_x_continuous(
    breaks = axis.AbundWoodyFlowers,
    labels =  labs.AbundWoodyFlowers) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=14)) +
  theme_classic() +
  y(0,1)



ggsave("figuresR2/bothaxes/parpath_bombus_AbundWoodyFlowers.pdf",
       height=4, width=5)



## **************************************************************
## parasite/path richness in bombus: garden size
## **************************************************************


p3.parasite  <- fit.parpathRich.bombus %>%
  spread_draws(b_ParPathRichness_Intercept,
               b_ParPathRichness_Size) %>%
  mutate(Size = 
           list(seq(range(par.and.path$Size)[1],
                    range(par.and.path$Size)[2],
                    0.01))) %>%
  unnest(Size) %>%
  mutate(pred = exp(b_ParPathRichness_Intercept +
                      b_ParPathRichness_Size*Size)/
           (1+exp(b_ParPathRichness_Intercept +
                    b_ParPathRichness_Size*Size))) %>%
  group_by(Size) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_90 = quantile(pred, prob = 0.05),
            pred_high_90 = quantile(pred, prob = 0.95),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = Size, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="blue4") +
  geom_ribbon(aes(ymin=pred_low_90, ymax = pred_high_90), alpha = 0.2,
              fill="blue4") +
  geom_ribbon(aes(ymin=pred_low_85, ymax = pred_high_85), alpha = 0.2,
              fill="blue4") +
  ylab("Parasite & Pathogen Richness \nin Bumble Bees") +
  xlab("Garden Size") +
  scale_x_continuous(
    breaks = axis.Size,
    labels =  labs.Size) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=14)) +
  theme_classic() +
  y(0,1)



ggsave("figuresR2/bothaxes/parpath_bombus_Size.pdf",
       height=4, width=5)




## **************************************************************
## parasite/path richness in bombus: apis parasite rates
## **************************************************************


p4.parasite  <- fit.parpathRich.bombus %>%
  spread_draws(b_ParPathRichness_Intercept,
               b_ParPathRichness_apis.parpath.rate) %>%
  mutate(apis.parpath.rate = 
           list(seq(range(par.and.path$apis.parpath.rate)[1],
                    range(par.and.path$apis.parpath.rate)[2],
                    0.01))) %>%
  unnest(apis.parpath.rate) %>%
  mutate(pred = exp(b_ParPathRichness_Intercept +
                      b_ParPathRichness_apis.parpath.rate*apis.parpath.rate)/
           (1+exp(b_ParPathRichness_Intercept +
                    b_ParPathRichness_apis.parpath.rate*apis.parpath.rate))) %>%
  group_by(apis.parpath.rate) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_90 = quantile(pred, prob = 0.05),
            pred_high_90 = quantile(pred, prob = 0.95),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = apis.parpath.rate, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="blue4") +
  geom_ribbon(aes(ymin=pred_low_90, ymax = pred_high_90), alpha = 0.2,
              fill="blue4") +
  geom_ribbon(aes(ymin=pred_low_85, ymax = pred_high_85), alpha = 0.2,
              fill="blue4") +
  ylab(" ") +
  xlab("Honey Bee Parasitism Rate") +
  scale_x_continuous(
    breaks = axis.apis.parpath.rate ,
    labels =  labs.apis.parpath.rate ) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=14)) +
  theme_classic() +
  y(0,1)

ggsave("figuresR2/bothaxes/parpath_bombus_apis.parpath.rate.pdf",
       height=4, width=5)

## **************************************************************
## panels combine
## **************************************************************

parasite.bombus.all <- grid.arrange(p1.parasite, p2.parasite, p3.parasite, p4.parasite, ncol=2)

ggsave(parasite.bombus.all, file="figuresR2/bothaxes/parpath_bombus.pdf",
       height=8, width=10)


## **************************************************************
## BEE DIVERSITY AND BEE ABUNDANCE FIGURES 
## **************************************************************

## make sure variables in the figures are scaled
# all.indiv.vars <- c("BeeDiversity", "BeeAbund",
#                     "natural1000m",
#                     "AbundWoodyFlowers",
#                     "AbundAnnualFlowers",
#                     "Size")
# xvar.par.path.apis <- c(all.indiv.vars,
#                         "bombus.parpath.rate")
# xvar.par.path.bombus <- c(all.indiv.vars,
#                           "apis.parpath.rate")

# scale x variables

#par.and.path[, c(all.indiv.vars, "apis.parpath.rate", "bombus.parpath.rate" )] <-
 # apply(par.and.path[, c(all.indiv.vars, "apis.parpath.rate", "bombus.parpath.rate" )], 2, scale)

## **************************************************************
## bee diversity: garden size
## **************************************************************


p1.beediv  <- fit.parpathRich.bombus %>%
  spread_draws(b_BeeDiversity_Intercept,
               b_BeeDiversity_Size) %>%
  mutate(Size = 
           list(seq(min(par.and.path$Size),
                    max(par.and.path$Size),
                    0.01))) %>%
  unnest(Size) %>%
  mutate(pred = exp(b_BeeDiversity_Intercept +
                       b_BeeDiversity_Size*Size)/
            (1+exp(b_BeeDiversity_Intercept +
                     b_BeeDiversity_Size*Size))) %>%
  group_by(Size) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_90 = quantile(pred, prob = 0.05),
            pred_high_90 = quantile(pred, prob = 0.95),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = Size, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="goldenrod3") +
  geom_ribbon(aes(ymin=pred_low_90, ymax = pred_high_90), alpha = 0.2,
              fill="goldenrod3") +
  geom_ribbon(aes(ymin=pred_low_85, ymax = pred_high_85), alpha = 0.2,
              fill="goldenrod3") +
  ylab("Bee Diversity") +
  xlab("Garden Size") +
  scale_x_continuous(
    breaks = axis.Size ,
    labels =  labs.Size ) +
  # scale_y_continuous(
  #   breaks = axis.BeeDiversity,
  #   labels =  labs.BeeDiversity) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=14)) +
  theme_classic() 

ggsave("figuresR2/BeeDiversity_Size.pdf",
       height=4, width=5)

## **************************************************************
## bee diversity: abund perennials
## **************************************************************

p2.beediv  <- fit.parpathRich.bombus %>%
  spread_draws(b_BeeDiversity_Intercept,
               b_BeeDiversity_AbundWoodyFlowers) %>%
  mutate(AbundWoodyFlowers = 
           list(seq(range(par.and.path$AbundWoodyFlowers)[1],
                    range(par.and.path$AbundWoodyFlowers)[2],
                    0.01))) %>%
  unnest(AbundWoodyFlowers) %>%
  mutate(pred = exp(b_BeeDiversity_Intercept +
                      b_BeeDiversity_AbundWoodyFlowers*AbundWoodyFlowers)/
           (1+exp(b_BeeDiversity_Intercept +
                    b_BeeDiversity_AbundWoodyFlowers*AbundWoodyFlowers))) %>%
  group_by(AbundWoodyFlowers) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_90 = quantile(pred, prob = 0.05),
            pred_high_90 = quantile(pred, prob = 0.95),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = AbundWoodyFlowers, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="goldenrod3") +
  geom_ribbon(aes(ymin=pred_low_90, ymax = pred_high_90), alpha = 0.2,
              fill="goldenrod3") +
  geom_ribbon(aes(ymin=pred_low_85, ymax = pred_high_85), alpha = 0.2,
              fill="goldenrod3") +
  ylab("Bee Diversity") +
  xlab("Abundance Flowering Perennials") +
  scale_x_continuous(
    breaks = axis.AbundWoodyFlowers ,
    labels =  labs.AbundWoodyFlowers ) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=14)) +
  theme_classic() +
  ylim(0,1)


ggsave("figuresR2/BeeDiversity_AbundWoodyFlowers.pdf",
       height=4, width=5)

## **************************************************************
## bee abundance: natural habitat
## **************************************************************


p3.beeabund  <- fit.parpathRich.bombus %>%
  spread_draws(b_BeeAbund_Intercept,
               b_BeeAbund_natural1000m) %>%
  mutate(natural1000m = 
           list(seq(range(par.and.path$natural1000m)[1],
                    range(par.and.path$natural1000m)[2],
                    0.01))) %>%
  unnest(natural1000m) %>%
  mutate(pred = exp(b_BeeAbund_Intercept +
                      b_BeeAbund_natural1000m*natural1000m)/
           (1+exp(b_BeeAbund_Intercept +
                    b_BeeAbund_natural1000m*natural1000m))) %>%
  group_by(natural1000m) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_90 = quantile(pred, prob = 0.05),
            pred_high_90 = quantile(pred, prob = 0.95),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = natural1000m, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="goldenrod3") +
  geom_ribbon(aes(ymin=pred_low_90, ymax = pred_high_90), alpha = 0.2,
              fill="goldenrod3") +
  geom_ribbon(aes(ymin=pred_low_85, ymax = pred_high_85), alpha = 0.2,
              fill="goldenrod3") +
  ylab("Bee Abundance") +
  xlab("Natural Habitat 1km") +
  scale_x_continuous(
    breaks = axis.natural1000m ,
    labels =  labs.natural1000m ) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=14)) +
  theme_classic() +
  ylim(0,1)


ggsave("figuresR2/BeeAbund_natural1000m.pdf",
       height=4, width=5)
## **************************************************************
## bee abundance: garden size
## **************************************************************

p4.beeabund  <- fit.parpathRich.bombus %>%
  spread_draws(b_BeeAbund_Intercept,
               b_BeeAbund_Size) %>%
  mutate(Size = 
           list(seq(range(par.and.path$Size)[1],
                    range(par.and.path$Size)[2],
                    0.01))) %>%
  unnest(Size) %>%
  mutate(pred = exp(b_BeeAbund_Intercept +
                      b_BeeAbund_Size*Size)/
           (1+exp(b_BeeAbund_Intercept +
                    b_BeeAbund_Size*Size))) %>%
  group_by(Size) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low_95 = quantile(pred, prob = 0.025),
            pred_high_95 = quantile(pred, prob = 0.975),
            pred_low_90 = quantile(pred, prob = 0.05),
            pred_high_90 = quantile(pred, prob = 0.95),
            pred_low_85 = quantile(pred, prob = 0.075),
            pred_high_85 = quantile(pred, prob = 0.925)) %>%
  ggplot(aes(x = Size, y=pred_m)) +
  geom_line() +
  geom_ribbon(aes(ymin=pred_low_95, ymax = pred_high_95), alpha = 0.2,
              fill="goldenrod3") +
  geom_ribbon(aes(ymin=pred_low_90, ymax = pred_high_90), alpha = 0.2,
              fill="goldenrod3") +
  geom_ribbon(aes(ymin=pred_low_85, ymax = pred_high_85), alpha = 0.2,
              fill="goldenrod3") +
  ylab("Bee Abundance") +
  xlab("Garden Size") +
  scale_x_continuous(
    breaks = axis.Size,
    labels =  labs.Size) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=14)) +
  theme_classic() +
  ylim(0,1)


ggsave("figuresR2/BeeAbund_Size.pdf",
       height=4, width=5)


## **************************************************************
## panels combine
## **************************************************************

beeabund.and.div <- grid.arrange(p1.beediv, p2.beediv, p4.beeabund, p3.beeabund, ncol=2)

ggsave(beeabund.and.div, file="figuresr2/beeabund_div.pdf",
       height=8, width=10)





## **************************************************************
## Parasite Summary Figure
## **************************************************************


load("../data/specimens-complete.Rdata")
spec <- par.and.path
parasites <- c("Apicystis", "Phorid", "Crithidia",
               "CBPV", "DWV_KV_VDV",
               "ABPV_KBV_IAPV", "BQCV", "SBPV", "SBV")

# summarize parasite data by 
makeStructParasite <- function(spec, parasites){
  parasite.pre.comm <- spec[, c("Site", parasites)]
  
  parasite.pre.comm <- parasite.pre.comm  %>%
    group_by(Site) %>%
    summarise_each(list(mean))
  
  
  Genus <- spec$Genus[match(parasite.pre.comm$Site,
                            spec$Site)]
  
  sites <- parasite.pre.comm$Site
  
  comm <- parasite.pre.comm
  comm$Site <- NULL
  comm <- as.matrix(comm)
  rownames(comm) <- parasite.pre.comm$Site
  list(comm=comm,
       ## sites=sites,
       ## site.type=site.type,
       Genus = Genus)
}

## make parasite matrix
parasite.comm <- makeStructParasite(spec, parasites)


## make a long table to plot parasitism by site type,parasite species
parasite.pre.long <- as.data.frame(parasite.comm$comm)
parasite.pre.long$Site <- rownames(parasite.pre.long)
parasite.long <- pivot_longer(parasite.pre.long, cols=parasites,
                              names_to = "Parasite",
                              values_to = "Parasitism")
parasite.long$Genus <- spec$Genus[match(parasite.long$Site,
                                        spec$Site)]

# make graph
parasite.long$Parasite <- factor(parasite.long$Parasite, levels=c("Apicystis", "Crithidia", "Phorid",
                                                                  "CBPV","SBPV", "SBV", "ABPV_KBV_IAPV","BQCV", "DWV_KV_VDV"))

p <- ggplot(parasite.long, aes(x=Parasite, y=Parasitism,
                               fill=Genus)) +
  geom_boxplot()
p <- p + scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/parasitism_byGenus.pdf", height=5, width=7)


# make graph
parasite.long$Parasite <- factor(parasite.long$Parasite, levels=c("Apicystis", "Crithidia", "Phorid",
                                                                  "CBPV","SBPV", "SBV", "ABPV_KBV_IAPV","BQCV", "DWV_KV_VDV"))

p <- ggplot(parasite.long, aes(x=Parasite, y=Parasitism,
                               fill=Genus)) +
  geom_boxplot()
p <- p + scale_colour_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("figures/parasitism_byGenus.pdf", height=5, width=7)





