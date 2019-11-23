## setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(nlme)
library(lme4)
library(lmerTest)

load("../data/specimens-complete.Rdata")
bombus <- par.path[par.path$Genus == "Bombus",]
apis <- par.path[par.path$Genus == "Apis",]

## *************************************************************
## parasite infection
## *************************************************************
##  in Bombus
bombus.par.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                         AnnualFlowerDensity,
                     data = site.char),
    AnyParasite = lme(AnyParasite ~ BeeDensity +
                               natural1000m +
                               WoodyFlowerDensity +
                               AnnualFlowerDensity,
                           random = ~ 1 | Site,
                           data = bombus))
summary(bombus.par.infected.mod)
rsquared(bombus.par.infected.mod)

## parasite infected in Apis

apis.par.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                         AnnualFlowerDensity,
                     data = site.char),
    AnyParasite = lme(AnyParasite ~ BeeDensity +
                               natural1000m +
                               WoodyFlowerDensity +
                               AnnualFlowerDensity,
                           random = ~ 1 | Site,
                           data = apis))
summary(apis.par.infected.mod)
rsquared(apis.par.infected.mod)

## *************************************************************
## parasite richness
## *************************************************************
bombus.par.richness.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                         AnnualFlowerDensity,
                     data = site.char),
    ParasiteRichness = lme(ParasiteRichness ~ BeeDensity +
                               natural1000m +
                               WoodyFlowerDensity +
                               AnnualFlowerDensity,
                           random = ~ 1 | Site,
                           data = bombus))
summary(bombus.par.richness.mod)
rsquared(bombus.par.richness.mod)

## parasite richness in Apis
apis <- par.path[par.path$Genus == "Apis",]
apis.par.richness.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                         AnnualFlowerDensity,
                     data = site.char),
    ParasiteRichness = lme(ParasiteRichness ~ BeeDensity +
                               natural1000m +
                               WoodyFlowerDensity +
                               AnnualFlowerDensity,
                           random = ~ 1 | Site,
                           data = apis))
summary(apis.par.richness.mod)
rsquared(apis.par.richness.mod)

## *************************************************************
## pathogen richness
## *************************************************************
## in Bombus
bombus.path.richness.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                         AnnualFlowerDensity,
                     data = site.char),
    PathogenRichness = lme(PathogenRichness ~ BeeDensity +
                               natural1000m +
                               WoodyFlowerDensity +
                               AnnualFlowerDensity,
                           random = ~ 1 | Site,
                           data = bombus[!is.na(bombus$PathogenRichness),]))
summary(bombus.path.richness.mod)
rsquared(bombus.path.richness.mod)

## pathogen richness in Apis
apis.path.richness.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                         AnnualFlowerDensity,
                     data = site.char),
    PathogenRichness = lme(PathogenRichness ~ BeeDensity +
                               natural1000m +
                               WoodyFlowerDensity +
                               AnnualFlowerDensity,
                           random = ~ 1 | Site,
                           data = apis[!is.na(apis$PathogenRichness),]))
summary(apis.path.richness.mod)
rsquared(apis.path.richness.mod)

## *************************************************************
## pathogen infection
## *************************************************************
##  in Bombus
bombus.path.infected.mod = psem(
    BeeDensity = lme(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                         AnnualFlowerDensity,
                     random = ~ 1 | Site,
                     data = site.char),
    AnyPathogen = lme(AnyPathogen ~ BeeDensity +
                               natural1000m +
                               WoodyFlowerDensity +
                               AnnualFlowerDensity,
                           random = ~ 1 | Site,
                           data = bombus[!is.na(bombus$AnyPathogen),]))
summary(bombus.path.infected.mod)
rsquared(bombus.path.infected.mod)

## pathogen infected in Apis
apis.path.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                         AnnualFlowerDensity,
                     data = site.char),
    AnyPathogen = lme(AnyPathogen ~ BeeDensity +
                               natural1000m +
                               WoodyFlowerDensity +
                               AnnualFlowerDensity,
                           random = ~ 1 | Site,
                           data = apis[!is.na(apis$AnyPathogen),]))
summary(apis.path.infected.mod)
rsquared(apis.path.infected.mod)

## *************************************************************
## sanity check using a linear models of Bee density
## *************************************************************
library(car)
bee.density.mod <- lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                          AnnualFlowerDensity  +
                          PercentBareSoil,
                      data=site.char)

vif(bee.density.mod)
summary(bee.density.mod)
