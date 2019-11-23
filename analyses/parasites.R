## setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(nlme)
library(lme4)
library(lmerTest)

load("../data/specimens-complete.Rdata")

bombus <- sick.totals[sick.totals$Genus == "Bombus",]
apis <- sick.totals[sick.totals$Genus == "Apis",]


## *************************************************************
## Phorid
## *************************************************************
##  in Bombus
bombus.Phorid.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    Phorid = glm(Phorid ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPar,
                 family="binomial",
                 data = bombus))

summary(bombus.Phorid.infected.mod)
rsquared(bombus.Phorid.infected.mod)

## in Apis
apis.Phorid.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    Phorid = glm(Phorid ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPar,
                 family="binomial",
                 data = apis))

summary(apis.Phorid.infected.mod)
rsquared(apis.Phorid.infected.mod)


## *************************************************************
## Crithidia
## *************************************************************
##  in Bombus
bombus.Crithidia.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    Crithidia = glm(Crithidia ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPar,
                 family="binomial",
                 data = bombus))

summary(bombus.Crithidia.infected.mod)
rsquared(bombus.Crithidia.infected.mod)

## in Apis
apis.Crithidia.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    Crithidia = glm(Crithidia ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPar,
                 family="binomial",
                 data = apis))

summary(apis.Crithidia.infected.mod)
rsquared(apis.Crithidia.infected.mod)



## *************************************************************
## Apicystis
## *************************************************************
##  in Bombus
bombus.Apicystis.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    Apicystis = glm(Apicystis ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPar,
                 family="binomial",
                 data = bombus))

summary(bombus.Apicystis.infected.mod)
rsquared(bombus.Apicystis.infected.mod)

## in Apis
apis.Apicystis.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    Apicystis = glm(Apicystis ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPar,
                 family="binomial",
                 data = apis))

summary(apis.Apicystis.infected.mod)
rsquared(apis.Apicystis.infected.mod)




## *************************************************************
## CBPV
## *************************************************************
##  in Bombus
bombus.CBPV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    CBPV = glm(CBPV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = bombus))

summary(bombus.CBPV.infected.mod)
rsquared(bombus.CBPV.infected.mod)

## in Apis
apis.CBPV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    CBPV = glm(CBPV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = apis))

summary(apis.CBPV.infected.mod)
rsquared(apis.CBPV.infected.mod)


## *************************************************************
## DWV_KV_VDV
## *************************************************************
##  in Bombus
bombus.DWV_KV_VDV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    DWV_KV_VDV = glm(DWV_KV_VDV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = bombus))

summary(bombus.DWV_KV_VDV.infected.mod)
rsquared(bombus.DWV_KV_VDV.infected.mod)

## in Apis
apis.DWV_KV_VDV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    DWV_KV_VDV = glm(DWV_KV_VDV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = apis))

summary(apis.DWV_KV_VDV.infected.mod)
rsquared(apis.DWV_KV_VDV.infected.mod)


## *************************************************************
## ABPV_KBV_IAPV
## *************************************************************
##  in Bombus
bombus.ABPV_KBV_IAPV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    ABPV_KBV_IAPV = glm(ABPV_KBV_IAPV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = bombus))

summary(bombus.ABPV_KBV_IAPV.infected.mod)
rsquared(bombus.ABPV_KBV_IAPV.infected.mod)

## in Apis
apis.ABPV_KBV_IAPV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    ABPV_KBV_IAPV = glm(ABPV_KBV_IAPV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = apis))

summary(apis.ABPV_KBV_IAPV.infected.mod)
rsquared(apis.ABPV_KBV_IAPV.infected.mod)



## *************************************************************
## BQCV
## *************************************************************
##  in Bombus
bombus.BQCV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    BQCV = glm(BQCV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = bombus))

summary(bombus.BQCV.infected.mod)
rsquared(bombus.BQCV.infected.mod)

## in Apis
apis.BQCV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    BQCV = glm(BQCV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = apis))

summary(apis.BQCV.infected.mod)
rsquared(apis.BQCV.infected.mod)


## *************************************************************
## SBPV
## *************************************************************
##  in Bombus
bombus.SBPV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    SBPV = glm(SBPV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = bombus))

summary(bombus.SBPV.infected.mod)
rsquared(bombus.SBPV.infected.mod)

## in Apis
apis.SBPV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    SBPV = glm(SBPV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = apis))

summary(apis.SBPV.infected.mod)
rsquared(apis.SBPV.infected.mod)




## *************************************************************
## SBV
## *************************************************************
##  in Bombus
bombus.SBV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    SBV = glm(SBV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = bombus))

summary(bombus.SBV.infected.mod)
rsquared(bombus.SBV.infected.mod)

## in Apis
apis.SBV.infected.mod = psem(
    BeeDensity = lm(BeeDensity ~ natural1000m + WoodyFlowerDensity +
                        AnnualFlowerDensity,
                    data = site.char),
    SBV = glm(SBV ~ BeeDensity +
                     natural1000m +
                     WoodyFlowerDensity +
                     AnnualFlowerDensity,
                 weights=ScreenedPath,
                 family="binomial",
                 data = apis))

summary(apis.SBV.infected.mod)
rsquared(apis.SBV.infected.mod)
