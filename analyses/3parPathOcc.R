## setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)
library(lmerTest)

load("../data/specimens-complete.Rdata")
bombus <- sick.totals[sick.totals$Genus == "Bombus",]
apis <- sick.totals[sick.totals$Genus == "Apis",]

## *************************************************************
## make formulas for path analyses
## *************************************************************

## formula for site effects on the bee community
formula.bee <- formula(BeeRichnessArea~ natural1000m +
                                   WoodyFlowerDensity +
                                   AnnualFlowerDensity)

## formulas for the site effects on parasites and pathogens

parasites <- c("Phorid", "Crithidia", "Apicystis")
pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV",
               "BQCV","SBPV", "SBV")

xvar.par.path <- c("BeeRichnessArea",
                   "natural1000m",
                   "WoodyFlowerDensity",
                   "AnnualFlowerDensity")

formulas.par <-lapply(parasites, function(x) {
    as.formula(paste(x, "~",
                     paste(xvar.par.path, collapse="+")))
})

formulas.path <-lapply(pathogens, function(x) {
    as.formula(paste(x, "~",
                     paste(xvar.par.path, collapse="+")))
})

## *************************************************************

calcMods <- function(this.formula,
                     formula.bee,
                     col.trials,
                     dats,
                     site.char){
    colnames(dats)[colnames(dats) == col.trials] <- "trials"
    mod = psem(
        BeeDensity = lm(formula.bee,
                        data = site.char),
        ParPath = glm(this.formula,
                       data = dats,
                       weights=trials,
                       family="binomial"))
}

## *************************************************************
## bombus
## *************************************************************

bombus.mods.par <- lapply(formulas.par, calcMods,
                          formula.bee=formula.bee,
                          col.trials= "ScreenedPar",
                          dats=bombus,
                          site.char=site.char)

bombus.mods.path <- lapply(formulas.path, calcMods,
                          formula.bee,
                          col.trials= "ScreenedPath",
                          dats=bombus,
                          site.char=site.char)

names(bombus.mods.par) <- parasites
names(bombus.mods.path) <- pathogens

print("bombus parasites")
lapply(bombus.mods.par, summary)
lapply(bombus.mods.par, rsquared)

print("bombus pathogens")
lapply(bombus.mods.path, summary)
lapply(bombus.mods.path, rsquared)

## *************************************************************
## apis
## *************************************************************

apis.mods.par <- lapply(formulas.par, calcMods,
                          formula.bee=formula.bee,
                          col.trials= "ScreenedPar",
                          dats=apis,
                          site.char=site.char)

apis.mods.path <- lapply(formulas.path, calcMods,
                          formula.bee,
                          col.trials= "ScreenedPath",
                          dats=apis,
                          site.char=site.char)

names(apis.mods.par) <- parasites
names(apis.mods.path) <- pathogens

print("apis parasites")
lapply(apis.mods.par, summary)
lapply(apis.mods.par, rsquared)

print("apis pathogens")
lapply(apis.mods.path, summary)
lapply(apis.mods.path, rsquared)


