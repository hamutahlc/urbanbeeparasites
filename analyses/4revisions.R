## *************************************************************************************
## REV 1.5 make formulas for path analyses, with honey bee/bumble bee instead of bee richness
## *************************************************************************************

setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)

load("../data/specimens-complete.Rdata")

##*********** HONEY BEE ABUNDANCE************
## formula for site effects on the bee community
formula.bee <- formula(HBAbund~ natural1000m +
                         WoodyFlowerDensity +
                         AnnualFlowerDensity)

ys <- c("AnyParasite", "ParasiteRichness", "AnyPathogen",
        "PathogenRichness")
xvar.par.path <- c("HBAbund",
                   "natural1000m",
                   "AbundWoodyFlowers",
                   "AbundAnnualFlowers")

formulas.par.path <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.path, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})


## *************************************************************

calcMods <- function(this.formula, formula.bee,
                     dats,
                     site.char){
  mod = psem(
    BeeDensity = lm(formula.bee,
                    data = site.char),
    ParPath = lmer(this.formula,
                   data = dats))
}

## *************************************************************
bombus <- par.path[par.path$Genus == "Bombus",]
apis <- par.path[par.path$Genus == "Apis",]

## honey bees
apis.mods <- lapply(formulas.par.path, calcMods,
                    formula.bee, apis, site.char)

#names(apis.mods) <- names(bombus.mods) <- ys

print("apis")
lapply(apis.mods, summary)
lapply(apis.mods, rsquared)

##*********** BUMBLE BEE ABUNDANCE************

## formula for site effects on the bee community
formula.bee <- formula(BBAbund~ natural1000m +
                         WoodyFlowerDensity +
                         AnnualFlowerDensity)

ys <- c("AnyParasite", "ParasiteRichness", "AnyPathogen",
        "PathogenRichness")
xvar.par.path <- c("BBAbund",
                   "natural1000m",
                   "AbundWoodyFlowers",
                   "AbundAnnualFlowers")

formulas.par.path <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.path, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})


## *************************************************************

calcMods <- function(this.formula, formula.bee,
                     dats,
                     site.char){
  mod = psem(
    BeeDensity = lm(formula.bee,
                    data = site.char),
    ParPath = lmer(this.formula,
                   data = dats))
}

## *************************************************************
bombus <- par.path[par.path$Genus == "Bombus",]
apis <- par.path[par.path$Genus == "Apis",]

## bumble bees
bombus.mods <- lapply(formulas.par.path, calcMods,
                    formula.bee, apis, site.char)

#names(apis.mods) <- names(bombus.mods) <- ys

print("bombus")
lapply(bombus.mods, summary)
lapply(bombus.mods, rsquared)

