## *************************************************************************************
## REV 1.5 make formulas for path analyses, with honey bee/bumble bee abund instead of bee richness
## *************************************************************************************
## setwd("/Volumes/Mac 2/Dropbox/urbanbeeparasites/analyses")

setwd("~/Dropbox/urbanbeeparasites/analyses")

rm(list=ls())
library(piecewiseSEM)
library(lme4)

load("../data/specimens-complete.Rdata")

##*********** HONEY BEE ABUNDANCE************
## formula for site effects on the bee community
formula.bee <- formula(ApisAbund~ natural1000m +
                         AbundWoodyFlowers +
                         AbundAnnualFlowers)

ys <- c( "ParasiteRichness", 
        "PathogenRichness")

xvar.par.path <- c("ApisAbund",
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

calcMods <- function(this.formula,
                     formula.bee,
                     col.trials,
                     dats,
                     site.char){
  colnames(dats)[colnames(dats) == col.trials] <- "trials"
  mod = psem(
    BeeDensity = do.call(lm,
                         list(formula=formula.bee,
                              data=site.char)),
    ParPath = do.call(glm,
                      list(formula=this.formula,
                           data = dats,
                           weights=dats$trials,
                           family="binomial"))
  )
  print(summary(mod))
  return(mod)
}

## *************************************************************
#dont need this one
#bombus <- par.path[par.path$Genus == "Bombus",]
apis <- par.path[par.path$Genus == "Apis",]

## honey bees
lapply(formulas.par, calcMods,
       formula.bee=formula.bee,
       col.trials= "ScreenedPar",
       dats=apis,
       site.char=site.char)

apis.mods <- lapply(formulas.par.path, calcMods,
                    formula.bee, apis, site.char)

names(apis.mods) <- ys


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
#apis <- par.path[par.path$Genus == "Apis",]

## bumble bees
bombus.mods <- lapply(formulas.par.path, calcMods,
                    formula.bee, apis, site.char)

#names(apis.mods) <- names(bombus.mods) <- ys
 names(bombus.mods) <- ys

print("bombus")
lapply(bombus.mods, summary)
lapply(bombus.mods, rsquared)

