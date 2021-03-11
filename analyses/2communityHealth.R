## setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)

load("../data/specimens-complete.Rdata")

## *************************************************************
## make formulas for path analyses
## *************************************************************

## formula for site effects on the bee community
formula.bee <- formula(BeeRichness~ natural1000m +
                                   WoodyFlowerDensity +
                                   AnnualFlowerDensity)

## formulas for the site effects on parasites and pathogens

ys <- c("AnyParasite", "ParasiteRichness", "AnyPathogen",
        "PathogenRichness")
xvar.par.path <- c("BeeRichness",
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

## bumbles
bombus.mods <- lapply(formulas.par.path, calcMods,
                      formula.bee, bombus, site.char)

## honey bees
apis.mods <- lapply(formulas.par.path, calcMods,
                    formula.bee, apis, site.char)

names(apis.mods) <- names(bombus.mods) <- ys

print("bombus")
lapply(bombus.mods, summary)
lapply(bombus.mods, rsquared)

print("apis")
lapply(apis.mods, summary)
lapply(apis.mods, rsquared)

## *************************************************************
## sanity check using a linear models of Bee density
## *************************************************************

library(car)
#this says density, but we arent using density anymore
bee.density.mod <- lm(BeeAbund ~ natural1000m +
                         AbundWoodyFlowers +
                         AbundAnnualFlowers + Size,
                       data=site.char)

AIC(bee.density.mod)
vif(bee.density.mod)
summary(bee.density.mod)
plot(density(bee.density.mod$resid))
