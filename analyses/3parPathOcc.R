setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
source("src/writeResultsTables.R")
source("src/init_bayes.R")
library(piecewiseSEM)
library(lme4)

# brms for bayseian regression models 
# install.packages("R2admb")
# install.packages("glmmADMB",
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")
library(brms)

load("../data/specimens-complete.Rdata")

## make weights. take site.char and pass it to the function that lauren wrote
makeDataMultiLevel <- function(indiv.data){
    site.ids <- unlist(tapply(indiv.data$Site,
                              indiv.data$Site,
                              function(x) 1:length(x)))
    names(site.ids) <- NULL
    indiv.data$SiteIDs <- site.ids
    indiv.data$Weights <- indiv.data$SiteIDs
    indiv.data$Weights[indiv.data$Weights > 1] <- 0
    return(indiv.data)
}

#we want to run the pathogen models from path.only and parasites from par.only cause par.and.path has lots of samples dropped 
par.only <- par.only[order(par.only$Site),]
path.only <- path.only[order(path.only$Site),]


path.only <- makeDataMultiLevel(path.only)
par.only <- makeDataMultiLevel(par.only)



#site effects on bee community
formula.bee.abund <- formula(BeeAbund| weights(Weights)~ natural1000m +
                                 AbundWoodyFlowers +
                                 AbundAnnualFlowers  + Size)

formula.bee.div <- formula(BeeDiversity| weights(Weights)~ natural1000m +
                               AbundWoodyFlowers +
                               AbundAnnualFlowers  + Size)


## formulas for the site effects on parasites and pathogens
all.indiv.vars <- c("BeeDiversity", "BeeAbund",
                    "natural1000m",
                    "AbundWoodyFlowers",
                    "AbundAnnualFlowers",
                    "Size")

xvar.par.apis <- c(all.indiv.vars, 
                   "bombus.par.rate")

xvar.par.bombus <- c(all.indiv.vars,
                     "apis.par.rate")

xvar.path.apis <- c(all.indiv.vars,
                    "bombus.path.rate")

xvar.path.bombus <- c(all.indiv.vars,
                      "apis.path.rate")




# scale x variables

par.only[, c(all.indiv.vars, "apis.par.rate", "bombus.par.rate" )] <-
    apply(par.only[, c(all.indiv.vars, "apis.par.rate", "bombus.par.rate" )], 2, scale)

path.only[, c(all.indiv.vars, "apis.path.rate", "bombus.path.rate" )] <-
    apply(path.only[, c(all.indiv.vars, "apis.path.rate", "bombus.path.rate" )], 2, scale)



# define y variables. 


parasites <- c("Phorid", "Crithidia", "Apicystis")
pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV",
               "BQCV","SBPV", "SBV")

# make formulas

formulas.par.apis <-lapply(parasites, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.par.apis, collapse="+"),
                           "(1|Site)",
                           sep="+")))
})


formulas.par.bombus <-lapply(parasites, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.par.bombus, collapse="+"),
                           "(1|Site)",
                           sep="+")))
})


formulas.path.apis <-lapply(pathogens, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.path.apis, collapse="+"),
                           "(1|Site)",
                           sep="+")))
})


formulas.path.bombus <-lapply(pathogens, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.path.bombus, collapse="+"),
                           "(1|Site)",
                           sep="+")))
})



## psem 

# split apis and bombus data
bombusPath <- path.only[path.only$Genus == "Bombus",]
apisPath <- path.only[path.only$Genus == "Apis",]

bombusPara <- par.only[par.only$Genus == "Bombus",]
apisPara <- par.only[par.only$Genus == "Apis",]


# models for site effects on community, make compatible with bayesian analyses
bf.bee.abund <- bf(formula.bee.abund)
bf.bee.div <- bf(formula.bee.div)


bf.Phorid.apis <- bf(formulas.par.apis[[1]],  family = "bernoulli")
bf.Crithidia.apis <- bf(formulas.par.apis[[2]],  family="bernoulli")
bf.Apicystis.apis <- bf(formulas.par.apis[[2]],  family="bernoulli")

bf.Phorid.bombus <- bf(formulas.par.bombus[[1]],  family = "bernoulli")
bf.Crithidia.bombus <- bf(formulas.par.bombus[[2]],  family="bernoulli")
bf.Apicystis.bombuys <- bf(formulas.par.bombus[[2]],  family="bernoulli")

#next do pathogens









## *************************************************************
## manuscript R.0. below
## *************************************************************



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
formula.bee <- formula(BeeRichness~ natural1000m +
                         AbundWoodyFlowers +
                         AbundAnnualFlowers + Size)

## formulas for the site effects on parasites and pathogens

parasites <- c("Phorid", "Crithidia", "Apicystis")
pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV",
               "BQCV","SBPV", "SBV")

xvar.par.path <- c("BeeRichness",
                   "natural1000m",
                   "AbundWoodyFlowers",
                   "AbundAnnualFlowers")

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


