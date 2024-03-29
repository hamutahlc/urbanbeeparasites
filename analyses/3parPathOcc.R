setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
source("src/writeResultsTables.R")
source("src/init_bayes.R")

# library(piecewiseSEM)
# library(lme4)
# brms for bayseian regression models 
# install.packages("R2admb")
# install.packages("glmmADMB",
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")

# packageurl <- "http://cran.r-project.org/src/contrib/Archive/rstan/rstan_2.19.3.tar.gz"
# install.packages(packageurl, repos = NULL, type = "source")

library(rstan)
auto_write = TRUE
library(brms)


load("../data/specimens-complete.Rdata")



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
bombusPara <- bombusPara[order(bombusPara$Site),]
apisPara <- apisPara[order(apisPara$Site),]
bombusPath <- bombusPath[order(bombusPath$Site),]
apisPath <- apisPath[order(apisPath$Site),]

bombusPara <- makeDataMultiLevel(bombusPara)
apisPara <- makeDataMultiLevel(apisPara)
bombusPath <- makeDataMultiLevel(bombusPath)
apisPath <- makeDataMultiLevel(apisPath)


# models for site effects on community, make compatible with bayesian analyses
bf.bee.abund <- bf(formula.bee.abund)
bf.bee.div <- bf(formula.bee.div)

# parasites
bf.Phorid.apis <- bf(formulas.par.apis[[1]],  family = "bernoulli")
bf.Crithidia.apis <- bf(formulas.par.apis[[2]],  family="bernoulli")
bf.Apicystis.apis <- bf(formulas.par.apis[[3]],  family="bernoulli")

bf.Phorid.bombus <- bf(formulas.par.bombus[[1]],  family = "bernoulli")
bf.Crithidia.bombus <- bf(formulas.par.bombus[[2]],  family="bernoulli")
bf.Apicystis.bombus <- bf(formulas.par.bombus[[3]],  family="bernoulli")


bf.CBPV.apis <- bf(formulas.path.apis[[1]],  family = "bernoulli")
bf.DWV_KV_VDV.apis <- bf(formulas.path.apis[[2]],  family="bernoulli")
bf.ABPV_KBV_IAPV.apis <- bf(formulas.path.apis[[3]],  family="bernoulli")
bf.BQCV.apis <- bf(formulas.path.apis[[4]],  family = "bernoulli")
bf.SBPV.apis <- bf(formulas.path.apis[[5]],  family="bernoulli")
bf.SBV.apis <- bf(formulas.path.apis[[6]],  family="bernoulli")

bf.CBPV.bombus <- bf(formulas.path.bombus[[1]],  family = "bernoulli")
bf.DWV_KV_VDV.bombus <- bf(formulas.path.bombus[[2]],  family="bernoulli")
bf.ABPV_KBV_IAPV.bombus <- bf(formulas.path.bombus[[3]],  family="bernoulli")
bf.BQCV.bombus <- bf(formulas.path.bombus[[4]],  family = "bernoulli")
bf.SBPV.bombus <- bf(formulas.path.bombus[[5]],  family="bernoulli")
bf.SBV.bombus <- bf(formulas.path.bombus[[6]],  family="bernoulli")

#pathogens


## ****************** HONEY BEE PARASITES *****************

# Apis with Phorid, bayesian
bform.Phorid.apis <- bf.bee.abund + bf.bee.div + bf.Phorid.apis + 
    set_rescor(FALSE)

prior <- c(set_prior("normal(0,1)", class="b"))
fit.Phorid.apis  <- brm(bform.Phorid.apis , apisPara,
                        cores=1,
                        iter = 10^5,
                        chains = 4,
                        control = list(adapt_delta = 0.99))

summary(fit.Phorid.apis)
write.ms.table(fit.Phorid.apis, "Phorid.apis")

save.dir.git <- "~/Dropbox/urbanbeeparasites/data/indivparAndpath"
save(fit.Phorid.apis,
     file=file.path(save.dir.git, "PhoridApis.Rdata"))


## Apis with Crithidia
bform.Crithidia.apis <- bf.bee.abund + bf.bee.div + bf.Crithidia.apis + 
    set_rescor(FALSE)
fit.Crithidia.apis  <- brm(bform.Crithidia.apis , apisPara,
                        cores=1,
                        iter = 10^5,
                        chains = 4,
                        control = list(adapt_delta = 0.99))

summary(fit.Crithidia.apis)
write.ms.table(fit.Crithidia.apis, "Crithidia.apis")

save.dir.git <- "~/Dropbox/urbanbeeparasites/data/indivparAndpath"
save(fit.Crithidia.apis,
     file=file.path(save.dir.git, "CrithidiaApis.Rdata"))

## Apis with Apicystis
bform.Apicystis.apis <- bf.bee.abund + bf.bee.div + bf.Apicystis.apis + 
    set_rescor(FALSE)
fit.Apicystis.apis  <- brm(bform.Apicystis.apis , apisPara,
                           cores=1,
                           iter = 10^5,
                           chains = 4,
                           control = list(adapt_delta = 0.99))

summary(fit.Apicystis.apis)
write.ms.table(fit.Apicystis.apis, "Apicystis.apis")

save(fit.Apicystis.apis,
     file=file.path(save.dir.git, "Apicystisapis.Rdata"))

## ****************** BUMBLE BEE PARASITES *****************
# Bombus with Phorid, bayesian
bform.Phorid.bombus <- bf.bee.abund + bf.bee.div + bf.Phorid.bombus + 
    set_rescor(FALSE)

fit.Phorid.bombus  <- brm(bform.Phorid.bombus , bombusPara,
                        cores=1,
                        iter = 10^5,
                        chains = 4,
                        control = list(adapt_delta = 0.99))

summary(fit.Phorid.bombus)
write.ms.table(fit.Phorid.bombus, "Phorid.bombus")

save(fit.Phorid.bombus,
     file=file.path(save.dir.git, "Phoridbombus.Rdata"))


## bombus with Crithidia
bform.Crithidia.bombus <- bf.bee.abund + bf.bee.div + bf.Crithidia.bombus + 
    set_rescor(FALSE)

prior <- c(set_prior("normal(0,1)", class="b"))
fit.Crithidia.bombus  <- brm(bform.Crithidia.bombus , bombusPara,
                           cores=1,
                           iter = 10^5,
                           chains = 4,
                           prior = prior,
                           control = list(adapt_delta = 0.99))

summary(fit.Crithidia.bombus)
write.ms.table(fit.Crithidia.bombus, "Crithidia.bombus")

save(fit.Crithidia.bombus,
    file=file.path(save.dir.git, "Crithidiabombus.Rdata"))


## Bombus with Apicystis
bform.Apicystis.bombus <- bf.bee.abund + bf.bee.div + bf.Apicystis.bombus + 
    set_rescor(FALSE)

prior <- c(set_prior("normal(0,1)", class="b"))
fit.Apicystis.bombus  <- brm(bform.Apicystis.bombus , bombusPara,
                           cores=1,
                           iter = 10^5,
                           chains = 4,
                           prior = prior,
                           control = list(adapt_delta = 0.99))

summary(fit.Apicystis.bombus)
write.ms.table(fit.Apicystis.bombus, "Apicystis.bombus")

save(fit.Apicystis.bombus,
     file=file.path(save.dir.git, "Apicystisbombus.Rdata"))



load("../data/indivparAndpath/Phoridbombus.Rdata")
load("../data/indivparAndpath/Crithidiabombus.Rdata")
load("../data/indivparAndpath/Apicystisbombus.Rdata")
load("../data/indivparAndpath/PhoridApis.Rdata")
load("../data/indivparAndpath/CrithidiaApis.Rdata")
load("../data/indivparAndpath/Apicystisapis.Rdata")
     
 save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
 save(fit.Phorid.apis, fit.Crithidia.apis, fit.Apicystis.apis, 
      fit.Phorid.bombus, fit.Crithidia.bombus, fit.Apicystis.bombus,
      file=file.path(save.dir.git, "IndivParBayes.Rdata"))




## **********************APIS PATHOGENS********************


prior <- c(set_prior("normal(0,1)", class="b"))

## apis with CBPV
bform.CBPV.apis <- bf.bee.abund + bf.bee.div + bf.CBPV.apis + 
    set_rescor(FALSE)
fit.CBPV.apis  <- brm(bform.CBPV.apis , apisPath,
                        cores=1,
                        iter = 10^5,
                        chains = 4,
                        prior = prior,
                        control = list(adapt_delta = 0.99))

summary(fit.CBPV.apis)
write.ms.table(fit.CBPV.apis, "CBPV.apis")

save.dir.git <- "~/Dropbox/urbanbeeparasites/data/indivparAndpath"
save(fit.CBPV.apis,
     file=file.path(save.dir.git, "CBPVapis.Rdata"))

## apis with DWV_KV_VDV
bform.DWV_KV_VDV.apis <- bf.bee.abund + bf.bee.div + bf.DWV_KV_VDV.apis + 
    set_rescor(FALSE)
fit.DWV_KV_VDV.apis  <- brm(bform.DWV_KV_VDV.apis , apisPath,
                              cores=1,
                              iter = 10^5,
                              chains = 4,
                              prior = prior,
                              control = list(adapt_delta = 0.99))

summary(fit.DWV_KV_VDV.apis)
write.ms.table(fit.DWV_KV_VDV.apis, "DWV_KV_VDV.apis")

save(fit.DWV_KV_VDV.apis,
     file=file.path(save.dir.git, "DWV_KV_VDVapis.Rdata"))


## apis with ABPV_KBV_IAPV
bform.ABPV_KBV_IAPV.apis <- bf.bee.abund + bf.bee.div + bf.ABPV_KBV_IAPV.apis + 
    set_rescor(FALSE)
fit.ABPV_KBV_IAPV.apis  <- brm(bform.ABPV_KBV_IAPV.apis , apisPath,
                                 cores=1,
                                 iter = 10^5,
                                 chains = 4,
                                  prior = prior,
                                 control = list(adapt_delta = 0.99))

summary(fit.ABPV_KBV_IAPV.apis)
write.ms.table(fit.ABPV_KBV_IAPV.apis, "ABPV_KBV_IAPV.apis")

save(fit.ABPV_KBV_IAPV.apis,
     file=file.path(save.dir.git, "ABPV_KBV_IAPVapis.Rdata"))

## apis with BQCV
bform.BQCV.apis <- bf.bee.abund + bf.bee.div + bf.BQCV.apis + 
    set_rescor(FALSE)

prior <- c(set_prior("normal(0,1)", class="b"))
fit.BQCV.apis  <- brm(bform.BQCV.apis , apisPath,
                        cores=1,
                        iter = 10^5,
                        chains = 4,
                        prior = prior,
                        control = list(adapt_delta = 0.99))

summary(fit.BQCV.apis)
write.ms.table(fit.BQCV.apis, "BQCV.apis")

save(fit.BQCV.apis,
     file=file.path(save.dir.git, "BQCVapis.Rdata"))

## apis with SBPV
bform.SBPV.apis <- bf.bee.abund + bf.bee.div + bf.SBPV.apis + 
    set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))
fit.SBPV.apis  <- brm(bform.SBPV.apis , apisPath,
                        cores=1,
                        iter = 10^5,
                        chains = 4, 
                        prior = prior,
                        control = list(adapt_delta = 0.99))

summary(fit.SBPV.apis)
write.ms.table(fit.SBPV.apis, "SBPV.apis")

save(fit.SBPV.apis,
     file=file.path(save.dir.git, "SBPVapis.Rdata"))

## apis with SBV
bform.SBV.apis <- bf.bee.abund + bf.bee.div + bf.SBV.apis + 
    set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))
fit.SBV.apis  <- brm(bform.SBV.apis , apisPath,
                       cores=1,
                       iter = 10^5,
                       chains = 4, 
                       prior = prior,
                       control = list(adapt_delta = 0.99))

summary(fit.SBV.apis)
write.ms.table(fit.SBV.apis, "SBV.apis")

save(fit.SBV.apis,
     file=file.path(save.dir.git, "SBVapis.Rdata"))

## **********************BOMBUS PATHOGENS********************

## Bombus with CBPV
bform.CBPV.bombus <- bf.bee.abund + bf.bee.div + bf.CBPV.bombus + 
    set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))
fit.CBPV.bombus  <- brm(bform.CBPV.bombus , bombusPath,
                        cores=1,
                        iter = 10^5,
                        chains = 4, 
                        prior = prior,
                        control = list(adapt_delta = 0.99))

summary(fit.CBPV.bombus)
write.ms.table(fit.CBPV.bombus, "CBPV.bombus")

save.dir.git <- "~/Dropbox/urbanbeeparasites/data/indivparAndpath"
save(fit.CBPV.bombus,
     file=file.path(save.dir.git, "CBPVbombus.Rdata"))


## Bombus with DWV_KV_VDV
bform.DWV_KV_VDV.bombus <- bf.bee.abund + bf.bee.div + bf.DWV_KV_VDV.bombus + 
    set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))
fit.DWV_KV_VDV.bombus  <- brm(bform.DWV_KV_VDV.bombus , bombusPath,
                        cores=1,
                        iter = 10^5,
                        chains = 4, 
                        prior = prior,
                        control = list(adapt_delta = 0.99))

summary(fit.DWV_KV_VDV.bombus)
write.ms.table(fit.DWV_KV_VDV.bombus, "DWV_KV_VDV.bombus")

save(fit.DWV_KV_VDV.bombus,
     file=file.path(save.dir.git, "DWV_KV_VDVbombus.Rdata"))


## Bombus with ABPV_KBV_IAPV
bform.ABPV_KBV_IAPV.bombus <- bf.bee.abund + bf.bee.div + bf.ABPV_KBV_IAPV.bombus + 
    set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))
fit.ABPV_KBV_IAPV.bombus  <- brm(bform.ABPV_KBV_IAPV.bombus , bombusPath,
                        cores=1,
                        iter = 10^5,
                        chains = 4, 
                        prior = prior,
                        control = list(adapt_delta = 0.99))

summary(fit.ABPV_KBV_IAPV.bombus)
write.ms.table(fit.ABPV_KBV_IAPV.bombus, "ABPV_KBV_IAPV.bombus")

save(fit.ABPV_KBV_IAPV.bombus,
     file=file.path(save.dir.git, "ABPV_KBV_IAPVbombus.Rdata"))


## Bombus with BQCV
bform.BQCV.bombus <- bf.bee.abund + bf.bee.div + bf.BQCV.bombus + 
    set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))
fit.BQCV.bombus  <- brm(bform.BQCV.bombus , bombusPath,
                        cores=1,
                        iter = 10^5,
                        chains = 4, 
                        prior = prior,
                        control = list(adapt_delta = 0.99))

summary(fit.BQCV.bombus)
write.ms.table(fit.BQCV.bombus, "BQCV.bombus")

save(fit.BQCV.bombus,
     file=file.path(save.dir.git, "BQCVbombus.Rdata"))

## Bombus with SBPV
bform.SBPV.bombus <- bf.bee.abund + bf.bee.div + bf.SBPV.bombus + 
    set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))
fit.SBPV.bombus  <- brm(bform.SBPV.bombus , bombusPath,
                       cores=1,
                       iter = 10^5,
                       chains = 4, 
                       prior = prior,
                       control = list(adapt_delta = 0.99))

summary(fit.SBPV.bombus)
write.ms.table(fit.SBPV.bombus, "SBPV.bombus")

save(fit.SBPV.bombus,
     file=file.path(save.dir.git, "SBPVbombus.Rdata"))

## Bombus with SBV
bform.SBV.bombus <- bf.bee.abund + bf.bee.div + bf.SBV.bombus + 
    set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))
fit.SBV.bombus  <- brm(bform.SBV.bombus , bombusPath,
                             cores=1,
                             iter = 10^5,
                             chains = 4, 
                             prior = prior,
                             control = list(adapt_delta = 0.99))

summary(fit.SBV.bombus)
write.ms.table(fit.SBV.bombus, "SBV.bombus")

save(fit.SBV.bombus,
     file=file.path(save.dir.git, "SBVbombus.Rdata"))



load("../data/indivparAndpath/ABPV_KBV_IAPVapis.Rdata")
load("../data/indivparAndpath/ABPV_KBV_IAPVbombus.Rdata")
load("../data/indivparAndpath/BQCVapis.Rdata")
load("../data/indivparAndpath/BQCVbombus.Rdata")
load("../data/indivparAndpath/CBPVapis.Rdata")
load("../data/indivparAndpath/CBPVbombus.Rdata")
load("../data/indivparAndpath/DWV_KV_VDVapis.Rdata")
load("../data/indivparAndpath/DWV_KV_VDVbombus.Rdata")
load("../data/indivparAndpath/SBPVbombus.Rdata")
load("../data/indivparAndpath/SBPVapis.Rdata")
load("../data/indivparAndpath/SBVbombus.Rdata")
load("../data/indivparAndpath/SBVapis.Rdata")



save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(fit.CBPV.apis, fit.SBV.apis, fit.SBPV.apis, fit.ABPV_KBV_IAPV.apis, fit.BQCV.apis, fit.DWV_KV_VDV.apis, 
     fit.CBPV.bombus, fit.SBV.bombus, fit.SBPV.bombus, fit.ABPV_KBV_IAPV.bombus, fit.BQCV.bombus, fit.DWV_KV_VDV.bombus,
     file=file.path(save.dir.git, "IndivPathBayes.Rdata"))
















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


