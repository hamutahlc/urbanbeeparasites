## *************************************************************************************
## REV 1.5 make formulas for path analyses, with honey bee/bumble bee abund instead of bee richness
## *************************************************************************************
setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
source("src/writeResultsTables.R")
source("src/init_bayes.R")
library(piecewiseSEM)
library(lme4)

# brms for bayseian regression models 
install.packages("R2admb")
install.packages("glmmADMB",
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")
library(brms)

load("../data/specimens-complete.Rdata")


##***********Prep data for bayesian analyses************

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


par.and.path <- par.and.path [order(par.and.path $Site),]
par.and.path <- makeDataMultiLevel(par.and.path)




## formula for site effects on the bee community
formula.apis.abund <- formula(ApisAbund| weights(Weights)~ natural1000m +
                               AbundWoodyFlowers +
                               AbundAnnualFlowers  + Size)

formula.bombus.abund <- formula(BombusAbund| weights(Weights)~ natural1000m +
                             AbundWoodyFlowers +
                             AbundAnnualFlowers  + Size)


## formulas for the site effects on parasites and pathogens



xvar.par.path.apis <- c("ApisAbund",
                        "natural1000m",
                        "AbundWoodyFlowers",
                        "AbundAnnualFlowers",
                        "Size",
                        "bombus.parpath.rate")

xvar.par.path.bombus <- c("BombusAbund",
                          "natural1000m",
                          "AbundWoodyFlowers",
                          "AbundAnnualFlowers",
                          "Size",
                          "apis.parpath.rate")

# scale x variables

par.and.path[, c("ApisAbund", "BombusAbund",
                 "natural1000m", "AbundWoodyFlowers",
                 "AbundAnnualFlowers","Size",
                 "bombus.parpath.rate", "apis.parpath.rate")] <-
apply(par.and.path[, c("ApisAbund", "BombusAbund",
                       "natural1000m", "AbundWoodyFlowers",
                       "AbundAnnualFlowers","Size",
                       "bombus.parpath.rate", "apis.parpath.rate" )], 2, scale)

# define y variables. 

ys3 <- c("ParPathRichness", "AnyParPath")


# make formulas

# apis parpath
formulas.parpath.apis <-lapply(ys3, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.path.apis, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})

# bombus parpath
formulas.parpath.bombus <-lapply(ys3, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.path.bombus, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})

## **********************psem ***************************************


# split apis and bombus data

bombusParaPath <- par.and.path[par.and.path$Genus == "Bombus",]
apisParaPath <- par.and.path[par.and.path$Genus == "Apis",]

# models for site effects on community, make compatible with bayesian analyses
bf.apis.abund <- bf(formula.apis.abund)
bf.bombus.abund <- bf(formula.bombus.abund)

#models

bf.parpathRich.apis <- bf(formulas.parpath.apis[[1]],  family = binomial)
bf.parpathRich.bombus <- bf(formulas.parpath.bombus[[1]],  family = binomial)

#psem for par and path richness in APIS

bform.parpathRich.apis <- bf.apis.abund  + bf.parpathRich.apis + 
  set_rescor(FALSE)

prior <- c(set_prior("normal(0,1)", class="b"))

fit.parpathRich.apis.HostAbund <- brm(bform.parpathRich.apis, apisParaPath,
                            cores=1,
                            iter = 10^5,
                            chains = 4,
                            prior=prior,
                            control = list(adapt_delta = 0.99))

summary(fit.parpathRich.apis.HostAbund)

write.ms.table(fit.parpathRich.apis.HostAbund,"parpathRich.apis.HostAbund")

save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(fit.parpathRich.apis.HostAbund, 
     file=file.path(save.dir.git, "CommunityHealthHostAbundResults.Rdata"))


#psem for par and path richness in Bombus

bform.parpathRich.bombus <- bf.bombus.abund  + bf.parpathRich.bombus + 
  set_rescor(FALSE)

prior <- c(set_prior("normal(0,1)", class="b"))

fit.parpathRich.bombus.HostAbund <- brm(bform.parpathRich.bombus, bombusParaPath,
                                      cores=1,
                                      iter = 10^5,
                                      chains = 4,
                                      prior=prior,
                                      control = list(adapt_delta = 0.99))

summary(fit.parpathRich.bombus.HostAbund)

write.ms.table(fit.parpathRich.bombus.HostAbund,"parpathRich.bombus.HostAbund")

save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(fit.parpathRich.apis.HostAbund, fit.parpathRich.bombus.HostAbund, 
     file=file.path(save.dir.git, "CommunityHealthHostAbundResults.Rdata"))


##***********Compare host abund models to original models************

load("~/Dropbox/urbanbeeparasites/data/CommunityHealthresults.Rdata")
load("~/Dropbox/urbanbeeparasites/data/CommunityHealthHostAbundResults.Rdata")

waic(fit.parpathRich.bombus.HostAbund, fit.parpathRich.bombus, compare=TRUE,pointwise = FALSE)

waic(fit.parpathRich.apis.HostAbund, fit.parpathRich.apis, compare=TRUE)

loo(fit.parpathRich.bombus.HostAbund)
loo(fit.parpathRich.bombus)
loo(fit.parpathRich.apis.HostAbund)
loo(fit.parpathRich.apis)