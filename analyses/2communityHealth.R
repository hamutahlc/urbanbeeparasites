setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
source("src/writeResultsTables.R")
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

## *************************************************************
## make formulas for path analyses
## *************************************************************
## pearsons correlation test to confirm that variables of interst in our models not collinear 
cor.test(site.char$BeeDiversity, site.char$BeeAbund, method=c("pearson"))
cor.test(site.char$BeeRichness, site.char$BeeAbund, method=c("pearson"))


# par.only and path.only are individual-level specimen data, site.char is site-level summary data
# we are trying to put everythign in one dataset since it needs to be together for the psem
# we pass in multiple datasets into psem with weights, which drops duplicate data for each site
# smash site.char into individual-level parasite data (par.only)
# then again into pathogen (par.only) data (need to do it twice)/ 

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

par.only <- par.only[order(par.only$Site),]
path.only <- path.only[order(path.only$Site),]

path.only <- makeDataMultiLevel(path.only)
par.only <- makeDataMultiLevel(par.only)


## formula for site effects on the bee community
## because we use the weights, we use site characterstic data from one site data point 

formula.bee.rich <- formula(BeeRichness| weights(Weights)~ natural1000m +
                                   AbundWoodyFlowers +
                                   AbundAnnualFlowers  + Size)

formula.bee.abund <- formula(BeeAbund| weights(Weights)~ natural1000m +
                         AbundWoodyFlowers +
                         AbundAnnualFlowers  + Size)

formula.bee.div <- formula(BeeDiversity| weights(Weights)~ natural1000m +
                          AbundWoodyFlowers +
                          AbundAnnualFlowers  + Size)


## formulas for the site effects on parasites and pathogens
## each formula should have the same set of site and bee community variables, "all.indiv.variables"
## in addition,  need to add other.par.rate/other.path.rate

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

# ys1 <- c("ParasiteRichness",
#          "PathogenRichness")
# 
# ys2 <- c("AnyParasite", "AnyPathogen")

ys1 <- c("ParasiteRichness",
         "AnyParasite")

ys2 <- c("PathogenRichness", "AnyPathogen")

# make formulas


# formulas.par.path1 <-lapply(ys1, function(x) {
#     as.formula(paste(x, "~",
#                      paste(paste(xvar.par.path, collapse="+"),
#                            "(1|Site)",
#                            sep="+")))
# })

#apis parasites
formulas.par.apis1 <-lapply(ys1, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.apis, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})

#apis pathogen
formulas.par.apis2 <-lapply(ys2, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.path.apis, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})


#bombus parasites
formulas.par.bombus1 <-lapply(ys1, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.bombus, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})

#bombus pathogen
formulas.par.bombus2 <-lapply(ys2, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.path.bombus, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})



## *************************************************************
## psem 

# split apis and bombus data
bombusPath <- path.only[path.only$Genus == "Bombus",]
apisPath <- path.only[path.only$Genus == "Apis",]

bombusPara <- par.only[par.only$Genus == "Bombus",]
apisPara <- par.only[par.only$Genus == "Apis",]


# models for site effects on community, make compatible with bayesian analyses
bf.bee.abund <- bf(formula.bee.abund)
bf.bee.div <- bf(formula.bee.div)
bf.bee.rich <- bf(formula.bee.rich)


#models for parasites and pathogens, make compatible with baysian analyses
# 1 = richness parasite or pathogen, 2 = any parasite or pathogen (rate)

# family is beta binomial. wiht binomial, there is one probability. estimating prob of sucess. since prob varies between parasites
# so use beta binmoial sine probability of sucess varies between trials

beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
}
int beta_binomial2_rng(real mu, real phi, int T) {
return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
}
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

bf.parRich.apis <- bf(formulas.par.apis1[[1]],  family = beta_binomial2, stanvars=stanvars)
bf.parAny.apis <- bf(formulas.par.apis1[[2]],  family="bernoulli")

bf.pathRich.apis <- bf(formulas.par.apis2[[1]],  family = beta_binomial2)
bf.pathAny.apis <- bf(formulas.par.apis2[[2]],  family="bernoulli")


bf.parRich.bombus <- bf(formulas.par.bombus1[[1]],  family = beta_binomial2)
bf.parAny.bombus <- bf(formulas.par.bombus1[[2]],  family="bernoulli")

bf.pathRich.bombus <- bf(formulas.par.bombus2[[1]],  family = beta_binomial2)
bf.pathAny.bombus <- bf(formulas.par.bombus2[[2]],  family="bernoulli")


## PSEMS HONEY BEES: Parasite rich, Any parasite, Pathogen Richen, Any Pathogen

# 


# apis with parasite richness, bayesian
bform.parRich.apis <- bf.bee.abund + bf.bee.div + bf.parRich.apis + 
  set_rescor(FALSE)

fit.parRich.apis <- brm(bform.parRich.apis, apisPara,
           cores=1,
           iter = 10^4,
           chains = 3,
           control = list(adapt_delta = 0.99))
           

summary(fit.parRich.apis)

write.ms.table(fit.parRich.apis, "parRichnessApis")
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens marginally -> bee diversity, nat habitat -> bee abund
# apis parasite richenss impacted by: nothing


#apis with any parasite, bayesian
bform.parAny.apis <- bf.bee.abund + bf.bee.div + bf.parAny.apis + 
  set_rescor(FALSE)

fit.parAny.apis <- brm(bform.parAny.apis, apisPara,
                        cores=1,
                        iter = 10^4,
                        chains = 3,
                        control = list(adapt_delta = 0.99))


summary(fit.parAny.apis)
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens -> bee diversity, nat habitat -> bee abund
# apis parasite rate impacted by: nothing


## work on this one now
#apis with pathogen richness, bayesian
bform.pathRich.apis <- bf.bee.abund + bf.bee.div + bf.pathRich.apis + 
  set_rescor(FALSE)

fit.pathRich.apis <- brm(bform.pathRich.apis, apisPara,
                        cores=1,
                        iter = 10^4,
                        chains = 3,
                        control = list(adapt_delta = 0.99))


summary(fit.pathRich.apis)
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens -> bee diversity, nat habitat -> bee abund
# apis pathogen richenss impacted by: 


#apis with any pathogen, bayesian
bform.pathAny.apis <- bf.bee.abund + bf.bee.div + bf.pathAny.apis + 
  set_rescor(FALSE)

fit.pathAny.apis <- brm(bform.pathAny.apis, apisPara,
                       cores=1,
                       iter = 10^4,
                       chains = 3,
                       control = list(adapt_delta = 0.99))


summary(fit.pathAny.apis)
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens -> bee diversity, nat habitat -> bee abund
# apis pathogen rate impacted by: 


## BOMBUS: Parasite rich, Any parasite, Pathogen Richen, Any Pathogen

#bombus with parasite richness, bayesian
bform.parRich.bombus <- bf.bee.abund + bf.bee.div + bf.parRich.bombus + 
  set_rescor(FALSE)

fit.parRich.bombus <- brm(bform.parRich.bombus, bombusPara,
                        cores=1,
                        iter = 10^4,
                        chains = 3,
                        control = list(adapt_delta = 0.99))


summary(fit.parRich.bombus)
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens -> bee diversity, nat habitat -> bee abund
# bombus parasite richenss impacted by: 


#bombus with any parasite, bayesian
bform.parAny.bombus <- bf.bee.abund + bf.bee.div + bf.parAny.bombus + 
  set_rescor(FALSE)

fit.parAny.bombus <- brm(bform.parAny.bombus, bombusPara,
                       cores=1,
                       iter = 10^4,
                       chains = 3,
                       control = list(adapt_delta = 0.99))


summary(fit.parAny.bombus)
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens -> bee diversity, nat habitat -> bee abund
# apis parasite rate impacted by: 



#apis with pathogen richness, bayesian
bform.pathRich.bombus <- bf.bee.abund + bf.bee.div + bf.pathRich.bombus + 
  set_rescor(FALSE)

fit.pathRich.bombus <- brm(bform.pathRich.bombus, bombusPara,
                         cores=1,
                         iter = 10^4,
                         chains = 3,
                         control = list(adapt_delta = 0.99))


summary(fit.pathRich.bombus)
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens -> bee diversity, nat habitat -> bee abund
# bombus pathogen richenss impacted by: 


#apis with any pathogen, bayesian
bform.pathAny.bombus <- bf.bee.abund + bf.bee.div + bf.pathAny.bombus + 
  set_rescor(FALSE)

fit.pathAny.bombus <- brm(bform.pathAny.bombus, bombusPara,
                        cores=1,
                        iter = 10^4,
                        chains = 3,
                        control = list(adapt_delta = 0.99))


summary(fit.pathAny.bombus)
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens -> bee diversity, nat habitat -> bee abund
# bombus pathogen rate impacted by: 


###******************* EXPORT results as tables  *******************



# results: 
summary(fit.parRich.apis)
summary(fit.parAny.apis)
summary(fit.pathRich.apis)
summary(fit.pathAny.apis)
summary(fit.parRich.bombus)
summary(fit.parAny.bombus)
summary(fit.pathRich.bombus)
summary(fit.pathAny.bombus)


# to do: 
# show LCP the conceptual figure
# run all these models, output them in a nice csv or smething easy to read
# make conceptual figures of them and a table - which goes into the ms? 
# how do i report statistics to evaluate model fit? d-test shows if we could improve models with exclusion of hypothesized paths
  #or inclusion of non hypothesized paths. 
  
# make scatterplots/effect plots of major findings

# rerun with honey bee abund or with bumble bee abund instead of bee diversity or bee abund. 
    # how do i evaluate which model is better?










## bumble parasite richness 
bombus.mods.Parrichness  = psem(
  BeeDensity = do.call(lm,
                       list(formula=formula.bee,
                            data=site.char)),
  ParPath = do.call(glmer,
                    list(formula=formulas.par.path1[[1]],
                         data = bombusPara,
                         weights=bombusPara$PossibleParasite,
                         family="binomial"))
)

summary(bombus.mods.Parrichness)

site.mod <- do.call(lm,
        list(formula=formula.bee,
             data=site.char))


bombus.mods.Parrichness <- summary(glmer(formula=formulas.par.path1[[1]], 
                                         weights=bombusPara$PossibleParasite, 
                                         family="binomial",
                                         data = bombusPara))




## bumble pathogen richness model
bombus.mods.Pathrichness <- summary(glmer(formula=formulas.par.path1[[2]], 
                                         weights=bombusPath$PossiblePathogen, 
                                         family="binomial",
                                         data = bombusPath))




## apis parasite richness 
apis.mods.Parrichness  = psem(
  BeeDensity = do.call(lm,
                       list(formula=formula.bee,
                            data=site.char)),
  ParPath = do.call(glmer,
                    list(formula=formulas.par.path1[[1]],
                         data = apisPara,
                         weights=apisPara$PossibleParasite,
                         family="binomial"))
)

site.mod <- do.call(lm,
                    list(formula=formula.bee,
                         data=site.char))


apis.mods.Parrichness <- summary(glmer(formula=formulas.par.path1[[1]], 
                                         weights=apisPara$PossibleParasite, 
                                         family="binomial",
                                         data = apisPara))




## apis pathogen richness model
bombus.mods.Pathrichness <- summary(glmer(formula=formulas.par.path1[[2]], 
                                          weights=apisPath$PossiblePathogen, 
                                          family="binomial",
                                          data = apisPath))





###not finished below

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
bee.density.mod <- glm(BeeAbund ~ natural1000m +
                         AbundWoodyFlowers +
                         AbundAnnualFlowers + Size,
                       family="poisson",
                       data=site.char)

AIC(bee.density.mod)
vif(bee.density.mod)
summary(bee.density.mod)
plot(density(bee.density.mod$resid))



