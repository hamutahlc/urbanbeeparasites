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

## *************************************************************
## make formulas for path analyses
## *************************************************************
## pearsons correlation test to confirm that variables of interst in our models not collinear 
cor.test(site.char$BeeDiversity, site.char$BeeAbund, method=c("pearson"))
cor.test(site.char$BeeRichness, site.char$BeeAbund, method=c("pearson"))


## formula for site effects on the bee community


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
xvar.par.path.apis <- c(all.indiv.vars,
                        "bombus.parpath.rate")
xvar.par.path.bombus <- c(all.indiv.vars,
                          "apis.parpath.rate")

# scale x variables

par.only[, c(all.indiv.vars, "apis.par.rate", "bombus.par.rate" )] <-
  apply(par.only[, c(all.indiv.vars, "apis.par.rate", "bombus.par.rate" )], 2, scale)

path.only[, c(all.indiv.vars, "apis.path.rate", "bombus.path.rate" )] <-
  apply(path.only[, c(all.indiv.vars, "apis.path.rate", "bombus.path.rate" )], 2, scale)

par.and.path[, c(all.indiv.vars, "apis.parpath.rate", "bombus.parpath.rate" )] <-
  apply(par.and.path[, c(all.indiv.vars, "apis.parpath.rate", "bombus.parpath.rate" )], 2, scale)

# define y variables. 

ys1 <- c("ParasiteRichness", "AnyParasite")
ys2 <- c("PathogenRichness", "AnyPathogen")
ys3 <- c("ParPathRichness", "AnyParPath")

# make formulas


#apis parasites
formulas.par.apis <-lapply(ys1, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.apis, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})

#apis pathogen
formulas.path.apis <-lapply(ys2, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.path.apis, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})

#bombus parasites
formulas.par.bombus <-lapply(ys1, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.bombus, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})

#bombus pathogen
formulas.path.bombus <-lapply(ys2, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.path.bombus, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})

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

## *************************************************************
## psem 
## *************************************************************

# split apis and bombus data
bombusPath <- path.only[path.only$Genus == "Bombus",]
apisPath <- path.only[path.only$Genus == "Apis",]

bombusPara <- par.only[par.only$Genus == "Bombus",]
apisPara <- par.only[par.only$Genus == "Apis",]

bombusParaPath <- par.and.path[par.and.path$Genus == "Bombus",]
apisParaPath <- par.and.path[par.and.path$Genus == "Apis",]

### ***************8working with multi-level data****************

# we have individual-level specimen data, site.char is site-level summary data
# we are trying to put everything in one dataset since it needs to be together for the psem
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

bombusParaPath <- bombusParaPath[order(bombusParaPath$Site),]
apisParaPath <- apisParaPath[order(apisParaPath$Site),]
bombusPara <- bombusPara[order(bombusPara$Site),]
apisPara <- apisPara[order(apisPara$Site),]
bombusPath <- bombusPath[order(bombusPath$Site),]
apisPath <- apisPath[order(apisPath$Site),]

bombusParaPath<- makeDataMultiLevel(bombusParaPath)
apisParaPath<- makeDataMultiLevel(apisParaPath)
bombusPara<- makeDataMultiLevel(bombusPara)
apisPara<- makeDataMultiLevel(apisPara)
bombusPath<- makeDataMultiLevel(bombusPath)
apisPath<- makeDataMultiLevel(apisPath)


# models for site effects on community, make compatible with bayesian analyses
bf.bee.abund <- bf(formula.bee.abund)
bf.bee.div <- bf(formula.bee.div)


bf.parRich.apis <- bf(formulas.par.apis[[1]],  family = binomial)
bf.parAny.apis <- bf(formulas.par.apis[[2]],  family="bernoulli")

bf.pathRich.apis <- bf(formulas.path.apis[[1]],  family = binomial)
bf.pathAny.apis <- bf(formulas.path.apis[[2]],  family="bernoulli")

bf.parRich.bombus <- bf(formulas.par.bombus[[1]],  family = binomial)
bf.parAny.bombus <- bf(formulas.par.bombus[[2]],  family="bernoulli")

bf.pathRich.bombus <- bf(formulas.path.bombus[[1]],  family = binomial)
bf.pathAny.bombus <- bf(formulas.path.bombus[[2]],  family="bernoulli")

# for par and path models im skipping "any" since they mostly all have at least 1 
bf.parpathRich.apis <- bf(formulas.parpath.apis[[1]],  family = binomial)
bf.parpathRich.bombus <- bf(formulas.parpath.bombus[[1]],  family = binomial)


# FINAL MODELS for manuscript R1: combining parasites and pathogens together....
# apis with par-path richness, bayesian

bform.parpathRich.apis <- bf.bee.abund + bf.bee.div + bf.parpathRich.apis + 
  set_rescor(FALSE)

prior <- c(set_prior("normal(0,1)", class="b"))
fit.parpathRich.apis <- brm(bform.parpathRich.apis, apisParaPath,
                        cores=1,
                        iter = 10^5,
                        chains = 4,
                        prior=prior,
                        control = list(adapt_delta = 0.99))

summary(fit.parpathRich.apis)

# diagnostic plot is the trace plot, which is a time series plot of the Markov chains. 
# That is, a trace plot shows the evolution of parameter vector over the iterations of one or many Markov chains.
# confirm that all chains explore same region of parameter values
mcmc_trace(fit.parpathRich.apis)
ggsave("figures/bayesMods/ParPathRichnessApis.pdf",
       height=11, width=8.5)


write.ms.table(fit.parpathRich.apis, "parpathRichnessApis")
# look at rhat. should be around 1. indicates convergence


# bombus with par-path richness, bayesian
bform.parpathRich.bombus <- bf.bee.abund + bf.bee.div + bf.parpathRich.bombus + 
  set_rescor(FALSE)
prior <- c(set_prior("normal(0,1)", class="b"))

fit.parpathRich.bombus <- brm(bform.parpathRich.bombus, bombusParaPath,
                            cores=1,
                            iter = 10^5,
                            chains = 4,
                            prior = prior,
                            control = list(adapt_delta = 0.99))

summary(fit.parpathRich.bombus)

mcmc_trace(fit.parpathRich.bombus)
ggsave("figures/bayesMods/ParPathRichnessBombus.pdf",
       height=11, width=8.5)

write.ms.table(fit.parpathRich.bombus, "parpathRichnessBombus")


## save model results as .Rdata to call in later
save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(fit.parpathRich.bombus, fit.parpathRich.apis,
     file=file.path(save.dir.git, "CommunityHealthresults.Rdata"))


# ###########################################################################

# ## PSEMS HONEY BEES: Parasite rich, Any parasite, Pathogen Richen, Any Pathogen
# 
# # apis with parasite richness, bayesian
# bform.parRich.apis <- bf.bee.abund + bf.bee.div + bf.parRich.apis + 
#   set_rescor(FALSE)
# 
# fit.parRich.apis <- brm(bform.parRich.apis, apisPara,
#            cores=1,
#            iter = 10^4,
#            chains = 3,
#            control = list(adapt_delta = 0.99))
#            
# summary(fit.parRich.apis)
# 
# write.ms.table(fit.parRich.apis, "parRichnessApis")
# # look at rhat. should be around 1. indicates convergence
# # then look at confidence interval. if centered right at 0, not sig
# # bigger gardens marginally -> bee diversity, nat habitat -> bee abund
# # apis parasite richenss impacted by: nothing
# # ran it
# 
# 
# #apis with any parasite, bayesian
# bform.parAny.apis <- bf.bee.abund + bf.bee.div + bf.parAny.apis + 
#   set_rescor(FALSE)
# 
# fit.parAny.apis <- brm(bform.parAny.apis, apisPara,
#                         cores=1,
#                         iter = 10^4,
#                         chains = 3,
#                         control = list(adapt_delta = 0.99))
# 
# 
# summary(fit.parAny.apis)
# write.ms.table(fit.parAny.apis, "parAnyApis")
# # ran it
# 
# 
# #apis with pathogen richness, bayesian
# bform.pathRich.apis <- bf.bee.abund + bf.bee.div + bf.pathRich.apis + 
#   set_rescor(FALSE)
# 
# fit.pathRich.apis <- brm(bform.pathRich.apis, apisPath,
#                         cores=1,
#                         iter = 10^4,
#                         chains = 3,
#                         control = list(adapt_delta = 0.99))
# 
# 
# summary(fit.pathRich.apis)
# write.ms.table(fit.pathRich.apis, "pathRichnessApis")
# # ran it
# 
# 
# #apis with any pathogen, bayesian
# bform.pathAny.apis <- bf.bee.abund + bf.bee.div + bf.pathAny.apis + 
#   set_rescor(FALSE)
# 
# fit.pathAny.apis <- brm(bform.pathAny.apis, apisPath,
#                        cores=1,
#                        iter = 10^4,
#                        chains = 3,
#                        control = list(adapt_delta = 0.99))
# 
# 
# summary(fit.pathAny.apis)
# write.ms.table(fit.pathAny.apis, "pathAnyApis")
# # ran it
# 
# 
# ## BOMBUS: Parasite rich, Any parasite, Pathogen Richen, Any Pathogen
# 
# #bombus with parasite richness, bayesian
# bform.parRich.bombus <- bf.bee.abund + bf.bee.div + bf.parRich.bombus + 
#   set_rescor(FALSE)
# 
# prior <- c(set_prior("normal(0,1)", class="b"))
# fit.parRich.bombus <- brm(bform.parRich.bombus, bombusPara,
#                         cores=1,
#                         iter = 10^5,
#                         chains = 4,
#                         prior=prior,
#                         control = list(adapt_delta = 0.99))
# 
# summary(fit.parRich.bombus)
# 
# write.ms.table(fit.parRich.bombus, "parRichnessBombus")
# # ran it
# 
# #bombus with any parasite, bayesian
# bform.parAny.bombus <- bf.bee.abund + bf.bee.div + bf.parAny.bombus + 
#   set_rescor(FALSE)
# 
# prior <- c(set_prior("normal(0,1)", class="b"))
# fit.parAny.bombus <- brm(bform.parAny.bombus, bombusPara,
#                        cores=1,
#                        iter = 10^5,
#                        chains = 4,
#                        prior=prior,
#                        control = list(adapt_delta = 0.99))
# 
# 
# summary(fit.parAny.bombus)
# write.ms.table(fit.parAny.bombus, "parAnyBombus4")
# # ran it
# 
# 
# #bombus with pathogen richness, bayesian
# bform.pathRich.bombus <- bf.bee.abund + bf.bee.div + bf.pathRich.bombus + 
#   set_rescor(FALSE)
# 
# prior <- c(set_prior("normal(0,1)", class="b"))
# fit.pathRich.bombus <- brm(bform.pathRich.bombus, bombusPath,
#                          cores=1,
#                          iter = 10^5,
#                          chains = 4,
#                          prior=prior,
#                          control = list(adapt_delta = 0.99))
# 
# 
# summary(fit.pathRich.bombus)
# # look at rhat. should be around 1. indicates convergence
# write.ms.table(fit.pathRich.bombus, "pathRichnessBombus")
# # ran it
# 
# 
# #bombus with any pathogen, bayesian
# bform.pathAny.bombus <- bf.bee.abund + bf.bee.div + bf.pathAny.bombus + 
#   set_rescor(FALSE)
# 
# prior <- c(set_prior("normal(0,1)", class="b"))
# fit.pathAny.bombus <- brm(bform.pathAny.bombus, bombusPath,
#                         cores=1,
#                         iter = 10^5,
#                         chains = 4,
#                         prior=prior,
#                         control = list(adapt_delta = 0.99))
# 
# 
# summary(fit.pathAny.bombus)
# write.ms.table(fit.pathAny.bombus, "pathAnyBombuschain")
# # ran it


#
# 
# ## bumble parasite richness 
# bombus.mods.Parrichness  = psem(
#   BeeDensity = do.call(lm,
#                        list(formula=formula.bee,
#                             data=site.char)),
#   ParPath = do.call(glmer,
#                     list(formula=formulas.par.path1[[1]],
#                          data = bombusPara,
#                          weights=bombusPara$PossibleParasite,
#                          family="binomial"))
# )
# 
# summary(bombus.mods.Parrichness)
# 
# site.mod <- do.call(lm,
#         list(formula=formula.bee,
#              data=site.char))
# 
# 
# bombus.mods.Parrichness <- summary(glmer(formula=formulas.par.path1[[1]], 
#                                          weights=bombusPara$PossibleParasite, 
#                                          family="binomial",
#                                          data = bombusPara))
# 
# 
# 
# 
# ## bumble pathogen richness model
# bombus.mods.Pathrichness <- summary(glmer(formula=formulas.par.path1[[2]], 
#                                          weights=bombusPath$PossiblePathogen, 
#                                          family="binomial",
#                                          data = bombusPath))
# 
# 
# 
# 
# ## apis parasite richness 
# apis.mods.Parrichness  = psem(
#   BeeDensity = do.call(lm,
#                        list(formula=formula.bee,
#                             data=site.char)),
#   ParPath = do.call(glmer,
#                     list(formula=formulas.par.path1[[1]],
#                          data = apisPara,
#                          weights=apisPara$PossibleParasite,
#                          family="binomial"))
# )
# 
# site.mod <- do.call(lm,
#                     list(formula=formula.bee,
#                          data=site.char))
# 
# 
# apis.mods.Parrichness <- summary(glmer(formula=formulas.par.path1[[1]], 
#                                          weights=apisPara$PossibleParasite, 
#                                          family="binomial",
#                                          data = apisPara))
# 
# 
# 
# 
# ## apis pathogen richness model
# bombus.mods.Pathrichness <- summary(glmer(formula=formulas.par.path1[[2]], 
#                                           weights=apisPath$PossiblePathogen, 
#                                           family="binomial",
#                                           data = apisPath))
# 
# 
# 
# 
# 
# ###not finished below
# 
# ## honey bees
# apis.mods <- lapply(formulas.par.path, calcMods,
#                     formula.bee, apis, site.char)
# 
# names(apis.mods) <- names(bombus.mods) <- ys
# 
# print("bombus")
# lapply(bombus.mods, summary)
# lapply(bombus.mods, rsquared)
# 
# print("apis")
# lapply(apis.mods, summary)
# lapply(apis.mods, rsquared)
# 
# ## *************************************************************
# ## sanity check using a linear models of Bee density
# ## *************************************************************
# 
library(car)
#bee abundance
bee.abund.mod <- glm(BeeAbund ~ natural1000m +
                         AbundWoodyFlowers +
                         AbundAnnualFlowers + Size,
                       family="poisson",
                       data=site.char)

AIC(bee.abund.mod)
vif(bee.abund.mod)
summary(bee.abund.mod)
plot(density(bee.abund.mod$resid))


#bee div
bee.div.mod <- glm(BeeDiversity ~ natural1000m +
                         AbundWoodyFlowers +
                         AbundAnnualFlowers + Size,
                       family="poisson",
                       data=site.char)

AIC(bee.div.mod)
vif(bee.div.mod)
summary(bee.div.mod)
plot(density(bee.div.mod$resid))


#***** notes
#*
#*#models for parasites and pathogens, make compatible with baysian analyses
# 1 = richness parasite or pathogen, 2 = any parasite or pathogen (rate)

# family is beta binomial. wiht binomial, there is one probability. estimating prob of sucess. since prob varies between parasites
# so use beta binmoial sine probability of sucess varies between trials

# beta_binomial2 <- custom_family(
#   "beta_binomial2", dpars = c("mu", "phi"),
#   links = c("logit", "log"), lb = c(NA, 0),
#   type = "int", vars = "vint1[n]"
# )
# 
# stan_funs <- "
#   real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
# return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
# }
# int beta_binomial2_rng(real mu, real phi, int T) {
# return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
# }
# "
# stanvars <- stanvar(scode = stan_funs, block = "functions")

# bf.parRich.apis <- bf(formulas.par.apis1[[1]],  family = beta_binomial2, stanvars=stanvars)

# i cant get this work, but want it here for future reference
