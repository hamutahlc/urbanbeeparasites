## setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)
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

# NEXT: need to a variable of bombus parasite richness on the apis data, and vice versa.
# call it other.bee.rate in each dataset baesd on any parasiterate. should be like 0, .5, etc. 
#make apis.par.rate, bombus.par.rate, apis.path.rate, bombus.path.rate

# take path.only$anypathogen and average across sites for just honey bees, and for just bumble bees
# take par.only$any parasite and average across sites for just honey bees, and for just bumble bees

par.rates <- par.only %>%
  group_by(Site, Genus) %>%
  summarise(par.rates = mean(AnyParasite))

path.rates <- path.only %>%
  group_by(Site, Genus) %>%
  summarise(path.rates = mean(AnyPathogen))

apis.par.rate <- par.rates[par.rates$Genus == "Apis",]
bombus.par.rate <- par.rates[par.rates$Genus == "Bombus",]
apis.path.rate <- path.rates[path.rates$Genus == "Apis",]
bombus.path.rate <- path.rates[path.rates$Genus == "Bombus",]

colnames(apis.par.rate)[colnames(apis.par.rate) == "par.rates"]  <-
  "apis.par.rate"
colnames(bombus.par.rate)[colnames(bombus.par.rate) == "par.rates"]  <-
  "bombus.par.rate"
colnames(apis.path.rate)[colnames(apis.path.rate) == "path.rates"]  <-
  "apis.path.rate"
colnames(bombus.path.rate)[colnames(bombus.path.rate) == "path.rates"]  <-
  "bombus.path.rate"

apis.par.rate$Genus <- NULL
bombus.par.rate$Genus <- NULL
apis.path.rate$Genus <- NULL
bombus.path.rate$Genus <- NULL

site.char <- merge(site.char, apis.par.rate,  all.x=TRUE)
site.char <- merge(site.char, bombus.par.rate,  all.x=TRUE)
site.char <- merge(site.char, apis.path.rate,  all.x=TRUE)
site.char <- merge(site.char, bombus.path.rate,  all.x=TRUE)

site.char$apis.par.rate[is.na(site.char$apis.par.rate)] <- 0
site.char$bombus.par.rate[is.na(site.char$bombus.par.rate)] <- 0
site.char$apis.path.rate[is.na(site.char$apis.path.rate)] <- 0
site.char$bombus.path.rate[is.na(site.char$bombus.path.rate)] <- 0


## formula for site effects on the bee community
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

ys1 <- c("ParasiteRichness",
        "PathogenRichness")

ys2 <- c("AnyParasite", "AnyPathogen")

## do we need size as a covariate here as well? add honey bee parsitism/bumble bee parsitism
xvar.par.path <- c("BeeDiversity", "BeeAbund",
                   "natural1000m",
                   "AbundWoodyFlowers",
                   "AbundAnnualFlowers",
                   "Size",
                   "other.bee.rate")

par.path[, xvar.par.path] <- apply(par.path[, xvar.par.path], 2, scale)
par.only[, xvar.par.path] <- apply(par.only[, xvar.par.path], 2, scale)
site.char[, xvar.par.path] <- apply(site.char[, xvar.par.path], 2, scale)

formulas.par.path1 <-lapply(ys1, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.par.path, collapse="+"),
                           "(1|Site)",
                           sep="+")))
})

formulas.par.path2 <-lapply(ys2, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.par.path, collapse="+"),
                         "(1|Site)",
                         sep="+")))
})




## *************************************************************

# NEXT: need to a variable of bombus parasite richness on the apis data, and vice versa.
# call it other.bee.rate in each dataset baesd on any parasiterate. should be like 0, .5, etc. 


# split apis and bombus data
bombusPath <- par.path[par.path$Genus == "Bombus",]
apisPath <- par.path[par.path$Genus == "Apis",]

bombusPara <- par.only[par.only$Genus == "Bombus",]
apisPara <- par.only[par.only$Genus == "Apis",]


# pass in multiple datasets into psem with weights, which drops duplicate data
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

bombusPara <- makeDataMultiLevel(bombusPara)
apisPara <- makeDataMultiLevel(apisPara)
bombusPath <- makeDataMultiLevel(bombusPath)
apisPath <- makeDataMultiLevel(apisPath)

#this is bayesian version of psem
bf.bee.abund <- bf(formula.bee.abund)
bf.bee.div <- bf(formula.bee.div)

bf.par <- bf(formulas.par.path2[[1]],  family="bernoulli")

bf.path <- bf(formulas.par.path2[[2]], family="bernoulli")


bf.par.rich <- bf(formulas.par.path1[[1]],  family="binomial")

bf.path.rich <- bf(formulas.par.path1[[2]], family="binomial")

#bumble bee with any parasites, bayesian
bform.para <- bf.bee.abund + bf.bee.div + bf.par + 
  set_rescor(FALSE)


fit.para <- brm(bform.para, bombusPara,
           cores=1,
           iter = 10^4,
           chains = 3,
           control = list(adapt_delta = 0.99))
           

summary(fit.para)
# look at rhat. should be around 1. indicates convergence
# then look at confidence interval. if centered right at 0, not sig
# bigger gardens -> bee diversity, nat habitat -> bee abund
# no impact to any parasite in bumbles

## bumble parasite richness, bayseian
bform.para.rich <- bf.bee.abund + bf.bee.div + bf.par.rich + 
  set_rescor(FALSE)

fit.para.rich <- brm(bform.para.rich, bombusPara,
                cores=1,
                iter = 10^4,
                chains = 3,
                control = list(adapt_delta = 0.99))


summary(fit.para.rich)
















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



