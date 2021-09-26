rm(list=ls())
## prepares raw data and creates dataset for analyses
library(tidyverse)
library(vegan)

save.dir <- "~/Dropbox/urbanbeeparasites_saved"
## save.dir <- "/Volumes/Mac\ 2/Dropbox/urbanbeeparasites_saved"

bees <- read.csv(file.path(save.dir,
                           "BeeDiversity/BeeIDs2per2015.csv"),
                 stringsAsFactors=FALSE)
veg <- read.csv(file.path(save.dir,
                          "VegData/VegData2Per2015.csv"),
                stringsAsFactors=FALSE)
apisPP <- read.csv(file.path(save.dir,
                             "Parasitism/ApisParasitism.csv"))
vosPP <- read.csv(file.path(save.dir,
                            "Parasitism/BombusVosParasitism.csv"),
                  stringsAsFactors=FALSE)

## *************************************************************
## standardize columns between apis and bombus data
## two different ID columns
apisPP$IDNum <- NULL
## date column missing
vosPP$Date <- NA

vosPP$Genus  <- "Bombus"
vosPP$Species  <- "vosnesenskii"

## line up thecolumns
vosPP <- vosPP[, colnames(apisPP)]

## create a single parasite and pathogen dataset for both species
path.only <- rbind(apisPP, vosPP)

## fix site names to match between datasets
path.only$Site <- as.character(path.only$Site)
path.only$Site[path.only$Site ==
              "Coyote"]  <- "CoyoteCreek"
path.only$Site[path.only$Site ==
              "Beach"]  <- "BeachFlats"
path.only$Site[path.only$Site ==
              "Beach Flats"]  <- "BeachFlats"
path.only$Site[path.only$Site ==
              "La Colina"]  <- "LaColina"
path.only$Site[path.only$Site ==
              "Laguna "]  <- "LagunaSeca"
path.only$Site[path.only$Site ==
              "Laguna"]  <- "LagunaSeca"
path.only$Site[path.only$Site ==
              "Prusch"]  <- "PruschPark"

## drop all the controls
path.only <- path.only[!path.only$Site == "Control",]

## i think judt want to make a separate parasite datasheet to use with the parasite data, 
## since the way the screenign occured is that a subset of bees sampled for parasites were sampled for pathogens
## path.only.either <- path.only[!is.na(path.only$B.actin) | !is.na(path.only$Phorid),]

#combo par and path
par.and.path <- path.only

#par only and path only
par.only <- path.only
par.only$CBPV <- NULL
par.only$DWV_KV_VDV <- NULL
par.only$ABPV_KBV_IAPV <- NULL
par.only$BQCV <- NULL
par.only$SBPV <- NULL
par.only$SBV <- NULL
par.only$B.actin <- NULL

path.only$Phorid <- NULL
path.only$Crithidia <- NULL
path.only$Apicystis <- NULL

## make sure each dataframe has samples screened for every single par/path, 
path.only <- path.only[!is.na(path.only$B.actin),]
path.only <- path.only[!is.na(path.only$CBPV),]
path.only <- path.only[!is.na(path.only$DWV_KV_VDV),]
path.only <- path.only[!is.na(path.only$ABPV_KBV_IAPV),]
path.only <- path.only[!is.na(path.only$BQCV),]
path.only <- path.only[!is.na(path.only$SBPV),]
path.only <- path.only[!is.na(path.only$SBV),]

par.only <- par.only[!is.na(par.only$Phorid),]
par.only <- par.only[!is.na(par.only$Crithidia),]
par.only <- par.only[!is.na(par.only$Apicystis),]

par.and.path$B.actin <- NULL
par.and.path <- par.and.path[!is.na(par.and.path$CBPV),]
par.and.path <- par.and.path[!is.na(par.and.path$DWV_KV_VDV),]
par.and.path <- par.and.path[!is.na(par.and.path$ABPV_KBV_IAPV),]
par.and.path <- par.and.path[!is.na(par.and.path$BQCV),]
par.and.path <- par.and.path[!is.na(par.and.path$SBPV),]
par.and.path <- par.and.path[!is.na(par.and.path$SBV),]
par.and.path<- par.and.path[!is.na(par.and.path$Phorid),]
par.and.path <- par.and.path[!is.na(par.and.path$Crithidia),]
par.and.path <- par.and.path[!is.na(par.and.path$Apicystis),]


## drop B.actin equals zero, i.e. 
path.only <- path.only[!path.only$B.actin == 0,]

## community health metrics
parasites <- c("Phorid", "Crithidia", "Apicystis")

pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV", "BQCV","SBPV", "SBV")

parandpath <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV", "BQCV","SBPV", "SBV","Phorid", "Crithidia", "Apicystis")

par.only$ParasiteRichness <- rowSums(par.only[,parasites],na.rm=TRUE)
par.only$PossibleParasite <- apply(par.only[,parasites],1,function(x) sum(!is.na(x)))
par.only$AnyParasite <- (par.only$ParasiteRichness > 0)*1

path.only$PathogenRichness <- rowSums(path.only[,pathogens],na.rm=TRUE)
path.only$PossiblePathogen <- apply(path.only[,pathogens],1,function(x) sum(!is.na(x)))
path.only$AnyPathogen <- (path.only$PathogenRichness > 0)*1

par.and.path$ParPathRichness <- rowSums(par.and.path[,parandpath ],na.rm=TRUE)
par.and.path$PossibleParPath <- apply(par.and.path[,parandpath ],1,function(x) sum(!is.na(x)))
par.and.path$AnyParPath <- (par.and.path$ParPathRichnes > 0)*1


#potentially useful variable
par.only$ParasiteRichnessRate <- par.only$ParasiteRichness/par.only$PossibleParasite
path.only$PathogenRichnessRate <- path.only$PathogenRichness/path.only$PossiblePathogen


## *************************************************************


print("Bombus par/path richness")
table(par.and.path$ParPathRichness[par.and.path$Genus == "Bombus"])/nrow(par.and.path[par.and.path$Genus == "Bombus",])

print("Apis par/path richness")
table(par.and.path$ParPathRichness[par.and.path$Genus == "Bombus"])/nrow(par.and.path[par.and.path$Genus == "Apis",])

print("ParPath Bombus")
colSums(par.and.path[par.and.path$Genus ==  "Bombus", parandpath])/
  sum(par.and.path$Genus ==  "Bombus")

print("ParPath Apis")
colSums(par.and.path[par.and.path$Genus ==  "Apis", parandpath])/
  sum(par.and.path$Genus ==  "Apis")

print("Bombus parasite richness")
table(par.only$ParasiteRichness[par.only$Genus == "Bombus"])/nrow(par.only[par.only$Genus == "Bombus",])

print("Apis parasite richness")
table(par.only$ParasiteRichness[par.only$Genus == "Apis"])/nrow(par.only[par.only$Genus == "Apis",])

print("Parasites Bombus")
colSums(par.only[par.only$Genus ==  "Bombus", parasites])/
    sum(par.only$Genus ==  "Bombus")

print("Parasites Apis")
colSums(par.only[par.only$Genus ==  "Apis", parasites])/
    sum(par.only$Genus ==  "Apis")


print("Bombus pathogen richness")
table(path.only$PathogenRichness[path.only$Genus == "Bombus"])/
    nrow(path.only[path.only$Genus == "Bombus",])

print("Apis pathogen richness")
table(path.only$PathogenRichness[path.only$Genus == "Apis"])/
    nrow(path.only[path.only$Genus == "Apis",])

print("Pathogens Bombus")
colSums(path.only[path.only$Genus ==  "Bombus", pathogens], na.rm=TRUE)/
    sum(path.only$Genus ==  "Bombus" & !is.na(path.only$CBPV))

print("Pathogens Apis")
colSums(path.only[path.only$Genus ==  "Apis", pathogens], na.rm=TRUE)/
    sum(path.only$Genus ==  "Apis" & !is.na(path.only$CBPV))


## *************************************************************
## fix IDs
bees$genus_sub_sp[bees$genus_sub_sp ==
                  "Agapostemon coloradinus"]  <- "Agapostemon texanus"
bees$genus_sub_sp[bees$genus_sub_sp ==
                  "Peponapis sp."]  <- "Peponapis pruinosa"
bees$genus_sub_sp[bees$genus_sub_sp ==
                  "Halictus sp."]  <- "Halictus tripartitus"
bees$genus_sub_sp[bees$genus_sub_sp ==
                  "Anthophora sp."]  <- "Anthophora urbana"

## fix site names to match between datasets
bees$site[bees$site ==
          "Coyote Creek"]  <- "CoyoteCreek"
bees$site[bees$site ==
          "Beach Flats"]  <- "BeachFlats"
bees$site[bees$site ==
          "La Colina"]  <- "LaColina"
bees$site[bees$site ==
          "MEarth"]  <- "Mearth"
bees$site[bees$site ==
          "Laguna Seca"]  <- "LagunaSeca"
bees$site[bees$site ==
          "Prusch"]  <- "PruschPark"

## *************************************************************
## calculate site level characteristics for bees
## abundance (average between sample rounds)


abund.SR.bees <- bees %>%
    group_by(site, sample_pd) %>%
    summarise(BeeAbund = sum(no_individuals),
              BeeRichness = length(unique(genus_sub_sp)),
              BeeDiversity = diversity(table(genus_sub_sp)))

abund.bees <- abund.SR.bees %>%
    group_by(site) %>%
    summarise(BeeAbund = mean(BeeAbund),
              BeeRichness = mean(BeeRichness),
              BeeDiversity = mean(BeeDiversity))

abund.SR.bees.genus <- bees %>%
    group_by(site, sample_pd, genus) %>%
    summarise(BeeAbund = sum(no_individuals))

abund.bees.genus <- abund.SR.bees.genus %>%
    group_by(site, genus) %>%
    summarise(BeeAbund = mean(BeeAbund))


abund.bees.apis <- abund.bees.genus[abund.bees.genus$genus == "Apis",]
colnames(abund.bees.apis)[colnames(abund.bees.apis) == "BeeAbund"]  <-
    "ApisAbund"
abund.bees.bombus <- abund.bees.genus[abund.bees.genus$genus ==
                                      "Bombus",]
colnames(abund.bees.bombus)[colnames(abund.bees.bombus) == "BeeAbund"]  <-
    "BombusAbund"

abund.bees.apis$genus <- NULL
abund.bees.bombus$genus <- NULL


site.bees  <- abund.bees
site.bees <- merge(site.bees, abund.bees.apis,  all.x=TRUE)
site.bees <- merge(site.bees, abund.bees.bombus, all.x=TRUE)
site.bees$BombusAbund[is.na(site.bees$BombusAbund)] <- 0
site.bees$ApisAbund[is.na(site.bees$ApisAbund)] <- 0


## *************************************************************
## calculate site level characteristics for veg and merge with bee data

veg.site <- aggregate(list(AbundWoodyFlowers=veg$NoTreeShrubsFlower,
                           AbundAnnualFlowers=veg$NumTotalFlowers,
                           RichAnnualFlowers=veg$NumHerbPlantSpp,
                           RichWoodyFlowers=veg$NoTreeShrubSpp,
                           PlantRichness=veg$NumHerbPlantSpp,
                           PercentBareSoil=veg$PercentBareSoil,
                           Size=veg$Size,
                           natural1000m=veg$natural1000m),
                      list(site=veg$Site),
                      mean)

site.char <- merge(veg.site, site.bees)


## *************************************************************
## calculate densities to control for garden size (removed from analyses upon reviews)

# site.char$BeeDensity <- site.char$BeeAbund/site.char$Size
# site.char$BeeRichnessArea <- site.char$BeeRichness/site.char$Size
# 
# site.char$WoodyFlowerDensity <- site.char$AbundWoodyFlowers/site.char$Size
# site.char$AnnualFlowerDensity <-
#     site.char$AbundAnnualFlowers/site.char$Size
# site.char$PlantRichnessArea <- site.char$PlantRichness/site.char$Size

## *************************************************************
## calculate total sick individuals for each site with parasites, bad thing

par.totals <- aggregate(par.only[c(parasites)],
                         list(Site=par.only$Site,
                              Genus=par.only$Genus),
                         sum, na.rm=TRUE)

tested.totals <- aggregate(par.only[c(parasites)],
                           list(Site=par.only$Site,
                                Genus=par.only$Genus),
                           function(x) sum(!is.na(x)))

par.totals$ScreenedPar <- tested.totals$Phorid

par.totals[,parasites] <-
    par.totals[,parasites]/par.totals$ScreenedPar


##  calculate total sick individuals for each site with pathogens, bad thing

sick.totals <- aggregate(path.only[c(pathogens)],
                         list(Site=path.only$Site,
                              Genus=path.only$Genus),
                         sum, na.rm=TRUE)

tested.totals <- aggregate(path.only[c(pathogens)],
                           list(Site=path.only$Site,
                                Genus=path.only$Genus),
                           function(x) sum(!is.na(x)))

sick.totals$ScreenedPath <- tested.totals$CBPV


sick.totals[,pathogens] <-
  sick.totals[,pathogens]/sick.totals$ScreenedPath


## calculate total sick individuals for each site with parasites and pathogens, bad thing

parpath.totals <- aggregate(par.and.path[c(parandpath)],
                        list(Site=par.and.path$Site,
                             Genus=par.and.path$Genus),
                        sum, na.rm=TRUE)

tested.totals <- aggregate(par.and.path[c(parandpath)],
                           list(Site=par.and.path$Site,
                                Genus=par.and.path$Genus),
                           function(x) sum(!is.na(x)))

parpath.totals$ScreenedParPath <- tested.totals$Phorid

parpath.totals[,parandpath] <-
  parpath.totals[,parandpath]/parpath.totals$ScreenedParPath
## *************************************************************
## calculate parasite/pathogen rates in apis and parsaite/pathogen rates in bombus

par.path.rates <- par.and.path%>%
  group_by(Site, Genus) %>%
  summarise(parpath.rates = mean(AnyParPath))

par.rates <- par.only %>%
  group_by(Site, Genus) %>%
  summarise(par.rates = mean(AnyParasite))

path.rates <- path.only %>%
  group_by(Site, Genus) %>%
  summarise(path.rates = mean(AnyPathogen))

names(par.rates)[names(par.rates) == "Site"] <- "site"
names(path.rates)[names(path.rates) == "Site"] <- "site"
names(par.path.rates)[names(par.path.rates) == "Site"] <- "site"

apis.parpath.rate <- par.path.rates[par.path.rates$Genus == "Apis",]
bombus.parpath.rate <- par.path.rates[par.path.rates$Genus == "Bombus",]
apis.par.rate <- par.rates[par.rates$Genus == "Apis",]
bombus.par.rate <- par.rates[par.rates$Genus == "Bombus",]
apis.path.rate <- path.rates[path.rates$Genus == "Apis",]
bombus.path.rate <- path.rates[path.rates$Genus == "Bombus",]

colnames(apis.parpath.rate)[colnames(apis.parpath.rate) == "parpath.rates"]  <-
  "apis.parpath.rate"
colnames(bombus.parpath.rate)[colnames(bombus.parpath.rate) == "parpath.rates"]  <-
  "bombus.parpath.rate"
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
apis.parpath.rate$Genus <- NULL
bombus.parpath.rate$Genus <- NULL

site.char <- merge(site.char, apis.par.rate,  all.x=TRUE)
site.char <- merge(site.char, bombus.par.rate,  all.x=TRUE)
site.char <- merge(site.char, apis.path.rate,  all.x=TRUE)
site.char <- merge(site.char, bombus.path.rate,  all.x=TRUE)
site.char <- merge(site.char, apis.parpath.rate,  all.x=TRUE)
site.char <- merge(site.char, bombus.parpath.rate,  all.x=TRUE)

site.char$apis.par.rate[is.na(site.char$apis.par.rate)] <- 0
site.char$bombus.par.rate[is.na(site.char$bombus.par.rate)] <- 0
site.char$apis.path.rate[is.na(site.char$apis.path.rate)] <- 0
site.char$bombus.path.rate[is.na(site.char$bombus.path.rate)] <- 0
site.char$apis.parpath.rate[is.na(site.char$apis.parpath.rate)] <- 0
site.char$bombus.parpath.rate[is.na(site.char$bombus.parpath.rate)] <- 0


## *************************************************************

## We ended up redoing this in 2CommunityHealth so dont need it here

## standardize varaibles
# path.variables <- c("WoodyFlowerDensity", "AbundWoodyFlowers",
#                     "AnnualFlowerDensity", "AbundAnnualFlowers",
#                     "PlantRichness",
#                     "PlantRichnessArea",
#                     "natural1000m",
#                     "natural2000m",
#                     "PercentBareSoil",
#                     "BeeDensity",
#                     "BeeAbund",
#                     "ApisAbund",
#                     "BombusAbund",
#                     "BeeRichness",
#                     "Size")
# 
# standardize <- function(x)
# (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
# 
# site.char[,path.variables] <- apply(site.char[,path.variables], 2,
#                                     standardize)

## *************************************************************
## merge parasite and pathogen and site data into specimens.complete

par.and.path$Site[!par.and.path$Site %in% site.char$Site]
par.and.path$Date <- NULL

## merge pathogen and site data first
path.only$Site[!path.only$Site %in% site.char$Site]

dim(path.only)
path.only$Date <- NULL
dim(path.only)

## need to change "site" to "Site" in site.char
colnames(site.char)[colnames(site.char) == 'site'] <- 'Site'

path.only <- merge(path.only, site.char)
dim(path.only)

par.and.path <- merge(par.and.path, site.char)
dim(par.and.path)

## merge parasite and site data 
par.only$Site[!par.only$Site %in% site.char$Site]

dim(par.only)
par.only$Date <- NULL
par.only <- merge(par.only, site.char)
dim(par.only)

sick.totals <- merge(sick.totals, site.char)

## write out pathogen data
write.csv(path.only, file=file.path(save.dir,
                                   "specimens-completePathogen.csv"), row.names=FALSE)
## write out parasite data
write.csv(par.only, file=file.path(save.dir,
                                   "specimens-completeParasite.csv"), row.names=FALSE)

write.csv(par.and.path, file=file.path(save.dir,
                                   "specimens-completeParPath.csv"), row.names=FALSE)

## save.dir.git <- "/Volumes/Mac\ 2/Dropbox/urbanbeeparasites/data"

save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(path.only, par.only, par.and.path,
     site.char, sick.totals,
     file=file.path(save.dir.git, "specimens-complete.Rdata"))



## *************************************************************
## ## species accumulation curves by site
## library(vegan)

## ## by site make bee comm matrix, with columns = taxa, rows = site, use
## ## tapply from base r

## beecomm <- with(bees, tapply(no_individuals,
##                              list(site, genus_sub_sp),
##                              FUN = mean))
## beecomm[is.na(beecomm)] <- 0

## beeaccum <- specaccum(beecomm, method = "random",
##                       permutations = 999, conditioned =TRUE,
##                       gamma = "jack1",  w = NULL)

## plot(beeaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
## boxplot(beeaccum, col="yellow", add=TRUE, pch="+")

## ********************************************************

## ## species accumulation curves by period
## library(vegan)

## ## by site make bee comm matrix, with columns = taxa, rows = period, use
## ## tapply from base r

## beecomm <- with(bees, tapply(no_individuals,
##                              list(sample_pd, genus_sub_sp),
##                              FUN = mean))
## beecomm[is.na(beecomm)] <- 0

## beeaccum <- specaccum(beecomm, method = "random",
##                       permutations = 999, conditioned =TRUE,
##                       gamma = "jack1",  w = NULL)

## plot(beeaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
## boxplot(beeaccum, col="yellow", add=TRUE, pch="+")


## *************************************************************
## are par and path in apis related to par and path in bombus
## *************************************************************

# bombusPath <- path.only[path.only$Genus == "Bombus",]
# apisPath <- path.only[path.only$Genus == "Apis",]
# 
# bombusPara <- par.only[par.only$Genus == "Bombus",]
# apisPara <- par.only[par.only$Genus == "Apis",]
# 
# 
# bombusPath.site <- aggregate(list(bombusPathRichness=bombusPath$PathogenRichness,
#                                  bombusAnyPathgen=bombusPath$AnyPathogen),
#                       list(site=bombusPath$Site),
#                       mean)
# 
# bombusPara.site <- aggregate(list(bombusParRichness=bombusPara$ParasiteRichness,
#                                  bombusAnyParasite=bombusPara$AnyParasite),
#                             list(site=bombusPara$Site),
#                             mean)
# 
# apisPath.site <- aggregate(list( apisPathRichness=apisPath$PathogenRichness,
#                                   apisAnyPathgen=apisPath$AnyPathogen),
#                              list(site=apisPath$Site),
#                              mean)
# 
# apisPara.site <- aggregate(list( apisParRichness=apisPara$ParasiteRichness,
#                                   apisAnyParasite=apisPara$AnyParasite),
#                              list(site=apisPara$Site),
#                              mean)
# 
# apisPath.site <- apisPath.site[-c(14), ] 
# apisPara.site <- apisPara.site[-c(14), ] 
# 
# PathogenLM <- lm(bombusPath.site$bombusPathRichness ~ apisPath.site$apisPathRichness)
# summary(PathogenLM)
# 
# ParasiteLM <- lm(bombusPara.site$bombusParRichness ~ apisPara.site$apisParRichness)
# summary(ParasiteLM)
# 
# PathogenLM2 <- lm(bombusPath.site$bombusAnyPathgen ~ apisPath.site$apisAnyPathgen)
# summary(PathogenLM2)
# 
# ParasiteLM2 <- lm(bombusPara.site$bombusAnyParasite ~ apisPara.site$apisAnyParasite)
# summary(ParasiteLM2)
# 
# 
