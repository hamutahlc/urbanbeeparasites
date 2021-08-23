rm(list=ls())
## prepares raw data and creates dataset for analyses

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
par.path <- rbind(apisPP, vosPP)

## fix site names to match between datasets
par.path$Site <- as.character(par.path$Site)
par.path$Site[par.path$Site ==
              "Coyote"]  <- "CoyoteCreek"
par.path$Site[par.path$Site ==
              "Beach"]  <- "BeachFlats"
par.path$Site[par.path$Site ==
              "Beach Flats"]  <- "BeachFlats"
par.path$Site[par.path$Site ==
              "La Colina"]  <- "LaColina"
par.path$Site[par.path$Site ==
              "Laguna "]  <- "LagunaSeca"
par.path$Site[par.path$Site ==
              "Laguna"]  <- "LagunaSeca"
par.path$Site[par.path$Site ==
              "Prusch"]  <- "PruschPark"

## drop all the controls
par.path <- par.path[!par.path$Site == "Control",]

## i think judt want to make a separate parasite datasheet to use with the parasite data, 
## since the way the screenign occured is that a subset of bees sampled for parasites were sampled for pathogens
## par.path.either <- par.path[!is.na(par.path$B.actin) | !is.na(par.path$Phorid),]

par.only <- par.path
par.only$CBPV <- NULL
par.only$DWV_KV_VDV <- NULL
par.only$ABPV_KBV_IAPV <- NULL
par.only$BQCV <- NULL
par.only$SBPV <- NULL
par.only$SBV <- NULL
par.only$B.actin <- NULL

## par.path has only samples screened for every single par/path, 
par.path <- par.path[!is.na(par.path$B.actin),]
par.path <- par.path[!is.na(par.path$Phorid),]
par.path <- par.path[!is.na(par.path$Crithidia),]
par.path <- par.path[!is.na(par.path$Apicystis),]
par.path <- par.path[!is.na(par.path$CBPV),]
par.path <- par.path[!is.na(par.path$DWV_KV_VDV),]
par.path <- par.path[!is.na(par.path$ABPV_KBV_IAPV),]
par.path <- par.path[!is.na(par.path$BQCV),]
par.path <- par.path[!is.na(par.path$SBPV),]
par.path <- par.path[!is.na(par.path$SBV),]

## drop B.actin equals zero, i.e. 
par.path <- par.path[!par.path$B.actin == 0,]

## community health metrics
parasites <- c("Phorid", "Crithidia", "Apicystis")

pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV", "BQCV","SBPV", "SBV")

par.only$ParasiteRichness <- rowSums(par.only[,parasites],na.rm=TRUE)
par.only$PossibleParasite <- apply(par.only[,parasites],1,function(x) sum(!is.na(x)))
par.only$AnyParasite <- (par.only$ParasiteRichness > 0)*1

par.path$PathogenRichness <- rowSums(par.path[,pathogens],na.rm=TRUE)
par.path$PossiblePathogen <- apply(par.path[,pathogens],1,function(x) sum(!is.na(x)))
par.path$AnyPathogen <- (par.path$PathogenRichness > 0)*1

## *************************************************************
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

par.path$PathogenRichness <- rowSums(par.path[, pathogens])
par.path$AnyPathogen <- (par.path$PathogenRichness > 0)*1

print("Bombus pathogen richness")
table(par.path$PathogenRichness[par.path$Genus == "Bombus"])/
    nrow(par.path[par.path$Genus == "Bombus",])

print("Apis pathogen richness")
table(par.path$PathogenRichness[par.path$Genus == "Apis"])/
    nrow(par.path[par.path$Genus == "Apis",])

print("Pathogens Bombus")
colSums(par.path[par.path$Genus ==  "Bombus", pathogens], na.rm=TRUE)/
    sum(par.path$Genus ==  "Bombus" & !is.na(par.path$CBPV))

print("Pathogens Apis")
colSums(par.path[par.path$Genus ==  "Apis", pathogens], na.rm=TRUE)/
    sum(par.path$Genus ==  "Apis" & !is.na(par.path$CBPV))


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
              BeeRichness = length(unique(genus_sub_sp)))

abund.bees <- abund.SR.bees %>%
    group_by(site) %>%
    summarise(BeeAbund = mean(BeeAbund),
              BeeRichness = mean(BeeRichness))

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

## I TOOK THESE AND ADDED TO THE ABOVE...

## site.char <- merge(site.char,
 #                  unique(veg[, c("Site", "Size", "natural1000m",
   #                               "natural2000m")]))


## *************************************************************
## calculate densities to control for garden size

site.char$BeeDensity <- site.char$BeeAbund/site.char$Size
site.char$BeeRichnessArea <- site.char$BeeRichness/site.char$Size

site.char$WoodyFlowerDensity <- site.char$AbundWoodyFlowers/site.char$Size
site.char$AnnualFlowerDensity <-
    site.char$AbundAnnualFlowers/site.char$Size
site.char$PlantRichnessArea <- site.char$PlantRichness/site.char$Size

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

sick.totals <- aggregate(par.path[c(pathogens)],
                         list(Site=par.path$Site,
                              Genus=par.path$Genus),
                         sum, na.rm=TRUE)

tested.totals <- aggregate(par.path[c(pathogens)],
                           list(Site=par.path$Site,
                                Genus=par.path$Genus),
                           function(x) sum(!is.na(x)))

sick.totals$ScreenedPath <- tested.totals$CBPV


sick.totals[,pathogens] <-
  sick.totals[,pathogens]/sick.totals$ScreenedPath


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

## merge pathogen and site data first
par.path$Site[!par.path$Site %in% site.char$Site]

dim(par.path)
par.path$Date <- NULL
dim(par.path)

## need to change "site" to "Site" in site.char
colnames(site.char)[colnames(site.char) == 'site'] <- 'Site'

par.path <- merge(par.path, site.char)
dim(par.path)

## merge parasite and site data 
par.only$Site[!par.only$Site %in% site.char$Site]

dim(par.only)
par.only$Date <- NULL
par.only <- merge(par.only , site.char)
dim(par.only)

sick.totals <- merge(sick.totals, site.char)

## write out pathogen data
write.csv(par.path, file=file.path(save.dir,
                                   "specimens-completePathogen.csv"), row.names=FALSE)
## write out parasite data
write.csv(par.only, file=file.path(save.dir,
                                   "specimens-completeParasite.csv"), row.names=FALSE)

## save.dir.git <- "/Volumes/Mac\ 2/Dropbox/urbanbeeparasites/data"

save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(par.path, par.only,
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

