rm(list=ls())
## prepares raw data and creates dataset for analyses

## save.dir <- "~/Dropbox/urbanbeeparasites_saved"
save.dir <- "/Volumes/Mac\ 2/Dropbox/urbanbeeparasites_saved"

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
apisPP$IDNum <- NULL
vosPP$Date <- NA
vosPP$Genus  <- "Bombus"
vosPP$Species  <- "vosnesenskii"
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

## community health metrics
parasites <- c("Phorid", "Crithidia", "Apicystis")

pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV", "BQCV","SBPV", "SBV")

par.path$ParasiteRichness <- rowSums(par.path[,parasites])
par.path$AnyParasite <- (par.path$ParasiteRichness > 0)*1

## *************************************************************
print("Bombus parasite richness")
table(par.path$ParasiteRichness[par.path$Genus == "Bombus"])/nrow(par.path[par.path$Genus == "Bombus",])

print("Apis parasite richness")
table(par.path$ParasiteRichness[par.path$Genus == "Apis"])/nrow(par.path[par.path$Genus == "Apis",])

print("Parasites Bombus")
colSums(par.path[par.path$Genus ==  "Bombus", parasites])/
    sum(par.path$Genus ==  "Bombus")

print("Parasites Apis")
colSums(par.path[par.path$Genus ==  "Apis", parasites])/
    sum(par.path$Genus ==  "Apis")

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
site.bees <- merge(site.bees, abund.bees.apis)
site.bees <- merge(site.bees, abund.bees.bombus, all.x=TRUE)
site.bees$BombusAbund[is.na(site.bees$BombusAbund)] <- 0

## *************************************************************
## ## species accumulation curves
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

## *************************************************************
## calculate site level characteristics for veg and merge with bee data

veg.site <- aggregate(list(AbundWoodyFlowers=veg$NoTreeShrubsFlower,
                           AbundAnnualFlowers=veg$NumTotalFlowers,
                           PlantRichness=veg$NumHerbPlantSpp,
                           PercentBareSoil=veg$PercentBareSoil),
                      list(Site=veg$Site),
                      mean)

site.char <- merge(veg.site, site.bees)

site.char <- merge(site.char,
                   unique(veg[, c("Site", "Size", "natural1000m",
                                  "natural2000m")]))

## *************************************************************
## calculate densities to control for garden size

site.char$BeeDensity <- site.char$BeeAbund/site.char$Size
site.char$BeeRichnessArea <- site.char$BeeRichness/site.char$Size

site.char$WoodyFlowerDensity <- site.char$AbundWoodyFlowers/site.char$Size
site.char$AnnualFlowerDensity <-
    site.char$AbundAnnualFlowers/site.char$Size
site.char$PlantRichnessArea <- site.char$PlantRichness/site.char$Size

## *************************************************************
## calculate total sick individuals for each site, bad thing

sick.totals <- aggregate(par.path[c(parasites, pathogens)],
                         list(Site=par.path$Site,
                              Genus=par.path$Genus),
                         sum, na.rm=TRUE)

tested.totals <- aggregate(par.path[c(parasites, pathogens)],
                           list(Site=par.path$Site,
                                Genus=par.path$Genus),
                           function(x) sum(!is.na(x)))

sick.totals$ScreenedPath <- tested.totals$CBPV
sick.totals$ScreenedPar <- tested.totals$Phorid

sick.totals[,parasites] <-
    sick.totals[,parasites]/sick.totals$ScreenedPar
sick.totals[,pathogens] <-
    sick.totals[,pathogens]/sick.totals$ScreenedPath

## *************************************************************
## standardize varaibles
path.variables <- c("WoodyFlowerDensity", "AbundWoodyFlowers",
                    "AnnualFlowerDensity", "AbundAnnualFlowers",
                    "PlantRichness",
                    "PlantRichnessArea",
                    "natural1000m",
                    "natural2000m",
                    "PercentBareSoil",
                    "BeeDensity",
                    "BeeAbund",
                    "ApisAbund",
                    "BombusAbund",
                    "BeeRichness")

standardize <- function(x)
(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

site.char[,path.variables] <- apply(site.char[,path.variables], 2,
                                    standardize)

## *************************************************************
## merge pathogen and site data

dim(par.path)
## okay to drop controls
par.path$Site[!par.path$Site %in% site.char$Site]

par.path <- merge(par.path, site.char)
par.path$Date <- NULL
dim(par.path)

sick.totals <- merge(sick.totals, site.char)

## write out final data
write.csv(par.path, file=file.path(save.dir,
                                   "specimens-complete.csv"), row.names=FALSE)

save.dir.git <- "/Volumes/Mac\ 2/Dropbox/urbanbeeparasites/data"

## save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(par.path,
     site.char, sick.totals,
     file=file.path(save.dir.git, "specimens-complete.Rdata"))

