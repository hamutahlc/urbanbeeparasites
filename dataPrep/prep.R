rm(list=ls())
## prepares raw data and creates dataset for analyses

save.dir <- "~/Dropbox/urbanbeeparasites_saved"
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
vosPP <- vosPP[,colnames(apisPP)]

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
abund.SR.bees <- aggregate(list(Abund=bees$no_individuals),
                           list(Site=bees$site,
                                SampleRound=bees$sample_pd),
                           sum)

abund.SR.bees <- tapply(abund.SR.bees$Abund,
                        abund.SR.bees$Site, mean)

## richness (total for a site)
site.bees <- aggregate(list(BeeRichness=bees$genus_sub_sp),
                       list(Site=bees$site),
                       function(x) length(unique(x)))

site.bees$BeeAbund <- abund.SR.bees[match(names(abund.SR.bees),
                                          site.bees$Site)]

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



dim(par.path)

## okay to drop controls
par.path$Site[!par.path$Site %in% site.char$Site]

par.path <- merge(par.path, site.char)
dim(par.path)

## *************************************************************
## calculate densities to control for garden size

par.path$BeeDensity <- par.path$BeeAbund/par.path$Size
par.path$BeeRichnessArea <- par.path$BeeRichness/par.path$Size

par.path$WoodyFlowerDensity <- par.path$AbundWoodyFlowers/par.path$Size
par.path$AnnualFlowerDenisty <-
    par.path$AbundAnnualFlowers/par.path$Size
par.path$PlantRichnessArea <- par.path$PlantRichness/par.path$Size


## write out final data
write.csv(par.path, file=file.path(save.dir,
                           "specimens-complete.csv"), row.names=FALSE)

save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(par.path, file=file.path(save.dir.git, "specimens-complete.Rdata"))
