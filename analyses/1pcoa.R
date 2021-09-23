## setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
library(vegan)
source('src/misc.R')
source('src/pcoa.R')

load("../data/specimens-complete.Rdata")

parasites <- c("Phorid", "Crithidia", "Apicystis")
pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV",
               "BQCV","SBPV", "SBV")

GenSp <- par.path$Genus
Sites <- par.path$Site

parasite.comms <- calcPcoa(par.path, parasites, nperm=1000, GenSp,
                           Sites)
parasite.comms$tests

pathogens.comms <- calcPcoa(par.path, pathogens, nperm=1000, GenSp,
                            Sites)
pathogens.comms$tests


#$ plotting
plotCommDist(parasite.comms$dist$dist, parasite.comms$dist$sites,
             parasite.comms$dist$genus, "parasite")


plotCommDist(pathogens.comms$dist$dist, pathogens.comms$dist$sites,
             pathogens.comms$dist$genus, "pathogen")



##**************************R1 combiningg par and path***************************************

load("../data/specimens-complete.Rdata")

parandpath <- c("Phorid", "Crithidia", "Apicystis","CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV",
              "BQCV","SBPV", "SBV")

GenSp <- par.and.path$Genus
Sites <- par.and.path$Site

parandpath.comms <- calcPcoa(par.and.path, parandpath, nperm=1000, GenSp,
                           Sites)
parandpath.comms$tests



#$ plotting
plotCommDist(parandpath.comms$dist$dist, parandpath.comms$dist$sites,
             parandpath.comms$dist$genus, "parandpath1")

plotCommDist(parandpath.comms$dist$dist, parandpath.comms$dist$genus,
             parandpath.comms$dist$site, "parandpath2")






