## setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
library(vegan)
source('src/misc.R')

load("../data/specimens-complete.Rdata")

parasites <- c("Phorid", "Crithidia", "Apicystis")
pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV",
               "BQCV","SBPV", "SBV")

GenSp <- par.path$Genus
Sites <- par.path$Site

calcPcoa <- function(par.path, paths.or.parasites, nperm,
                     genSp, sites){
    com <- par.path[, paths.or.parasites]
    no.p <- rowSums(com) == 0 | apply(com, 1, function(x) any(is.na(x)))
    com <- com[!no.p,]
    genSp <- as.character(genSp)[!no.p]
    sites <- as.character(sites)[!no.p]

    com.dist  <- vegdist(com, method="jaccard")

    ## are the species (bombus vs. apis) different in parasite composition?
    beta.disper.gen <- betadisper(com.dist, genSp,
                                  type="centroid")
    perm.test.gen <- permutest(beta.disper.gen,
                               control = permControl(nperm = nperm),
                               pairwise = TRUE)

    ## are the site different in parasite composition?
    beta.disper.site <- betadisper(com.dist, sites,
                                       type="centroid")
    perm.test.site <- permutest(beta.disper.site,
                                    control = permControl(nperm = nperm),
                                pairwise = TRUE)

    return(list(tests= list(site=perm.test.site,
                            species=perm.test.gen),
                dists=list(dist=com.dist, sites=sites,
                           genus=genSp)))
}

parasite.comms <- calcPcoa(par.path, parasites, nperm=100, GenSp,
                           Sites)
parasite.comms$tests

pathogens.comms <- calcPcoa(par.path, pathogens, nperm=100, GenSp,
                            Sites)
pathogens.comms$tests



plotCommDist  <- function(dist.mat, sites, genus, par.or.path){

    f.pcoa <- function(){
        cols <- rainbow(length(unique(sites)))
        names(cols) <- unique(sites)

        dist.mat <- as.matrix(dist.mat)
        pcoa.comm <- cmdscale(dist.mat)

        plot(NA, asp=1,  cex=1.5,
             ylim=range(pcoa.comm[,2]),
             xlim=range(pcoa.comm[,1]),
             xlab='',
             ylab='',
             xaxt='n',
             yaxt='n',
             cex.lab=1.5)
        for(site in unique(sites)){
            points(pcoa.comm[sites == site,],
                   col=cols[site], pch=16, cex=1.5)
            points(pcoa.comm[sites == site,],
                   col="black", pch=1, cex=1.5)
        }
        ordihull(pcoa.comm, sites)

        legend("topright", legend=unique(sites),
               bty="n", cex=0.6, col=cols[unique(sites)],
               pch=16)
        legend("topright", legend=unique(sites),
               bty="n", cex=0.6, col="black", pch=1)

        mtext('PCoA1', 1, line=2, cex=1.5)
        mtext('PCoA2', 2, line=2, cex=1.5)
    }
    ## function for plotting PcoA axes
    path <- 'figures'
    pdf.f(f.pcoa, file= file.path(path,
                                  sprintf("%s_pcoa.pdf", par.or.path)),
          width=7, height=7)
}


plotCommDist(parasite.comms$dist$dist, parasite.comms$dist$sites,
             parasite.comms$dist$genus, "parasite")


plotCommDist(pathogens.comms$dist$dist, pathogens.comms$dist$sites,
             pathogens.comms$dist$genus, "pathogen")

