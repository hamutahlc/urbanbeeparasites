## setwd("~/Dropbox/urbanbeeparasites")
setwd("analyses")
rm(list=ls())
library(vegan)

load("../data/specimens-complete.Rdata")

parasites <- c("Phorid", "Crithidia", "Apicystis")
pathogens <- c("CBPV", "DWV_KV_VDV", "ABPV_KBV_IAPV", "BQCV","SBPV", "SBV")

GenSp <- par.path$Genus
Sites <- par.path$Site


calcPcoa <- function(par.path, paths.or.parsites){
    com <- par.path[, parasites]
    no.p <- rowSums(com) == 0
    com.dist  <- vegdist(com[!no.p,], method="jaccard")

    ## are the species (bombus vs. apis) different in parasite composition?
    beta.disper.gen <- betadisper(com.dist, GenSp[!no.p],
                                  type="centroid")

    perm.test.gen <- permutest(beta.disper.gen,
                               control = permControl(nperm = 100),
                               pairwise = TRUE)

    ## are the site different in parasite composition?
    beta.disper.site <- betadisper(com.dist, Sites[!no.p],
                                       type="centroid")

    perm.test.site <- permutest(beta.disper.site,
                                    control = permControl(nperm = 100),
                                    pairwise = TRUE)
    return(list(tests= list(site=perm.test.site,
                            species=perm.test.gen),
                dists=list(dist=com.dist, sites=Sites[!no.p],  genus=GenSp[!no.p])))
}


parasite.comms <- calcPcoa(par.path, parasites)
parasite.comms$tests

pathogens.comms <- calcPcoa(par.path, pathogens)
pathogens.comms$tests



plotCommDist  <- function(dist.mat, sites, genus, par.or.path){

    f.pcoa <- function(){
        cols <- rainbow(length(unique(sites)))
        names(cols) <- unique(sites)

        dist.mat <- as.matrix(dist.mat)
        pcoa.comm <- cmdscale(dist.mat)

        pcoa.mod <- adonis(dist.mat~sites)
        plot(NA, asp=1,  cex=1.5,
             ylim=range(pcoa.comm[,2]),
             xlim=range(pcoa.comm[,1]),
             xlab='',
             ylab='',
             xaxt='n',
             yaxt='n',
             cex.lab=1.5)
        ## pcoa.comm.jit <- jitter(jitter(pcoa.comm))
        for(s in unique(sites)){
            ## all points sitting on top of eachother so triple jitter
            points(pcoa.comm[sites == s,],
                   col=cols[s], pch=16, cex=1.5)
            points(pcoa.comm[sites == s,],
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
        return(pcoa.mod)
    }
    ## function for plotting PcoA axes
    path <- 'figures'
    pdf.f(f.pcoa, file= file.path(path,
                                  sprintf("%s_pcoa.pdf", par.or.path)),
          width=7, height=7)
}

pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3],
                  alpha=alpha))
}



plotCommDist(parasite.comms$dist$dist, parasite.comms$dist$sites,
             parasite.comms$dist$genus, "parasite")

plotCommDist(pathogens.comms$dist$dist, pathogens.comms$dist$sites,
             pathogens.comms$dist$genus, "pathogen")

