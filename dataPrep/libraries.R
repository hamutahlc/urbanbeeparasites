rm(list=ls())
## install and load libraries for all analyses
install.packages("vegan")
install.packages("piecewiseSEM")
install.packages("lme4")
install.packages("lmerTest")


install.packages("R2admb")
install.packages("glmmADMB",
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")

install.packages("brms")

library(vegan)
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(brms)
