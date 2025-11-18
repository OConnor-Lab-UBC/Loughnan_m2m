# Started June 17, 2025 by D. Loughnan
# Aim of this code is to make plots for M2M analysis using the output from the joint Bayesian model
# muPlots and plots of site-level differences

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)
require(shinystan)

setwd("~/Documents/github/Calvert_O-Connor_eelgrass")

load("m2m/output/mdlShan16FourPara.Rdata")
load("m2m/output/mdlShan184Param.Rdata")

sumer18 <- summary(mdlShan18)$summary[c(
  "mu_grand","sigmaSt","muyear", "sigmayear", "sigma_microy", 
  "mu_grandM", "muyear2", "sigmayear2",
  "mutempSt","sigmatempSt","mubiomassSt","sigmabiomassSt",
  "muDepth", "sigmaDepth", "mubedA","sigmabedA",
  "sigmaMacroSt","sigmaMacro_y",
  "betaMicroxtemp","betaMicroxbiomass"),
  c("mean","25%", "75%","n_eff", "Rhat")]

fitSh18 <- data.frame(rstan::extract(mdlShan18))
fitSh18 <- fitSh18[, c("mu_grand","sigmaSt","muyear", "sigmayear", "sigma_microy", 
                       "mu_grandM", "muyear2", "sigmayear2",
                       "mutempSt","sigmatempSt","mubiomassSt","sigmabiomassSt",
                       "muDepth", "sigmaDepth", "mubedA","sigmabedA",
                       "sigmaMacroSt","sigmaMacro_y",
                       "betaMicroxtemp","betaMicroxbiomass")]

quantSh18 <- t(apply(fitSh18, 2, quantile2575) )

tabSh18 <- cbind(sumer18, quantSh18)

muSub18 <- tabSh18[c("mu_grand","muyear", 
                      "mu_grandM", "muyear2", 
                      "mutempSt","mubiomassSt",
                     "muDepth", "mubedA",
                     "betaMicroxtemp","betaMicroxbiomass"),]

pdf("..//figures/muPlot18S.pdf", width = 8, height = 5)
par(mfrow = c(1,1), mar = c(5,9, 2, 1))
plot(seq(-25,15,
         length.out = nrow(muSub18)), 
     1:nrow(muSub18),
     type = "n",
     xlab = "Parameter estimate",
     ylab = "",
     yaxt = "n")

points(muSub18[1:nrow(muSub18), 'mean'],
       1:nrow(muSub18),
       pch = 16,
       col = "cyan4",
       cex = 2)

axis(2, at = 1:nrow(muSub18), labels = rownames(muSub18), las = 1, cex.axis = 0.8)

arrows(muSub18[, "5%"], 1:nrow(muSub18), muSub18[, "95%"], 1:nrow(muSub18),
       len = 0, col = "black")
arrows(muSub18[, "25%"], 1:nrow(muSub18), muSub18[, "75%"], 1:nrow(muSub18),
       len = 0, col = "black", lwd = 3)

abline(v = 0, lty = 3)
dev.off()


colors <- c("#cc6a70ff", "maroon", "purple3","#593d9cff","#f9b641ff","#efe350ff","cyan4","#eb8055ff", "sienna")

# 