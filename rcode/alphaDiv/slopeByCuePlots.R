# Started June 17, 2025 by D. Loughnan
# Aim of this code is to make plots for M2M analysis using the output from the joint Bayesian model
# cue x slope plot 

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)
require(shinystan)

setwd("~/Documents/github/Calvert_O-Connor_eelgrass")

sitelist <- sort(unique(alphaFour$site))

load("m2m/output/mdlShan184Param.Rdata")

ModelFit <- rstan::extract(shan184P)

muStG <- data.frame(ModelFit$mu_grand) 
muSt <- data.frame(ModelFit$alphamuSt)
muGSt <- muStG[,1] + muSt

muStMean <- colMeans(muGSt)

betaTempSt <- data.frame(ModelFit$betatempSt)
betaTempStMean <- colMeans(betaTempSt)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaTempSt, 2, quantile2575) 
muSt_quan <- apply(muGSt, 2, quantile2575)

bfs <- rbind(betaTempStMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "Temp25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "Temp75"
bfs_df$Sites <- sitelist

mg<- rbind(muStMean, muSt_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "micro25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "micro75"
mg_df$Sites <- sitelist

muTempSt <- data.frame(ModelFit$mutempSt)
muTempStMean <- colMeans(muTempSt)
a5 <- quantile(muTempSt$ModelFit.mutempSt, prob = c(0.25))
a95 <- quantile(muTempSt$ModelFit.mutempSt, prob = c(0.75))
apoly <- subset(muTempSt, ModelFit.mutempSt > a5 & ModelFit.mutempSt < a95)

betaMicroxTemp<- data.frame(ModelFit$betaMicroxtemp)
b5 <- quantile(betaMicroxTemp$ModelFit.betaMicroxtemp, prob = c(0.25))
b95 <- quantile(betaMicroxTemp$ModelFit.betaMicroxtemp, prob = c(0.75))
bpolly <- subset(betaMicroxTemp, ModelFit.betaMicroxtemp > b5 & ModelFit.betaMicroxtemp < b95)
betaMicroxTempMean <- colMeans(betaMicroxTemp)

pdf("m2m/figures/cueSlopes18S.pdf", height = 4, width = 8)
par(mar = c(5, 5, 4, 2), mfrow = c(1,2))
plot( x= mg_df$muStMean, y = bfs_df$betaTempStMean, 
      type="n", xlim = c(min(mg_df$micro25), max(mg_df$micro75)), 
      ylim = c(-1,1), # c(min(bfs_df$Temp25), max(bfs_df$Temp75)), 
      ylab = "Site-level temperature slope", xlab = "Microbial diversity", cex.lab = 2.25, cex.axis = 2) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
abline(a= a5, b= b5, col = "pink")
abline(a=muTempStMean, b=betaMicroxTempMean, col = "black", lwd = 2)
abline(a= a95, b= b95, col = "pink")

arrows(
  mg_df[,"muStMean"], # x mean
  bfs_df[,"Temp25"], # y 25
  mg_df[,"muStMean"],
  bfs_df[,"Temp75"],
  length = 0, col= "#218380", lwd = 2
)

arrows(
  mg_df[,"micro25"], # x mean
  bfs_df[,"betaTempStMean"], # y 25
  mg_df[,"micro75"], # x mean
  bfs_df[,"betaTempStMean"],
  length = 0, col = "#218380", lwd = 2
)

#------------------------------------------------------------------------------#
betaBiomassSt <- data.frame(ModelFit$betabiomassSt)
betaBiomassStMean <- colMeans(betaBiomassSt)
bc_quan <- apply(betaBiomassSt, 2, quantile2575)

bcs <- rbind(betaBiomassStMean, bc_quan)
bcs_t <- t(bcs)
bcs_df <- data.frame(bcs_t)
colnames(bcs_df)[colnames(bcs_df) == "X25."] <- "Biomass25"
colnames(bcs_df)[colnames(bcs_df) == "X75."] <- "Biomass75"
bcs_df$Sites <- sitelist

muBiomassSt <- data.frame(ModelFit$mubiomassSt)
muBiomassStMean <- colMeans(muBiomassSt)
a5 <- quantile(muBiomassSt$ModelFit.mubiomassSt, prob = c(0.25))
a95 <- quantile(muBiomassSt$ModelFit.mubiomassSt, prob = c(0.75))
apoly <- subset(muBiomassSt, ModelFit.mubiomassSt > a5 & ModelFit.mubiomassSt < a95)

betaMicroxbiomass<- data.frame(ModelFit$betaMicroxbiomass)
b5 <- quantile(betaMicroxbiomass$ModelFit.betaMicroxbiomass, prob = c(0.25))
b95 <- quantile(betaMicroxbiomass$ModelFit.betaMicroxbiomass, prob = c(0.75))
bpolly <- subset(betaMicroxbiomass, ModelFit.betaMicroxbiomass > b5 & ModelFit.betaMicroxbiomass < b95)
betaMicroxbiomassMean <- colMeans(betaMicroxbiomass)


plot( x= mg_df$muStMean, y = bcs_df$betaBiomassStMean, type="n", xlim = c(min(mg_df$micro25), max(mg_df$micro75)), 
      ylim = c(-0.3, 0.3)#c(min(bcs_df$Biomass25), max(bcs_df$Biomass75))
      , ylab = "Sites-level Biomass slope", xlab = "Microbial diversity", cex.lab = 2.25, cex.axis = 2) # blank plot with x range 

# mtext(side = 3, text = "Biomass", adj = 0, cex = 3)
# for(j in 1:length(apoly[,1])){
#   abline(a = apoly[j,], b = bpolly[j,], col= "#E0E0E0")
# }
# for(j in 1:length(muBiomassSt[,1])){
#   abline(a = muBiomassSt[j,], b = betaMicroxbiomassMean, col=alpha("#73d2de", 0.085))
# }
abline(a=muBiomassStMean, b=betaMicroxbiomassMean, col = "black", lwd =2)
abline(a= a5, b= b5, col = "black")
abline(a= a95, b= b95, col = "black")

arrows(
  mg_df[,"muStMean"], # x mean
  bcs_df[,"Biomass25"], # y 25
  mg_df[,"muStMean"],
  bcs_df[,"Biomass75"],
  length = 0, col= "#218380", lwd = 2
)

arrows(
  mg_df[,"micro25"], # x mean
  bcs_df[,"betaBiomassStMean"], # y 25
  mg_df[,"micro75"], # x mean
  bcs_df[,"betaBiomassStMean"],
  length = 0, col = "#218380", lwd = 2
)

dev.off()
###############
load("m2m/output/mdlShan16FourParaSimple.Rdata")

ModelFit <- rstan::extract(shan164P)

muStG <- data.frame(ModelFit$mu_grand) 
muSt <- data.frame(ModelFit$alphamuSt)
muGSt <- muStG[,1] + muSt

muStMean <- colMeans(muGSt)

betaTempSt <- data.frame(ModelFit$betatempSt)
betaTempStMean <- colMeans(betaTempSt)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaTempSt, 2, quantile2575) 
muSt_quan <- apply(muGSt, 2, quantile2575)

bfs <- rbind(betaTempStMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "Temp25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "Temp75"
bfs_df$Sites <- sitelist

mg<- rbind(muStMean, muSt_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "micro25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "micro75"
mg_df$Sites <- sitelist

muTempSt <- data.frame(ModelFit$mutempSt)
muTempStMean <- colMeans(muTempSt)
a5 <- quantile(muTempSt$ModelFit.mutempSt, prob = c(0.25))
a95 <- quantile(muTempSt$ModelFit.mutempSt, prob = c(0.75))
apoly <- subset(muTempSt, ModelFit.mutempSt > a5 & ModelFit.mutempSt < a95)

betaMicroxTemp<- data.frame(ModelFit$betaMicroxtemp)
b5 <- quantile(betaMicroxTemp$ModelFit.betaMicroxtemp, prob = c(0.25))
b95 <- quantile(betaMicroxTemp$ModelFit.betaMicroxtemp, prob = c(0.75))
bpolly <- subset(betaMicroxTemp, ModelFit.betaMicroxtemp > b5 & ModelFit.betaMicroxtemp < b95)
betaMicroxTempMean <- colMeans(betaMicroxTemp)

pdf("m2m/figures/cueSlope16S.pdf", height = 4, width = 8)
par(mar = c(5, 5, 4, 2), mfrow = c(1,2))
plot( x= mg_df$muStMean, y = bfs_df$betaTempStMean, 
      type="n", xlim = c(min(mg_df$micro25), max(mg_df$micro75)), 
      ylim = c(-1,1), # c(min(bfs_df$Temp25), max(bfs_df$Temp75)), 
      ylab = "Site-level temperature slope", xlab = "Microbial diversity") # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
abline(a= a5, b= b5, col = "pink")
abline(a=muTempStMean, b=betaMicroxTempMean, col = "black", lwd =2)
abline(a= a95, b= b95, col = "pink")

arrows(
  mg_df[,"muStMean"], # x mean
  bfs_df[,"Temp25"], # y 25
  mg_df[,"muStMean"],
  bfs_df[,"Temp75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df[,"micro25"], # x mean
  bfs_df[,"betaTempStMean"], # y 25
  mg_df[,"micro75"], # x mean
  bfs_df[,"betaTempStMean"],
  length = 0, col = "maroon", lwd = 2
)

#------------------------------------------------------------------------------#
betaBiomassSt <- data.frame(ModelFit$betabiomassSt)
betaBiomassStMean <- colMeans(betaBiomassSt)
bc_quan <- apply(betaBiomassSt, 2, quantile2575)

bcs <- rbind(betaBiomassStMean, bc_quan)
bcs_t <- t(bcs)
bcs_df <- data.frame(bcs_t)
colnames(bcs_df)[colnames(bcs_df) == "X25."] <- "Biomass25"
colnames(bcs_df)[colnames(bcs_df) == "X75."] <- "Biomass75"
bcs_df$Sites <- sitelist

muBiomassSt <- data.frame(ModelFit$mubiomassSt)
muBiomassStMean <- colMeans(muBiomassSt)
a5 <- quantile(muBiomassSt$ModelFit.mubiomassSt, prob = c(0.25))
a95 <- quantile(muBiomassSt$ModelFit.mubiomassSt, prob = c(0.75))
apoly <- subset(muBiomassSt, ModelFit.mubiomassSt > a5 & ModelFit.mubiomassSt < a95)

betaMicroxbiomass<- data.frame(ModelFit$betaMicroxbiomass)
b5 <- quantile(betaMicroxbiomass$ModelFit.betaMicroxbiomass, prob = c(0.25))
b95 <- quantile(betaMicroxbiomass$ModelFit.betaMicroxbiomass, prob = c(0.75))
bpolly <- subset(betaMicroxbiomass, ModelFit.betaMicroxbiomass > b5 & ModelFit.betaMicroxbiomass < b95)
betaMicroxbiomassMean <- colMeans(betaMicroxbiomass)


plot( x= mg_df$muStMean, y = bcs_df$betaBiomassStMean, type="n", xlim = c(min(mg_df$micro25), max(mg_df$micro75)), 
      ylim = c(-0.3, 0.3)#c(min(bcs_df$Biomass25), max(bcs_df$Biomass75))
      , ylab = "Sites-level Biomass slope", xlab = "Microbial diversity") # blank plot with x range 

# mtext(side = 3, text = "Biomass", adj = 0, cex = 3)
# for(j in 1:length(apoly[,1])){
#   abline(a = apoly[j,], b = bpolly[j,], col= "#E0E0E0")
# }
# for(j in 1:length(muBiomassSt[,1])){
#   abline(a = muBiomassSt[j,], b = betaMicroxbiomassMean, col=alpha("#73d2de", 0.085))
# }
abline(a=muBiomassStMean, b=betaMicroxbiomassMean, col = "black", lwd =2)
abline(a= a5, b= b5, col = "black")
abline(a= a95, b= b95, col = "black")

arrows(
  mg_df[,"muStMean"], # x mean
  bcs_df[,"Biomass25"], # y 25
  mg_df[,"muStMean"],
  bcs_df[,"Biomass75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df[,"micro25"], # x mean
  bcs_df[,"betaBiomassStMean"], # y 25
  mg_df[,"micro75"], # x mean
  bcs_df[,"betaBiomassStMean"],
  length = 0, col = "maroon", lwd = 2
)

dev.off()