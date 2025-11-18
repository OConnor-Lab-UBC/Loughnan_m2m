# August 17, 2025 by D. Loughnan

# aim of this code is to create the code sourced in the manuscript:

library(rstan)
#library(shinystan)
#library(reshape2)
#library(bayesplot)
library(ggplot2)
# library(dplyr)
library(plyr)
library(stringr)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.25, 0.75,0.95))
  return(returnQuanilte)
}

### 16S prokaryotes

prok16 <- read.csv("..//..//Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
prok16 <- subset(prok16, site != "pruth_bay")
noProk <- ncol(prok16[,16:ncol(prok16)])
### 18S microeukaryotes ############

micro18 <- read.csv("..//..//Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
micro18 <- subset(micro18, site != "pruth_bay")

noEuk <- ncol(micro18[,10:ncol(micro18)])

### epifauna

epif <- read.csv("..//output/cleanedGrazer.csv")
epif <- subset(epif, site != "goose_north")

epiTot <- sum(epif$abundance)
noEpi <- length(unique(epif$taxonL))

epiYrSite <- aggregate(epif["abundance"], 
          epif[c("year","site")], 
          FUN = sum)
noEpiMin <- min(epiYrSite$abundance)
noEpiMax <- max(epiYrSite$abundance)

epiSiteMean <- aggregate(epiYrSite["abundance"], 
                         epiYrSite[c("site")], 
                       FUN = mean)

domSp <- aggregate(epif["abundance"], 
                   epif[c("taxonL")], 
                   FUN = sum)

rownames(domSp) <- domSp$taxonL
perDomTaxa <- round(((domSp["Caprella.californica",2] + 
                  domSp["Caprella.laeviuscula",2] +
                  domSp["Lacuna.vincta",2] +
                  domSp["Porcellidium",2] +
                  domSp["Gammaridea",2] 
                  )/epiTot)*100, 0)

noSpSite <- aggregate(epif["site"], 
                         epif[c("taxonL")], 
                         FUN = function(x)length(unique(x)))

noSpAll <- nrow(subset(noSpSite, site == "9"))
perSpMost  <- round((nrow(subset(noSpSite, site >4))/noEpi)*100,0)

## Model outputs:
load("..//output/mdlShan16ParamTLD_lumpedSp.Rdata")

# load("..//output/mdlShan16FourParaSimple.Rdata")
sumer16 <- summary(mdl4Shan16)$summary[c(
  "mu_grand","sigmaSt","muyear", "sigmayear", "sigma_microy", 
  "mu_grandM", "muyear2", "sigmayear2",
  "mutempSt","sigmatempSt","muLAISt","sigmaLAISt",
  "muDepth", 
  "sigmaMacroSt","sigmaMacro_y",
  "betaMicroxtemp","betaMicroxLAI"),
  c("mean","25%", "75%","n_eff", "Rhat")]

fitSh16 <- data.frame(rstan::extract(mdl4Shan16))
fitSh16 <- fitSh16[, c("mu_grand","sigmaSt","muyear", "sigmayear", "sigma_microy", 
                       "mu_grandM", "muyear2", "sigmayear2",
                       "mutempSt","sigmatempSt","muLAISt","sigmaLAISt",
                       "muDepth", 
                       "sigmaMacroSt","sigmaMacro_y",
                       "betaMicroxtemp","betaMicroxLAI")]


quantSh16 <- t(apply(fitSh16, 2, quantile2575) )

tabSh16 <- cbind(sumer16, quantSh16)

meanProkYr <- round(tabSh16["muyear","mean"],2)
lProkYr <- round(tabSh16["muyear","5%"],2)
uProkYr <- round(tabSh16["muyear","95%"],2)

meanEpiYr <- round(tabSh16["muyear2","mean"],2)
lEpiYr <- round(tabSh16["muyear2","5%"],2)
uEpiYr <- round(tabSh16["muyear2","95%"],2)

load("..//output/mdlShan18ParamTLD_lumpedSpp.Rdata")

sumer18 <- summary(mdl4Shan18)$summary[c(
  "mu_grand","sigmaSt","muyear", "sigmayear", "sigma_microy", 
  "mu_grandM", "muyear2", "sigmayear2",
  "mutempSt","sigmatempSt","muLAISt","sigmaLAISt",
  "muDepth", 
  "sigmaMacroSt","sigmaMacro_y",
  "betaMicroxtemp","betaMicroxLAI"),
  c("mean","25%", "75%","n_eff", "Rhat")]

fitSh18 <- data.frame(rstan::extract(mdl4Shan18))
fitSh18 <- fitSh18[, c( "mu_grand","sigmaSt","muyear", "sigmayear", "sigma_microy", 
                        "mu_grandM", "muyear2", "sigmayear2",
                        "mutempSt","sigmatempSt","muLAISt","sigmaLAISt",
                        "muDepth", 
                        "sigmaMacroSt","sigmaMacro_y",
                        "betaMicroxtemp","betaMicroxLAI")]

quantSh18 <- t(apply(fitSh18, 2, quantile2575) )

tabSh18 <- cbind(sumer18, quantSh18)


########################
m <- read.csv( "..//..//Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")

temp <-  aggregate(m["taxon"], m[c("year","site")], FUN = length)
temp <-  aggregate(temp["taxon"], temp[c("site")], FUN = mean)
# n= 193 

m <- m[,2:7]
m$taxon[which(m$taxon == "Porcellidium\xa0 ")] <- "Porcellidium"

m$taxon[which(m$taxon == "")] <- NA
m$taxon <- gsub( " ", ".", m$taxon )

m$taxon[which(m$taxon == "Gammaidean.Amphipod...1mm.")] <- "Gammaidean.Amphipod"
m$taxon[which(m$taxon == "Harpacticoid.copepod..multiple.species.")] <- "Harpacticoida"
m$taxon[which(m$taxon == "Lottia.alveus.paralella")] <- "Lottia.alveus"
m$taxon[which(m$taxon == "Phyllodocidae")] <- "Phyllodocida" #"Polychaeta"


# replace spaces with periods for consistency and merging names
# bring in the updated data 
mtaxa_update <- read.csv( "..//..//R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )
mtaxa_update$taxon <- gsub( " ", ".", mtaxa_update$taxon )

# merge
mTaxa <- merge( m,  mtaxa_update, by = "taxon", all.m = T )

# mTaxa <- mTaxa[complete.cases(mTaxa$taxon),]

rm <- c("eggs","Empty.polychaete.tubes", "Fluff","foil","Lacuna.spp..Eggs", 
        "larvae", "microfiber", "microplastic", "Pleocyemata", "Lacuna.spp..Eggs", "larvae", "microfiber", 
        "microplastic",  "Stalked.thing", "Stalked.thing","tin.foil", "Tiny.red.striped.snail", 
        "Unk.Fish.Eggs.", "Unk1..polychaete..", "Unk2..polychaete..", "unknown.squiggle",
        "Unknown.sea.star..I.think.only.1.sp.", "Unknown.red.Halacarid.mite", "Sea.anemone..perhaps.Epiactis.prolifera..difficult.to.ID.",
        "Sea.anemone..maybe.Diadumene.lineata.but.difficult.to.ID.", "Nematoda", "Nematode")

d <- mTaxa[!mTaxa$taxon %in% rm,] # nspp = 158
rmSp <- length(unique(subset(d, remove == "1")$taxon))
