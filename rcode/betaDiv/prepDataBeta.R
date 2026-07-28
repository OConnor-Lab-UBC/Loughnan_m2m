# aim of this file is to generate the data needed to run the m2m model:
# Started April 11, 2025 by D. Loughnan

# For the joint model need to have:
# 1. dataset of microbial data as alpha diversity metrics
# 2. dataset of macro data---diversity data, habitat data, climate data

# rm(list=ls())
# options(stringsAsFactors = FALSE)
# # 
 setwd("~/Documents/github/Loughnan_m2m")

require(vegan)
library(plyr)
library(dplyr)
library(reshape2)
library(shinystan)
library(tidyr)
library(stringr)

 setwd("~/Documents/github/Calvert_O-Connor_eelgrass/")
 m <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv" )
 
 st <- c("choked_inner","choked_sandspit","goose_south_east", "goose_south_west", "mcmullins_north","mcmullins_south","triquet_north", "triquet_south")

bc2014 <- list()

for(s in 1:length(st)){
 temp = m[m$site == st[s],]
 m.meta <- temp[1:7]
 
 # community data
 m.comm <- temp[8:ncol(temp)]
 
 # - save the bray-curtis distances
 for( i in 2014:2017 ){
    path <- "output/brayCurtis/plotYr/"
   metai <- m.meta %>%
     filter(year==i) #%>%
   # unite(sample, site, sample, sep="_")
   # filter zero columns and rows
   comm.tmp <- m.comm[m.meta$year==i,]
   keep.col <- which(colSums(comm.tmp)>0)
   keep.row <- which(rowSums(comm.tmp)>0)
   meta.use <- metai[ keep.row, ]
   comm.use <- comm.tmp[ keep.row,keep.col ]
   # write_csv( meta.use, paste0(path,i,"_macroeuk_metadata.csv") )
   sample.names <- make.cepnames(meta.use$ID)
   commdisti <- vegdist( comm.use, method = "bray" )
   commdisti <- as.matrix(commdisti)
   rownames(commdisti) <- sample.names
   colnames(commdisti) <- sample.names
   write_csv( data.frame(commdisti), paste0(path,st[s],"_",i,"_macroeuk_braycurtis_finest_27.csv") )
   
   rownames(commdisti) <- colnames(commdisti)
   temp <- matrixConvert(commdisti, colname = c("comm_id1", "comm_id2", "bray"))
   temp$year <- i
   temp$site <- st[s]
   write_csv( data.frame(temp), paste0(path,st[s],"_",i,"_betaBayes.csv") )
   
   # commdisti <- vegdist( comm.family[meta.family$year==i,], method = "bray" )
   # commdisti <- as.matrix(commdisti)
   # rownames(commdisti) <- sample.names
   # colnames(commdisti) <- sample.names
   # write_csv( data.frame(commdisti), paste0(path,i,"_macroeuk_braycurtis_family.csv") )
   # commdisti <- vegdist( comm.coarse[meta.coarse$year==i,], method = "bray" )
   # commdisti <- as.matrix(commdisti)
   # rownames(commdisti) <- sample.names
   # colnames(commdisti) <- sample.names
   # write_csv( data.frame(commdisti), paste0(path,i,"_macroeuk_braycurtis_coarse.csv") )
    }
}


#####################################################
files <- list.files(path = "output/brayCurtis/plotYr", pattern =".csv" )
files[1]

for(f in files){
  f <- 1
  d <- read.csv(paste0("output/brayCurtis/plotYr/",files[f]))
  rownames(d) <- colnames(d)
  temp <- matrixConvert(d, colname = c("comm_id1", "comm_id2", "bray"))
  temp$year <- i
  temp$site <- s
  
}
bray2014 <- bray2014[2:75]

rownames(bray2014) <- colnames(bray2014)
spec2014 <- matrixConvert(bray2014, colname = c("comm_id1", "comm_id2", "bray"))
## Micobial:
# based on code from R_Code_and_Analysis/alphadiversity/alpha_chao1_shannon_pielou.R

# prok16 <- read.csv("data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
# # names(prok16)[1:16]
# 
# ### Creating an object to store abundances only
# ### 16S prokaryotes
# abund16 <- prok16[,16:ncol(prok16)]
# # n= 165 prokaryotes
# # siteInfo <- prok16[, c("year", "region_year", "region", "site_quadrat_id", "site", "host_species", "host_type", "sample_type", "survey_type", "quadrat_id")]  
# 
# # Calculate beta diversity metrics
# 
# ## creating data frame with alpha metrics and metadata
# prok16Site <- prok16[,5:15]
# 
# ############ 18S microeukaryotes ############
# 
# micro18 <- read.csv("data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
# micro18Sites <- micro18[,2:9]
# 
# ### Creating an object to store abundances only
# abund18 <- micro18[,10:ncol(micro18)]
# 
# # Calculate beta diversity metrics
# 
# ## creating data frame with alpha metrics and metadata
# 
# ############ Inverts macroeukaryotes ############
# ### Read table metadata and abundances
# m <- read.csv( "output/cleanedGrazer.csv")
# 
# m.meta <- m[,1:5]
# 
# ### Creating an object to store abundances only
# m.comm <- m[,6:ncol(m)]
# 
# # - save the bray-curtis distances
# for( i in 2014:2017 ){
#   path <- "output/betaDiv/betaDivTemp"
#   metai <- m.meta %>%
#     filter(year==i) #%>%
#   # unite(sample, site, sample, sep="_")
#   # filter zero columns and rows
#   comm.tmp <- m.comm[m.meta$year==i,]
#   keep.col <- which(colSums(comm.tmp)>0)
#   keep.row <- which(rowSums(comm.tmp)>0)
#   meta.use <- metai[ keep.row, ]
#   comm.use <- comm.tmp[ keep.row,keep.col ]
#   write.csv( meta.use, paste0(path,i,"_macroeuk_metadata.csv") )
#   sample.names <- (meta.use$sampleID)
#   commdisti <- vegdist( comm.use, method = "bray" )
#   commdisti <- as.matrix(round(commdisti,2))
#   rownames(commdisti) <- sample.names
#   colnames(commdisti) <- sample.names
#   write.csv( data.frame(commdisti), paste0(path,i,"_macroeuk_braycurtis_finest.csv") )
# }

## creating data frame with alpha metrics and metadata

############ meadow traitvs ############
# sg <- read.csv("data/MASTER_seagrass_metrics_20200214.csv")
# sg <- sg[,2:ncol(sg)]
# 
# shoot <- read.csv("data/MASTER_shoots.csv")
# colnames(shoot)[colnames(shoot) == "id"] <- "quadrat_id"
# 
# seagrass <- merge(sg, shoot, by = c("site","quadrat_id","year"), all = T)
# 
# rmSite <- c("choked_wolf","pruth_bay", "wolf_choked_pass")  
# seagrass <- seagrass[!seagrass$site %in% rmSite,]
# seagrass$site[which(seagrass$site == "mcmullins_north")] <- "mcmullin_north"
# seagrass$site[which(seagrass$site == "mcmullins_south")] <- "mcmullin_south"
# 
# # Hakai 2014 biomass data:
# den2014 <- read.csv("data/Shoot_density.csv")
# den2014$region[which(den2014$region == "Triquet")] <- "triquet"
# den2014$region[which(den2014$region == "Goose")] <- "goose"
# den2014$region[which(den2014$region == "Choked")] <- "choked"
# den2014$region[which(den2014$region == "McMullin")] <- "mcmullin"
# 
# den2014$site[which(den2014$site == "L_Choked")] <- "choked_inner"
# den2014$site[which(den2014$site == "Sandspit")] <- "choked_sandspit"
# den2014$site[which(den2014$site == "Goose_E")] <- "goose_south_east"
# den2014$site[which(den2014$site == "Goose_W")] <- "goose_south_west"
# den2014$site[which(den2014$site == "McMullin_N")] <- "mcmullin_north"
# den2014$site[which(den2014$site == "McMullin_S")] <- "mcmullin_south"
# den2014$site[which(den2014$site == "Triquet_B")] <- "triquet_south"
# den2014$site[which(den2014$site == "N_Triquet")] <- "triquet_north"
# den2014$year <- "2014"
# colnames(den2014)[colnames(den2014) == "density"] <- "quadrat_shoot_density"
# colnames(den2014)[colnames(den2014) == "site_no"] <- "quadrat_id"
# 
# sd2014 <- den2014[,c("site","region","year","quadrat_id", "sample","depth_cor","quadrat_shoot_density")]
# 
# #c("quadrat_shoot_density", "depth_cor")
# # bio1 <- read.csv("m2m/data/Hakai2014/csvFiles/bio.tot.plot2.csv") # this data all appears duplicated in the Macrophyte_biomass_MHL.csv
# bio2014 <- read.csv("data/macro_plot_biomass_MHL.csv")
# bio2014$region[which(bio2014$region == "Triquet")] <- "triquet"
# bio2014$region[which(bio2014$region == "Goose")] <- "goose"
# bio2014$region[which(bio2014$region == "Choked")] <- "choked"
# bio2014$region[which(bio2014$region == "McMullin")] <- "mcmullin"
# 
# bio2014$site <- bio2014$site_id
# bio2014$site[which(bio2014$site == "CL")] <- "choked_inner"
# bio2014$site[which(bio2014$site == "CS")] <- "choked_sandspit"
# bio2014$site[which(bio2014$site == "GE")] <- "goose_south_east"
# bio2014$site[which(bio2014$site == "GW")] <- "goose_south_west"
# bio2014$site[which(bio2014$site == "MN")] <- "mcmullin_north"
# bio2014$site[which(bio2014$site == "MS")] <- "mcmullin_south"
# bio2014$site[which(bio2014$site == "TB")] <- "triquet_south"
# bio2014$site[which(bio2014$site == "TN")] <- "triquet_north"
# 
# bio2014$sample_id[which(bio2014$sample_id == "Macroalage")] <- "Macroalgae"
# bio2014$sample_id[which(bio2014$sample_id == "Maccroalgae")] <- "Macroalgae"
# 
# # bio2014$year <- "2014"
# # bio2014$quadrat_biomass_g <- bio2014$tot.bio*1000 # I think this is the same as quadrat_biomass_g in seagrass dataset
# bio2014$quadrat_biomass_g <- bio2014$final_dry_wgt_kg*1000 
# # hist(bio2014$final_dry_wgt_kg*1000, col = rgb(72 / 255, 38 / 255, 119 / 255, alpha = 0.5), xlim = c(0,60), ylim = c(0,100)); hist(seagrass$quadrat_biomass_g, col = rgb(149 / 255, 216 / 255, 64 / 255, alpha = 0.5),add = T)
# 
# # colnames(bio2014)[colnames(bio2014) == "tot.bio"] <- "shoot.biomass.g"
# # colnames(bio2014)[colnames(bio2014) == "tot.bio.msq"] <- "quadrat.shoot.biomass"
# # colnames(bio2014)[colnames(bio2014) == "sample_no"] <- "quadrat_id"
# bio2014Sub <- bio2014[,c("site","year","date","sample_no","mult_sample_no", "sample_id","quadrat_biomass_g")]
# 
# sg2014 <- subset(bio2014Sub, sample_id == "Seagrass")
# sg2014$mult_sample_no <- paste(sg2014$sample_no,sg2014$mult_sample_no, sep = "")
# 
# mAlgae2014 <- subset(bio2014Sub, sample_id == "Macroalgae"); colnames(mAlgae2014)[colnames(mAlgae2014) == "quadrat_biomass_g"] <- "quadrat_macroalgae_g"
# mAlgae2014 <- mAlgae2014[,c("site","year","sample_no","quadrat_macroalgae_g")]
# 
# bio2014Combo <- merge(sg2014, mAlgae2014, by = c("site", "year", "sample_no"), all = T)
# bio2014Combo$macroAm2 <- bio2014Combo$quadrat_macroalgae_g*16
# bio2014Combo <- bio2014Combo[,c("site","year","sample_no","mult_sample_no","quadrat_biomass_g","quadrat_macroalgae_g","macroAm2")]
# 
# lai2014 <- read.csv("data/Leaf_area.csv")
# 
# lai2014$region <- lai2014$site
# lai2014$region[which(lai2014$region == "N_Triquet")] <- "triquet"
# lai2014$region[which(lai2014$region == "Triquet_B")] <- "triquet"
# lai2014$region[which(lai2014$region == "Goose_E")] <- "goose"
# lai2014$region[which(lai2014$region == "Goose_W")] <- "goose"
# lai2014$region[which(lai2014$region == "L_Choked")] <- "choked"
# lai2014$region[which(lai2014$region == "Sandspit")] <- "choked"
# lai2014$region[which(lai2014$region == "McMullin_S")] <- "mcmullin"
# lai2014$region[which(lai2014$region == "McMullin_N")] <- "mcmullin"
# 
# 
# lai2014$site <- lai2014$site_id
# lai2014$site[which(lai2014$site == "LC")] <- "choked_inner"
# lai2014$site[which(lai2014$site == "SC")] <- "choked_sandspit"
# lai2014$site[which(lai2014$site == "GE")] <- "goose_south_east"
# lai2014$site[which(lai2014$site == "GW")] <- "goose_south_west"
# lai2014$site[which(lai2014$site == "MN")] <- "mcmullin_north"
# lai2014$site[which(lai2014$site == "MS")] <- "mcmullin_south"
# lai2014$site[which(lai2014$site == "TB")] <- "triquet_south"
# lai2014$site[which(lai2014$site == "NT")] <- "triquet_north"
# 
# lai2014$year <- "2014"
# 
# lai2014 <- lai2014[,c("site", "year","sample_no","leaf_area_index")]
# lai2014 <- lai2014[complete.cases(lai2014$leaf_area_index),]
# 
# bio2014Final <- merge(bio2014Combo, lai2014, by = c("site", "year", "sample_no"), all = T )
# 
# bio2014Fin <- aggregate(bio2014Final[c("quadrat_biomass_g","quadrat_macroalgae_g","leaf_area_index")], bio2014Final[c("year","site","sample_no")], FUN = mean)
# colnames(bio2014Fin)[colnames(bio2014Fin) == "sample_no"] <- "quadrat_id"
# colnames(bio2014Fin)[colnames(bio2014Fin) == "leaf_area_index"] <- "quadrat_lai"
# 
# 
# seagrassF <- rbind.fill(seagrass, bio2014Fin)
# 
# # shoot density, macroalgae, lai, microepiphyte, shoot length, shoot width, shoot biomass
# abio <- read.csv("data/MASTER_abiotic_20200214.csv") # so little data not worth including further---depth only for 2016-2018, bed area only for 
# 
# # get bed area 
# bedA <- subset(abio, !is.na(bed_area_m2))
# bedA <- bedA[,c("site","bed_area_m2")]
# 
# bedA$site[which(bedA$site == "mcmullins_north")] <- "mcmullin_north"
# bedA$site[which(bedA$site == "mcmullins_south")] <- "mcmullin_south"
# 
# # get depth data
# depthActual <- read.csv("data/Actual_depth_meso.csv")
# 
# depthActual$site <- depthActual$site_id
# depthActual$site[which(depthActual$site == "LC")] <- "choked_inner"
# depthActual$site[which(depthActual$site == "SC")] <- "choked_sandspit"
# depthActual$site[which(depthActual$site == "GE")] <- "goose_south_east"
# depthActual$site[which(depthActual$site == "GW")] <- "goose_south_west"
# depthActual$site[which(depthActual$site == "MN")] <- "mcmullin_north"
# depthActual$site[which(depthActual$site == "MS")] <- "mcmullin_south"
# depthActual$site[which(depthActual$site == "TB")] <- "triquet_south"
# depthActual$site[which(depthActual$site == "NT")] <- "triquet_north"
# 
# colnames(depthActual)[colnames(depthActual) == "sample_no"] <- "quadrat_id"
# depthActual <- depthActual[,c("site","quadrat_id","actual_depth")]
# depthSum <- aggregate(depthActual[c("actual_depth")], depthActual[c("site")], FUN = mean)


# depthNS <- read.csv("m2m/data/Hakai2014/NearshoreMasterSiteList_NearshoreSites.csv")
# depthSub <- depthNS[,c("Site.Name","Depth_m")]
# 
# siteSG <- unique(seagrassF$site)
# 
# depthSub$Site.Name[which(depthSub$Site.Name == "Choked Interior ")] <- "choked_inner"
# depthSub$Site.Name[which(depthSub$Site.Name == "Choked Interior")] <- "choked_inner"
# 
# depthSub$Site.Name[which(depthSub$Site.Name == "Choked Sandspit")] <- "choked_sandspit"
# depthSub$Site.Name[which(depthSub$Site.Name == "Goose SE")] <- "goose_south_east"
# depthSub$Site.Name[which(depthSub$Site.Name == "Goose SE ")] <- "goose_south_east"
# 
# depthSub$Site.Name[which(depthSub$Site.Name == "Goose SW ")] <- "goose_south_west"
# depthSub$Site.Name[which(depthSub$Site.Name == "Goose SW")] <- "goose_south_west"
# 
# depthSub$Site.Name[which(depthSub$Site.Name == "McMullin North")] <- "mcmullin_north"
# depthSub$Site.Name[which(depthSub$Site.Name == "McMullin North ")] <- "mcmullin_north"
# 
# depthSub$Site.Name[which(depthSub$Site.Name == "McMullin South")] <- "mcmullin_south"
# depthSub$Site.Name[which(depthSub$Site.Name == "McMullin South ")] <- "mcmullin_south"
# 
# depthSub$Site.Name[which(depthSub$Site.Name == "Triquet Bay")] <- "triquet_south"
# depthSub$Site.Name[which(depthSub$Site.Name == "Triquet Bay ")] <- "triquet_south"
# 
# depthSub$Site.Name[which(depthSub$Site.Name == "Triquet North")] <- "triquet_north"
# depthSub$Site.Name[which(depthSub$Site.Name == "Triquet North ")] <- "triquet_north"
# 
# depthSite <- depthSub[depthSub$Site.Name %in% siteSG,]
# depth2  <- subset(depthSite, Depth_m != "" )
# depth2  <- subset(depth2, Depth_m != "profile" )
# 
# colnames(depth2)[colnames(depth2) == "Site.Name"] <- "site"
#write.csv(depth2, "m2m/output/nearshoreDepths.csv", row.names = F)

# seagrassF <- merge(seagrassF, bedA, by = c("site"), all.x = T)
#seagrassF <- merge(seagrassF, depth2, by = c("site"), all.x = T)
# seagrassF <- merge(seagrassF, depthSum, by = c("site"), all.x = T)

# write.csv(bedA, "m2m/output/bedArea8site1yr.csv", row.names = F)
# alphaGraz$region[which(alphaGraz$region == "mcmullins")] <- "mcmullin"
# 
# alphaMacro <- merge(alphaGraz, seagrassF, by = c("site", "quadrat_id","year"), all = T)
# 
# alphaMacro <- alphaMacro[,c("site","region","quadrat_id","year","richness","shannon", "pielou","chao1", "quadrat_macroalgae_g", "quadrat_lai","quadrat_biomass_g")]
#  # aTemp <- alphaMacro[,c("shannon", "pielou", "chao1",  "quadrat_biomass_g")]
#  # temp <- seagrass[complete.cases(seagrass$quadrat_macroalgae_g),]
# 
# ############ climate variables ############
# 
# # the ideal climate data is the monthly average for the summer months using the daily temperature data
# # but if this is not available then using the most populated ctd data again averaged for that month
# 
# seasonal <- read.csv("input/seasonalMatTemp.csv")
# seasonal <- subset(seasonal, site!="PruthBay02")
# seasonal$site[which(seasonal$site == "chokedinner")] <- "choked_inner"
# seasonal$site[which(seasonal$site == "chokedsandspit")] <- "choked_sandspit"
# seasonal$site[which(seasonal$site == "gooseSE")] <- "goose_south_east"
# seasonal$site[which(seasonal$site == "gooseSW")] <- "goose_south_west"
# seasonal$site[which(seasonal$site == "McMullinN")] <- "mcmullin_north"
# seasonal$site[which(seasonal$site == "McMullinS")] <- "mcmullin_south"
# seasonal$site[which(seasonal$site == "PruthPocket")] <- "pruth_pocket"
# seasonal$site[which(seasonal$site == "TriquetBay02")] <- "triquet_south"
# seasonal$site[which(seasonal$site == "TriquetNSeagraa")] <- "triquet_north"
# 
# ctd <- read.csv("input/ctdSeasonalStn.csv")
# ctd$region[which(ctd$region == "mcmullin")] <- "mcmullin"
# 
# # mean temps that I use to fill in the gaps
# macroClim <- read.csv("input/macroSeasonal.csv")
# macroClim$region[which(macroClim$region == "mcmullin")] <- "mcmullin"
# 
# ctd$iter <- 1:nrow(ctd)
# ctd$year <- as.character(as.factor(ctd$year))
# ctd$sumTemp[which(is.na(ctd$sumTemp) & ctd$region == "goose" & ctd$year == "2015")] <- 13.25923
# ctd$sumTemp[which(is.na(ctd$sumTemp) & ctd$region == "mcmullin" & ctd$year == "2015")] <- 13.07650
# ctd$sumTemp[which(is.na(ctd$sumTemp) & ctd$region == "triquet" & ctd$year == "2015")] <- 13.02320
# ctd$sumTemp[which(is.na(ctd$sumTemp) & ctd$region == "pruth" & ctd$year == "2015")] <- 12.80765
# ctd$year <- as.numeric(ctd$year)
# 
# # # Now get a complete row of summer temp data:
# # ctd$comboTemp <- ctd$summerTemp
# # ctd$comboTemp <- ifelse(is.na(ctd$comboTemp), ctd$sumTemp, ctd$comboTemp)
# # ctd$tempType <- ctd$summerTemp
# # ctd$tempType <- ifelse(is.na(ctd$tempType), "ctd", ctd$comboTemp)
# 
# # #if(i in 1:nrow(alphaClim)) == NA then alphaClim$summerTemp[i, ] == alphaClim$sumTemp[i,]
# ctd$region[which(ctd$site == "triquet_south")] <- "triquet"
# ctd$region[which(ctd$site == "triquet_north")] <- "triquet"
# ctd$region[which(ctd$site == "goose_south_west")] <- "goose"
# ctd$region[which(ctd$site == "goose_south_east")] <- "goose"
# ctd$region[which(ctd$site == "pruth_pocket")] <- "pruth"
# ctd$region[which(ctd$site == "choked_sandspit")] <- "choked"
# ctd$region[which(ctd$site == "choked_inner")] <- "choked"
# ctd$region[which(ctd$site == "mcmullin_south")] <- "mcmullin"
# ctd$region[which(ctd$site == "mcmullin_north")] <- "mcmullin"
# 
# 
# head(lai2014)
# ctd <- merge(ctd, bedA, by = c("site"), all.x = T)
# ctd <- merge(alphaClim, depthSum, by = c("site"), all.x = T)

# # fix the missing depths: none for Goose north, but nearshore has one for pruth pocket:
# alphaClim$actual_depth[which(alphaClim$site == "pruth_pocket")] <- "-5"
# 
# 
# sort(unique(alphaClim$site))
# 
# 
