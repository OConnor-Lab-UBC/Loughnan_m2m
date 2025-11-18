# aim of this file is to generate the data needed to run the m2m model:
# Started April 11, 2025 by D. Loughnan

# For the joint model need to have:
# 1. dataset of microbial data as alpha diversity metrics
# 2. dataset of macro data---diversity data, habitat data, climate data

# rm(list=ls())
# options(stringsAsFactors = FALSE)

require(vegan)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)

## Microbial:
# based on code from R_Code_and_Analysis/alphadiversity/alpha_chao1_shannon_pielou.R

prok16 <- read.csv("data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
# names(prok16)[1:16]

### Creating an object to store abundances only
### 16S prokaryotes
abund16 <- prok16[,16:ncol(prok16)]
# n= 165 prokaryotes
# siteInfo <- prok16[, c("year", "region_year", "region", "site_quadrat_id", "site", "host_species", "host_type", "sample_type", "survey_type", "quadrat_id")]  

# Calculate alpha diversity metrics
shannon <- diversity(abund16, index = "shannon")
richness <- specnumber(abund16) 
pielou <- shannon/log(richness)
chao1 <- estimateR(abund16)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
prok16Site <- prok16[,5:15]
alpha16 <- data.frame(prok16Site,richness, shannon, pielou,chao1)

############ 18S microeukaryotes ############

micro18 <- read.csv("data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
micro18Sites <- micro18[,2:9]

### Creating an object to store abundances only
abund18 <- micro18[,10:ncol(micro18)]

#n = 141 
# Calculate alpha diversity metrics
shannon <- diversity(abund18, index = "shannon")
richness <- specnumber(abund18) 
pielou <- shannon/log(richness)
chao1 <- estimateR(abund18)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alpha18 <- data.frame(micro18Sites,richness, shannon, pielou,chao1)

############ Inverts macroeukaryotes ############
### Read table metadata and abundances
m <- read.csv( "data/MASTER_grazers.csv")

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
mtaxa_update <- read.csv( "data/O'Connor_hakai_seagrass_taxa_edit.csv" )
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
d <- subset(d, is.na(remove)) # nspp = 141

str(d)
# cleaning spp names and lumping:
d$taxon <- str_remove(d$taxon, ".sp.+$") # nspp = 139

# lumping: 
d$taxonL <- d$taxon
d$taxonL[which(d$taxonL == "Amphissa.columbiana")] <- "Amphissa"
d$taxonL[which(d$taxonL == "Ampithoe.lacertosa")] <- "Ampithoidae"#
d$taxonL[which(d$taxonL == "Ampithoe.sectimanus")] <- "Ampithoidae"
d$taxonL[which(d$taxonL == "Ampithoe")] <- "Ampithoidae"
d$taxonL[which(d$taxonL == "Aoroides.intermedia")] <- "Aoroides"
# d$taxonL[which(d$taxonL == "Caprella.californica")] 
# d$taxonL[which(d$taxonL == "Caprella.laeviuscula")] <- 
# d$taxonL[which(d$taxonL == "Caprella.mutica")] <-
d$taxonL[which(d$taxonL == "Ericthonius.rubricornis")] <- "Ericthonius"
d$taxonL[which(d$taxonL == "Gammaropsis.barnardi")] <- "Gammaropsis"
d$taxonL[which(d$taxonL == "Gammaropsis.thomsoni")] <- "Gammaropsis"
d$taxonL[which(d$taxonL == "Hesionura.coineaui.difficilis")] <- "Hesionidae"
d$taxonL[which(d$taxonL == "Idotea")] <- "Idoteidae" # how many
d$taxonL[which(d$taxonL == "Idotea.montereynensis")] <- "Idoteidae"
d$taxonL[which(d$taxonL == "Idotea.resecata")] <- "Idoteidae"
d$taxonL[which(d$taxonL == "Ischyrocerus")] <- "Ischyroceridae"
d$taxonL[which(d$taxonL == "Ischyrocerus.anguipes")] <- "Ischyroceridae"
# d$taxonL[which(d$taxonL == "Lacuna.vincta")] <- 
# d$taxonL[which(d$taxonL == "Lacuna.variegata")] <-
d$taxonL[which(d$taxonL == "Leptochelia.dubia")] <- "Leptochelia"
d$taxonL[which(d$taxonL == "Ischyrocerus.anguipes")] <- "Ischyroceridae"
d$taxonL[which(d$taxonL == "Lirularia.lirulata")] <- "Lirularia"
d$taxonL[which(d$taxonL == "Lottia.alveus")] <- "Lottidae"
d$taxonL[which(d$taxonL == "Monocorophium.acherusicum")] <- "Monocorophium"
d$taxonL[which(d$taxonL == "Mytilus.trossulus")] <- "Mytilidae"
d$taxonL[which(d$taxonL == "Pagurus")] <- "Paguridae"
d$taxonL[which(d$taxonL == "Pentidotea.wosnesenskii")] <- "Pentidotea"
d$taxonL[which(d$taxonL == "Pholoe.tuberculata")] <- "Pholoidae"
d$taxonL[which(d$taxonL == "Photis")] <- "Photis.brevipes"
d$taxonL[which(d$taxonL == "Phyllodoce")] <- "Phyllodocida"
d$taxonL[which(d$taxonL == "Platynereis.bicanaliculata")] <- "Platynereis"
d$taxonL[which(d$taxonL == "Pontogeneia.inermis")] <- "Pontogeneia"
d$taxonL[which(d$taxonL == "Pontogeneia.intermedia")] <- "Pontogeneia"
d$taxonL[which(d$taxonL == "Pontogeneia.rostrata")] <- "Pontogeneia"
d$taxonL[which(d$taxonL == "Protohyale.intermedia")] <- "Protohyale"
d$taxonL[which(d$taxonL == "Pugettia.richii")] <- "Pugettia"
d$taxonL[which(d$taxonL == "Pycnogonida")] <- "Pycnogonid"
d$taxonL[which(d$taxonL == "Thorlaksonius.subcarinatus")] <- "Thorlaksonius"

length(unique(d$taxonL)) # 107
# define regions as first word of site
d$region <- unlist(lapply( strsplit(d$site,split = "_"), function(z) z[1]))
# unique sample ID to differential samples from different sites
d$sampleID <- with(d, paste(site,sample,sep = "_"))

# change year to factor

# summarize taxon counts per sample
abund <- d %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, region, site, sample, sampleID, taxon5,taxonL ) %>% 
  summarize( abundance=length(size))

 write.csv(abund, "output/cleanedGrazer.csv", row.names = F)
# # summarinze mean abundance per site
# m.mean_finest <- m.sum_finest %>% 
#   # unite( "ID", year,site,sample, remove=FALSE ) %>% 
#   group_by( year, year0, yearF, region, site,  taxon  ) %>% 
#   summarize( abundance=mean(abundance)) 

# make a community dataset
abundG <- abund %>% 
  spread(taxonL, abundance, fill=0 ) 

abundGSite <- abundG[,1:5]

### Creating an object to store abundances only
abundMatrix <- abundG[,8:ncol(abundG)]

# Calculate alpha diversity metrics
shannon <- diversity(abundMatrix, index = "shannon")
richness <- specnumber(abundMatrix) 
pielou <- shannon/log(richness)
chao1 <- estimateR(abundMatrix)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alphaGraz <- data.frame(abundGSite,richness,shannon, pielou,chao1)

alphaGraz$quadrat_id <- alphaGraz$sample
alphaGraz$quadrat_id <- gsub('\\D+','', alphaGraz$quadrat_id)


alphaGraz$site[which(alphaGraz$site == "mcmullins_north")] <- "mcmullin_north"
alphaGraz$site[which(alphaGraz$site == "mcmullins_south")] <- "mcmullin_south"

############ meadow traitvs ############
sg <- read.csv("data/MASTER_seagrass_metrics_20200214.csv")
sg <- sg[,2:ncol(sg)]

shoot <- read.csv("data/MASTER_shoots.csv")
colnames(shoot)[colnames(shoot) == "id"] <- "quadrat_id"

seagrass <- merge(sg, shoot, by = c("site","quadrat_id","year"), all = T)

rmSite <- c("choked_wolf","pruth_bay", "wolf_choked_pass")  
seagrass <- seagrass[!seagrass$site %in% rmSite,]
seagrass$site[which(seagrass$site == "mcmullins_north")] <- "mcmullin_north"
seagrass$site[which(seagrass$site == "mcmullins_south")] <- "mcmullin_south"

# Hakai 2014 biomass data:
den2014 <- read.csv("data/Shoot_density.csv")
den2014$region[which(den2014$region == "Triquet")] <- "triquet"
den2014$region[which(den2014$region == "Goose")] <- "goose"
den2014$region[which(den2014$region == "Choked")] <- "choked"
den2014$region[which(den2014$region == "McMullin")] <- "mcmullin"

den2014$site[which(den2014$site == "L_Choked")] <- "choked_inner"
den2014$site[which(den2014$site == "Sandspit")] <- "choked_sandspit"
den2014$site[which(den2014$site == "Goose_E")] <- "goose_south_east"
den2014$site[which(den2014$site == "Goose_W")] <- "goose_south_west"
den2014$site[which(den2014$site == "McMullin_N")] <- "mcmullin_north"
den2014$site[which(den2014$site == "McMullin_S")] <- "mcmullin_south"
den2014$site[which(den2014$site == "Triquet_B")] <- "triquet_south"
den2014$site[which(den2014$site == "N_Triquet")] <- "triquet_north"
den2014$year <- "2014"
colnames(den2014)[colnames(den2014) == "density"] <- "quadrat_shoot_density"
colnames(den2014)[colnames(den2014) == "site_no"] <- "quadrat_id"

sd2014 <- den2014[,c("site","region","year","quadrat_id", "sample","depth_cor","quadrat_shoot_density")]

#c("quadrat_shoot_density", "depth_cor")
# bio1 <- read.csv("m2m/data/Hakai2014/csvFiles/bio.tot.plot2.csv") # this data all appears duplicated in the Macrophyte_biomass_MHL.csv
bio2014 <- read.csv("data/macro_plot_biomass_MHL.csv")
bio2014$region[which(bio2014$region == "Triquet")] <- "triquet"
bio2014$region[which(bio2014$region == "Goose")] <- "goose"
bio2014$region[which(bio2014$region == "Choked")] <- "choked"
bio2014$region[which(bio2014$region == "McMullin")] <- "mcmullin"

bio2014$site <- bio2014$site_id
bio2014$site[which(bio2014$site == "CL")] <- "choked_inner"
bio2014$site[which(bio2014$site == "CS")] <- "choked_sandspit"
bio2014$site[which(bio2014$site == "GE")] <- "goose_south_east"
bio2014$site[which(bio2014$site == "GW")] <- "goose_south_west"
bio2014$site[which(bio2014$site == "MN")] <- "mcmullin_north"
bio2014$site[which(bio2014$site == "MS")] <- "mcmullin_south"
bio2014$site[which(bio2014$site == "TB")] <- "triquet_south"
bio2014$site[which(bio2014$site == "TN")] <- "triquet_north"

bio2014$sample_id[which(bio2014$sample_id == "Macroalage")] <- "Macroalgae"
bio2014$sample_id[which(bio2014$sample_id == "Maccroalgae")] <- "Macroalgae"

# bio2014$year <- "2014"
# bio2014$quadrat_biomass_g <- bio2014$tot.bio*1000 # I think this is the same as quadrat_biomass_g in seagrass dataset
bio2014$quadrat_biomass_g <- bio2014$final_dry_wgt_kg*1000 
# hist(bio2014$final_dry_wgt_kg*1000, col = rgb(72 / 255, 38 / 255, 119 / 255, alpha = 0.5), xlim = c(0,60), ylim = c(0,100)); hist(seagrass$quadrat_biomass_g, col = rgb(149 / 255, 216 / 255, 64 / 255, alpha = 0.5),add = T)

# colnames(bio2014)[colnames(bio2014) == "tot.bio"] <- "shoot.biomass.g"
# colnames(bio2014)[colnames(bio2014) == "tot.bio.msq"] <- "quadrat.shoot.biomass"
# colnames(bio2014)[colnames(bio2014) == "sample_no"] <- "quadrat_id"
bio2014Sub <- bio2014[,c("site","year","date","sample_no","mult_sample_no", "sample_id","quadrat_biomass_g")]

sg2014 <- subset(bio2014Sub, sample_id == "Seagrass")
sg2014$mult_sample_no <- paste(sg2014$sample_no,sg2014$mult_sample_no, sep = "")

mAlgae2014 <- subset(bio2014Sub, sample_id == "Macroalgae"); colnames(mAlgae2014)[colnames(mAlgae2014) == "quadrat_biomass_g"] <- "quadrat_macroalgae_g"
mAlgae2014 <- mAlgae2014[,c("site","year","sample_no","quadrat_macroalgae_g")]

bio2014Combo <- merge(sg2014, mAlgae2014, by = c("site", "year", "sample_no"), all = T)
bio2014Combo$macroAm2 <- bio2014Combo$quadrat_macroalgae_g*16
bio2014Combo <- bio2014Combo[,c("site","year","sample_no","mult_sample_no","quadrat_biomass_g","quadrat_macroalgae_g","macroAm2")]

lai2014 <- read.csv("data/Leaf_area.csv")

lai2014$region <- lai2014$site
lai2014$region[which(lai2014$region == "N_Triquet")] <- "triquet"
lai2014$region[which(lai2014$region == "Triquet_B")] <- "triquet"
lai2014$region[which(lai2014$region == "Goose_E")] <- "goose"
lai2014$region[which(lai2014$region == "Goose_W")] <- "goose"
lai2014$region[which(lai2014$region == "L_Choked")] <- "choked"
lai2014$region[which(lai2014$region == "Sandspit")] <- "choked"
lai2014$region[which(lai2014$region == "McMullin_S")] <- "mcmullin"
lai2014$region[which(lai2014$region == "McMullin_N")] <- "mcmullin"


lai2014$site <- lai2014$site_id
lai2014$site[which(lai2014$site == "LC")] <- "choked_inner"
lai2014$site[which(lai2014$site == "SC")] <- "choked_sandspit"
lai2014$site[which(lai2014$site == "GE")] <- "goose_south_east"
lai2014$site[which(lai2014$site == "GW")] <- "goose_south_west"
lai2014$site[which(lai2014$site == "MN")] <- "mcmullin_north"
lai2014$site[which(lai2014$site == "MS")] <- "mcmullin_south"
lai2014$site[which(lai2014$site == "TB")] <- "triquet_south"
lai2014$site[which(lai2014$site == "NT")] <- "triquet_north"

lai2014$year <- "2014"

lai2014 <- lai2014[,c("site", "year","sample_no","leaf_area_index")]

bio2014Final <- merge(bio2014Combo, lai2014, by = c("site", "year", "sample_no"), all = T )

bio2014Fin <- aggregate(bio2014Final[c("quadrat_biomass_g","quadrat_macroalgae_g","leaf_area_index")], bio2014Final[c("year","site","sample_no")], FUN = mean)
colnames(bio2014Fin)[colnames(bio2014Fin) == "sample_no"] <- "quadrat_id"
colnames(bio2014Fin)[colnames(bio2014Fin) == "leaf_area_index"] <- "quadrat_lai"


seagrassF <- rbind.fill(seagrass, bio2014Fin)

# shoot density, macroalgae, lai, microepiphyte, shoot length, shoot width, shoot biomass
abio <- read.csv("data/MASTER_abiotic_20200214.csv") # so little data not worth including further---depth only for 2016-2018, bed area only for 

# get bed area 
bedA <- subset(abio, !is.na(bed_area_m2))
bedA <- bedA[,c("site","bed_area_m2")]

bedA$site[which(bedA$site == "mcmullins_north")] <- "mcmullin_north"
bedA$site[which(bedA$site == "mcmullins_south")] <- "mcmullin_south"

# get depth data
depthActual <- read.csv("data/Actual_depth_meso.csv")

depthActual$site <- depthActual$site_id
depthActual$site[which(depthActual$site == "LC")] <- "choked_inner"
depthActual$site[which(depthActual$site == "SC")] <- "choked_sandspit"
depthActual$site[which(depthActual$site == "GE")] <- "goose_south_east"
depthActual$site[which(depthActual$site == "GW")] <- "goose_south_west"
depthActual$site[which(depthActual$site == "MN")] <- "mcmullin_north"
depthActual$site[which(depthActual$site == "MS")] <- "mcmullin_south"
depthActual$site[which(depthActual$site == "TB")] <- "triquet_south"
depthActual$site[which(depthActual$site == "NT")] <- "triquet_north"

colnames(depthActual)[colnames(depthActual) == "sample_no"] <- "quadrat_id"
depthActual <- depthActual[,c("site","quadrat_id","actual_depth")]
depthSum <- aggregate(depthActual[c("actual_depth")], depthActual[c("site")], FUN = mean)


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
alphaGraz$region[which(alphaGraz$region == "mcmullins")] <- "mcmullin"

alphaMacro <- merge(alphaGraz, seagrassF, by = c("site", "quadrat_id","year"), all = T)

alphaMacro <- alphaMacro[,c("site","region","quadrat_id","year","richness","shannon", "pielou","chao1", "quadrat_macroalgae_g", "quadrat_lai","quadrat_biomass_g")]
 # aTemp <- alphaMacro[,c("shannon", "pielou", "chao1",  "quadrat_biomass_g")]
 # temp <- seagrass[complete.cases(seagrass$quadrat_macroalgae_g),]

############ climate variables ############

# the ideal climate data is the monthly average for the summer months using the daily temperature data
# but if this is not available then using the most populated ctd data again averaged for that month

seasonal <- read.csv("input/seasonalMatTemp.csv")
seasonal <- subset(seasonal, site!="PruthBay02")
seasonal$site[which(seasonal$site == "chokedinner")] <- "choked_inner"
seasonal$site[which(seasonal$site == "chokedsandspit")] <- "choked_sandspit"
seasonal$site[which(seasonal$site == "gooseSE")] <- "goose_south_east"
seasonal$site[which(seasonal$site == "gooseSW")] <- "goose_south_west"
seasonal$site[which(seasonal$site == "McMullinN")] <- "mcmullin_north"
seasonal$site[which(seasonal$site == "McMullinS")] <- "mcmullin_south"
seasonal$site[which(seasonal$site == "PruthPocket")] <- "pruth_pocket"
seasonal$site[which(seasonal$site == "TriquetBay02")] <- "triquet_south"
seasonal$site[which(seasonal$site == "TriquetNSeagraa")] <- "triquet_north"

alphaClim <- merge(alphaMacro, seasonal, by = c("site","year"), all = T)

ctd <- read.csv("input/ctdSeasonalStn.csv")
ctd$region[which(ctd$region == "mcmullin")] <- "mcmullin"

macroClim <- read.csv("input/macroSeasonal.csv")
macroClim$region[which(macroClim$region == "mcmullin")] <- "mcmullin"

alphaClim <- merge(alphaClim, ctd, by = c("region", "year"), all = T)
alphaClim$iter <- 1:nrow(alphaClim)
alphaClim$year <- as.character(as.factor(alphaClim$year))
alphaClim$sumTemp[which(is.na(alphaClim$sumTemp) & alphaClim$region == "goose" & alphaClim$year == "2015")] <- 13.25923
alphaClim$sumTemp[which(is.na(alphaClim$sumTemp) & alphaClim$region == "mcmullin" & alphaClim$year == "2015")] <- 13.07650
alphaClim$sumTemp[which(is.na(alphaClim$sumTemp) & alphaClim$region == "triquet" & alphaClim$year == "2015")] <- 13.02320
alphaClim$sumTemp[which(is.na(alphaClim$sumTemp) & alphaClim$region == "pruth" & alphaClim$year == "2015")] <- 12.80765
alphaClim$year <- as.numeric(alphaClim$year)

# Now get a complete row of summer temp data:
alphaClim$comboTemp <- alphaClim$summerTemp
alphaClim$comboTemp <- ifelse(is.na(alphaClim$comboTemp), alphaClim$sumTemp, alphaClim$comboTemp)
alphaClim$tempType <- alphaClim$summerTemp
alphaClim$tempType <- ifelse(is.na(alphaClim$tempType), "ctd", alphaClim$comboTemp)

#if(i in 1:nrow(alphaClim)) == NA then alphaClim$summerTemp[i, ] == alphaClim$sumTemp[i,]
alphaClim$region[which(alphaClim$site == "triquet_south")] <- "triquet"
alphaClim$region[which(alphaClim$site == "triquet_north")] <- "triquet"
alphaClim$region[which(alphaClim$site == "goose_south_west")] <- "goose"
alphaClim$region[which(alphaClim$site == "goose_south_east")] <- "goose"
alphaClim$region[which(alphaClim$site == "pruth_pocket")] <- "pruth"
alphaClim$region[which(alphaClim$site == "choked_sandspit")] <- "choked"
alphaClim$region[which(alphaClim$site == "choked_inner")] <- "choked"
alphaClim$region[which(alphaClim$site == "mcmullin_south")] <- "mcmullin"
alphaClim$region[which(alphaClim$site == "mcmullin_north")] <- "mcmullin"

alphaClim <- merge(alphaClim, bedA, by = c("site"), all.x = T)
alphaClim <- merge(alphaClim, depthSum, by = c("site"), all.x = T)

# fix the missing depths: none for Goose north, but nearshore has one for pruth pocket:
alphaClim$actual_depth[which(alphaClim$site == "pruth_pocket")] <- "-5"


sort(unique(alphaClim$site))

alphaSub <- alphaClim[,c("region", "year", "site", "quadrat_id","richness" ,"shannon","chao1","quadrat_macroalgae_g", "quadrat_lai", "quadrat_biomass_g","comboTemp","actual_depth","bed_area_m2")] 
# n = 86
alphaSix <- alphaSub[complete.cases(alphaSub),] # 101

alphaSub <- alphaClim[,c("region", "year", "site", "quadrat_id","richness" ,"shannon","chao1","quadrat_macroalgae_g", "quadrat_lai","comboTemp","actual_depth")]#,"bed_area_m2")] 
# n = 115
alphaFive <- alphaSub[complete.cases(alphaSub),]; dim(alphaFive)

alphaSub <- alphaClim[,c("region", "year", "site", "quadrat_id","richness" ,"shannon","chao1", "quadrat_biomass_g","comboTemp","actual_depth","bed_area_m2")] 
# n = 155
alphaFour <- alphaSub[complete.cases(alphaSub),]

alphaSub <- alphaClim[,c("region", "year", "site", "quadrat_id","richness" ,"shannon","chao1", "actual_depth","bed_area_m2")] # n = 174
alphatwo <- alphaSub[complete.cases(alphaSub),]

alphaSub <- alphaClim[,c("region", "year", "site", "quadrat_id","richness" ,"shannon","chao1","quadrat_lai", "quadrat_biomass_g","comboTemp","actual_depth","bed_area_m2")] 
# n = 115
alphaLAI <- alphaSub[complete.cases(alphaSub),]

 # write.csv(alphaSix, "m2m/output/alphaDivClim6Param.csv", row.names = F)
# write.csv(alphaMacro, "m2m/output/alphaDivClimMacro.csv", row.names = F)


### prep data for shoot density ananlysis:

# shootDen <- alphaClim[,c("region", "site", "year", "quadrat_id", "quadrat_shoot_density")]
# shootDen <- shootDen[complete.cases(shootDen),] # n = 210; lost all 2014 data
# 
# sd <- rbind.fill(sd2014, shootDen)
# sd$year <- as.numeric(sd$year)
# sd$year0 <- as.numeric((sd$year-2014))
# 
# # write.csv(sd, "m2m/output/shootDensity2014_2018.csv", row.names = F)
