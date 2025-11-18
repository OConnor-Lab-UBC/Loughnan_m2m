# Aim of this code is to make the plots relevant for the m2m ms:
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(ggplot2)
require(dplyr)
library(tidyr)
library(rnaturalearth)
library(sf)
library(viridis)
library(PNWColors)

setwd("~/Documents/github/Loughnan_m2m")

# map climate data vs diversity sites
# coord <- read.csv("m2m/output/MacroSiteLatLong.csv")
coord <- read.csv("input/m2mSiteLatLong.csv")

# param <- c("Sites","Latitude", "Longitude", "Source")
# coord <- coord[,param]
# coord <- coord[complete.cases(coord),]
# coord$Source <- as.factor(coord$Source)
canada_map <- ne_states(country = 'canada', returnclass = "sf")

ggplot(canada_map) + 
  geom_sf() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("lat") + xlab("long") +
  coord_sf(xlim = c(-129, -127.5), ylim = c(51.5, 52.3)) +
  geom_point(data = coord, mapping = aes(x = long, y = lat))


### making figure of macro diversity for the relative abundance of the 20 most abundant taxa

graz <- read.csv("output/cleanedGrazer.csv")

length(unique(graz$taxonL))

graz$regionYr <- as.factor(paste(graz$site, graz$year, sep = "_"))
grazAg <- aggregate(graz["abundance"], graz[c("site","year","regionYr", "taxon5", "taxonL")], FUN = sum)
rYr <- (unique(graz$regionYr))

fiveP <- subset(grazAg, abundance > 10)
lessP <- subset(grazAg, abundance < 10)
lessP$taxon5 <- "Other"

lessP <- aggregate(lessP["abundance"], lessP[c("site","year","regionYr", "taxon5", "taxonL")], FUN = sum)


fiveAg <- rbind(fiveP, lessP)

fiveAg2 <- vector()

for(i in 1:length(rYr)){
  temp <- subset(fiveAg, regionYr == rYr[i])
  temp <- temp[order(temp$abundance),]
  temp$rank <- rep(1:nrow(temp))
  temp15 <- tail(temp, n = 15)
  sumSite <- sum(temp15$abundance)
  temp15$relAbund <- temp15$abundance/sumSite
  fiveAg2 <- rbind(fiveAg2, temp15)
}

rmSite <- c("goose_south_east","goose_north", "mcmullins_north", "mcmullins_south", "pruth_pocket")
fiveAg22 <- fiveAg2[!fiveAg2$site %in% rmSite,]

colour_taxa <- c("#180f3d",
                 "#721f81",
                 "#cd4071",
                 "#24878e",
                 "#fd9668",
                 "#34618d",
                 "#25ac82",
                 "#67cc5c",
                 "#fde725",
                 "#2fb47c",
                 "#98d83e",
                 "pink3",
                 "yellow" ,
                 "goldenrod1",
                 "#3d4d8a",
                 "#feca8d",
                 "#f1605d",
                 "#9e2f7f",
                 "#440f76",
                 "#440154",
                 "gray38")  

# taxon5Name<- c("Balanidae" ,"Bivalvia","Buccinoidea","Caprellidae",
#                "Cerithioidea", "Cnidaria", "Corophiidae", "Decapoda",            
#                "Eunicida", "Gammaridean", "Halacaridae", "Harpacticoida",       
#                "Isopoda","Littorinoidea","Ostracoda", "Phyllaplysia.taylori",
#                "Phyllodocida","Rissoidae", "Tanaidacea","Trochoidea","Other")
temp <- subset(fiveAg2, site == "choked_sandspit")

pdf("figures/relAbundGrazer.pdf", height = 7, width = 12)
grazerRA <- ggplot(fiveAg22, aes(year, relAbund,fill = taxonL)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack") + 
  scale_y_continuous(expand = c(0,0)) + #height of the column equal the value, stacked
  facet_wrap( ~ site, nrow = 2) + 
  ylab("Relative Abundance") +
  # scale_fill_manual(values=colour_taxa, labels = taxonL) + theme_bw() + ylab("Relative abundance") +
  theme(strip.background = element_rect(fill="snow2")) + 
  theme(panel.spacing = unit(0.2, "lines")) + 
  theme (axis.title.x=element_blank(),
                  #axis.text.x=element_blank(),
                  axis.text.x=element_text(angle = 90, hjust = 1, size = 10),
                  axis.ticks.x=element_blank(),
                  axis.text.y=element_text(size = 14), #change font size of numbers
                  axis.title.y=element_text(size = 18), #change font size of y title
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  legend.text=element_text(size=11),
         legend.key.height = unit(0.55, "cm"),
        legend.title = element_blank(),
                  plot.title = element_text(hjust = 0.5, size = 18, face="bold")) +
  guides(fill = guide_legend(ncol = 2)) #center plot title and set font size
grazerRA
dev.off()

temp <- subset(fiveAg2, regionYr == "triquet_south_2015")

# microeukaryotes

micro18 <- read.csv("data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)

m18 <- micro18 %>%
  pivot_longer(
    cols = starts_with("ASV"),  # Columns to pivot (can also be a range like Score_Test1:Score_Test3 or a vector of names)
    names_to = "asv",            # Name for the new column holding the original column names
    values_to = "abund"           # Name for the new column holding the values
  )

m18 <- subset(m18, abund > 0)

m18$regionYr <- as.factor(paste(m18$site, m18$year, sep = "_"))
m18Ag <- aggregate(m18["abund"], m18[c("site","year","regionYr", "asv")], FUN = sum)

rYr <- (unique(m18Ag$regionYr))

m18T <- vector()

for(i in 1:length(rYr)){
  temp <- subset(m18Ag, regionYr == rYr[i])
  temp <- temp[order(temp$abund),]
  temp$rank <- rep(1:nrow(temp))
  temp15 <- tail(temp, n = 15)
  sumSite <- sum(temp15$abund)
  temp15$relAbund <- temp15$abund/sumSite
  m18T <- rbind(m18T, temp15)
}

rmSite <- c("goose_south_east","goose_south_west", "mcmullins_north", "mcmullins_south")
m18T2 <- m18T[!m18T$site %in% rmSite,]

# taxon5Name<- c("Balanidae" ,"Bivalvia","Buccinoidea","Caprellidae",
#                "Cerithioidea", "Cnidaria", "Corophiidae", "Decapoda",            
#                "Eunicida", "Gammaridean", "Halacaridae", "Harpacticoida",       
#                "Isopoda","Littorinoidea","Ostracoda", "Phyllaplysia.taylori",
#                "Phyllodocida","Rissoidae", "Tanaidacea","Trochoidea","Other")
temp <- subset(fiveAg2, site == "choked_sandspit")

pdf("figures/relAbundMicro18.pdf", height = 7, width = 12)
micro18RA <- ggplot(m18T2, aes(year, relAbund,fill = asv)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack") + 
  scale_y_continuous(expand = c(0,0)) + #height of the column equal the value, stacked
  facet_wrap( ~ site, nrow = 2) + 
  ylab("Relative Abundance") +
  # scale_fill_manual(values=colour_taxa, labels = taxonL) + theme_bw() + ylab("Relative abundance") +
  theme(strip.background = element_rect(fill="snow2")) + 
  theme(panel.spacing = unit(0.2, "lines")) + 
  theme (axis.title.x=element_blank(),
         #axis.text.x=element_blank(),
         axis.text.x=element_text(angle = 90, hjust = 1, size = 10),
         axis.ticks.x=element_blank(),
         axis.text.y=element_text(size = 14), #change font size of numbers
         axis.title.y=element_text(size = 18), #change font size of y title
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         legend.text=element_text(size=11),
         legend.key.height = unit(0.55, "cm"),
         legend.title = element_blank(),
         plot.title = element_text(hjust = 0.5, size = 18, face="bold")) +
  guides(fill = guide_legend(ncol = 2)) #center plot title and set font size
micro18RA
dev.off()

# prokaryotes 



prok16 <- read.csv("data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)

p16 <- prok16 %>%
  pivot_longer(
    cols = starts_with("ASV"),  # Columns to pivot (can also be a range like Score_Test1:Score_Test3 or a vector of names)
    names_to = "asv",            # Name for the new column holding the original column names
    values_to = "abund"           # Name for the new column holding the values
  )

p16 <- subset(p16, abund > 0)

p16$regionYr <- as.factor(paste(p16$site, p16$year, sep = "_"))
p16Ag <- aggregate(p16["abund"], p16[c("site","year","regionYr", "asv")], FUN = sum)

rYr <- (unique(p16Ag$regionYr))

p16T <- vector()

for(i in 1:length(rYr)){
  temp <- subset(p16Ag, regionYr == rYr[i])
  temp <- temp[order(temp$abund),]
  temp$rank <- rep(1:nrow(temp))
  temp15 <- tail(temp, n = 15)
  sumSite <- sum(temp15$abund)
  temp15$relAbund <- temp15$abund/sumSite
  p16T <- rbind(p16T, temp15)
}

rmSite <- c("goose_south_east","goose_south_west", "mcmullins_north", "mcmullins_south")
p16T2 <- p16T[!p16T$site %in% rmSite,]

# taxon5Name<- c("Balanidae" ,"Bivalvia","Buccinoidea","Caprellidae",
#                "Cerithioidea", "Cnidaria", "Corophiidae", "Decapoda",            
#                "Eunicida", "Gammaridean", "Halacaridae", "Harpacticoida",       
#                "Isopoda","Littorinoidea","Ostracoda", "Phyllaplysia.taylori",
#                "Phyllodocida","Rissoidae", "Tanaidacea","Trochoidea","Other")
temp <- subset(fiveAg2, site == "choked_sandspit")

pdf("figures/relAbundprok16.pdf", height = 7, width = 12)
prok16RA <- ggplot(p16T2, aes(year, relAbund,fill = asv)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack") + 
  scale_y_continuous(expand = c(0,0)) + #height of the column equal the value, stacked
  facet_wrap( ~ site, nrow = 2) + 
  ylab("Relative Abundance") +
  # scale_fill_manual(values=colour_taxa, labels = taxonL) + theme_bw() + ylab("Relative abundance") +
  theme(strip.background = element_rect(fill="snow2")) + 
  theme(panel.spacing = unit(0.2, "lines")) + 
  theme (axis.title.x=element_blank(),
         #axis.text.x=element_blank(),
         axis.text.x=element_text(angle = 90, hjust = 1, size = 10),
         axis.ticks.x=element_blank(),
         axis.text.y=element_text(size = 14), #change font size of numbers
         axis.title.y=element_text(size = 18), #change font size of y title
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         legend.text=element_text(size=11),
         legend.key.height = unit(0.55, "cm"),
         legend.title = element_blank(),
         plot.title = element_text(hjust = 0.5, size = 18, face="bold")) +
  guides(fill = guide_legend(ncol = 2)) #center plot title and set font size
prok16RA
dev.off()



## NMDS plots:
microbes_16S_ASV <- read.csv("data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_16S_ASV)[1:16]

### Creating an object to store abundances only
abundances_16S_NMDS <- microbes_16S_ASV %>% 
  dplyr::select(-(1:15))

### Get MDS stats
set.seed(2)
NMDS.16S.LOG <- metaMDS(log(abundances_16S_NMDS+1), distance = "bray", k=2)  
NMDS.16S.LOG 

stressplot(NMDS.16S.LOG)
plot(NMDS.16S.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.16S.LOG$points[,1]
MDS2 = NMDS.16S.LOG$points[,2]
NMDS_16S = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = microbes_16S_ASV$year, region = microbes_16S_ASV$region)
NMDS_16S

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS_16S, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))
NMDS_16S
NMDS_16S$region

# re-order the factor levels before the plot
NMDS_16S$region <- factor(NMDS_16S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

NMDS_16S$year <- factor(NMDS_16S$year, levels=c("2015", "2016", "2017","2018"))

pdf("figures/nmdsProkaryotes.pdf", width = 6, height = 5)
nmds_prokaryotes <- ggplot(NMDS_16S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Prokaryotes") + 
  annotate("text", label = "stress = 0.20", x = 1.1, y = -1.5, size = 4, colour = "black") +
  scale_colour_manual(values=c("#1d457f",
                               "#625a94",
                               "#cc5c76",
                               "#f57946",
                               "#f9ad2a")) +
  scale_shape_manual(values=c(19,8,17,18)) + theme_bw()  +
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold"))
nmds_prokaryotes 
dev.off()

###
microbes_18S_ASV <- read.csv("data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_18S_ASV)[1:16]

### Creating an object to store abundances only
abundances_18S_NMDS <- microbes_18S_ASV %>% 
  dplyr::select(-(1:9))

### Get MDS stats
set.seed(2)
NMDS.18S.LOG <- metaMDS(log(abundances_18S_NMDS+1), distance = "bray", k=2)  
NMDS.18S.LOG 

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.18S.LOG$points[,1]
MDS2 = NMDS.18S.LOG$points[,2]
NMDS_18S = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = microbes_18S_ASV$year, region = microbes_18S_ASV$region)

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))

setnames(NMDS_18S, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))

# re-order the factor levels before the plot
NMDS_18S$region <- factor(NMDS_18S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

NMDS_18S$year <- factor(NMDS_18S$year, levels=c("2015", "2016", "2017","2018"))

pdf("figures/nmdsMicroeukaryotes.pdf", width = 6, height = 5)
nmds_microeukaryotes <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Microeukaryotes") + 
  annotate("text", label = "stress = 0.15", x = 1.3, y = -2.1, size = 4, colour = "black") +
  scale_colour_manual(values=c("#1d457f",
                               "#625a94",
                               "#cc5c76",
                               "#f57946",
                               "#f9ad2a")) +
  scale_shape_manual(values=c(19,8,17,18)) +  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold")) #center plot title and set font size
nmds_microeukaryotes 
dev.off()

### grazers 

inverts_finest <- read.csv("output/cleanedGrazer.csv")

### Creating an object to store abundances only
abundances_inverts_finest <- inverts_finest %>% 
  dplyr::select(-(1:5))

### Get MDS stats
set.seed(2)
NMDS.inverts.LOG <- metaMDS(log(abundances_inverts_finest+1), distance = "bray", k=2)  
NMDS.inverts.LOG 

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.inverts.LOG$points[,1]
MDS2 = NMDS.inverts.LOG$points[,2]
NMDS_inverts = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = inverts_finest$year, region = inverts_finest$region)
NMDS_inverts

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))

setnames(NMDS_inverts, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))
NMDS_inverts
NMDS_inverts$region

# re-order the factor levels before the plot
NMDS_inverts$region <- factor(NMDS_inverts$region, levels=c("choked", "pruth", "triquet","goose","mcmullin"))

NMDS_inverts$year <- factor(NMDS_inverts$year, levels=c("2014", "2015", "2016","2017"))

pdf("figures/nmdsInverts.pdf", width = 6, height = 5)
nmds_macroeukaryotes <- ggplot(NMDS_inverts, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macroeukaryotes") + 
  annotate("text", label = "stress = 0.19", x = - 1.8, y = -1.8, size = 4, colour = "black") +
  scale_colour_manual(values=c("#1d457f",
                               "#625a94",
                               "#cc5c76",
                               "#f57946",
                               "#f9ad2a")) +
  scale_shape_manual(values=c(0,19,8,17)) +  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold")) #center plot title and set font size
nmds_macroeukaryotes 
dev.off()


pal <- pnw_palette(name="Bay",n= 5,type="discrete")
ggplot(NMDS_inverts, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macroeukaryotes") + 
  annotate("text", label = "stress = 0.19", x = - 1.8, y = -1.8, size = 4, colour = "black") +
  scale_colour_manual(values= pal)
