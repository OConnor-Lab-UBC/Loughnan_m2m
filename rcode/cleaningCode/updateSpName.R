# January 21, 2026 
# aim of this code is to better clean the grazer taxa names

setwd("~/Documents/github/Loughnan_m2m")

fa <- read.csv("data/Final_Species_Master_List_June 30_2025.csv", na.strings = c(""))
fa <- fa[,c("Original.Species.Name","Group.Classification","Final.ID")]

fa$Final.ID <- str_remove(fa$Final.ID, ".sp.+$") # nspp = 139
fa$Final.ID <- gsub( " ", ".", fa$Final.ID)

colnames(fa)[colnames(fa) == "Final.ID"] <- "taxon"


## The M2M data:
m <- read.csv( "data/MASTER_grazers.csv")

temp <-  aggregate(m["taxon"], m[c("year","site")], FUN = length)
temp <-  aggregate(temp["taxon"], temp[c("site")], FUN = mean)
# n= 193 

m <- m[,2:7]
m$taxon[which(m$taxon == "Porcellidium\xa0 ")] <- "Porcellidium"

m$taxon[which(m$taxon == "")] <- NA
m$taxon <- gsub( " ", ".", m$taxon )

m$taxon[which(m$taxon == "Gammaidean.Amphipod...1mm.")] <- "Gammaidean"
m$taxon[which(m$taxon == "Harpacticoid.copepod..multiple.species.")] <- "Harpacticoida"
m$taxon[which(m$taxon == "Lottia.alveus.paralella")] <- "Lottia alveus"
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

# cleaning spp names and lumping:
d$taxon <- str_remove(d$taxon, ".sp.+$") # nspp = 139
d <- unique(d[,c("year","site","taxon")])

# merge and see the differences:

taxaN <- merge(d, fa, by = "taxon", all.x = TRUE)
unmatch <- subset(taxaN, is.na(Original.Species.Name))
matched <- subset(taxaN, !is.na(Original.Species.Name))

missing <- unique(unmatch$taxon)
