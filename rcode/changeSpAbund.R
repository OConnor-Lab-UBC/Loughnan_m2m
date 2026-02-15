# Started Feb 5, 2026 by D. Loughnan

# aim of this code is to model changes in the abundance of the 
# 10 most abundant and ecoligcally importnat spp

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)
require(shinystan)
require(bayesplot)
require(brms)

setwd("~/Documents/github/Loughnan_m2m")

source("rcode/alphaDiv/prepData.R")

domSp <- c("Amphissa","Ampithoidae","Aoroides","Caprella","Exogone","Harpacticoida","Idoteidae",
            "Lacuna","Leptochelia","Mytilidae","Paguridae","Platynereis","Pugettia")

subAbund <- abund[abund$taxonL %in% domSp,]

abundData <- list(y = subAbund$abundance,
                    N = nrow(subAbund),
                    Nspp = length(unique(subAbund$taxonL)),
                    species = as.numeric(as.factor(subAbund$taxonL)),
                    year = subAbund$year-2014
                  )

mdlA <- stan("stan/spAbund.stan",
                   data = abundData,
                   iter = 4000, warmup = 3000, chains=4,
                   include = FALSE, pars = c("y_hat"))
            #, control = list(max_treedepth = 12, adapt_delta = 0.99))
save(mdlA, file = "output/mdlAbundChange.Rdata")

sumA <- summary(mdlA)$summary
y <- subAbund$abundance

post <- rstan::extract(mdlA)
y.ext <- post$y_rep

# pdf(file = "figures/mdl_densityplot_temp.pdf", width = 4, height = 4)
# par(mfrow = c(1,2))
ppc_dens_overlay(y, y.ext[1:500, ])

# dev.off()

#### plotting:
quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.25, 0.75,0.95))
  return(returnQuanilte)
}


bSp <- data.frame(sumA[grep("beta_group\\[", rownames(sumA)), c("mean","2.5%", "97.5%")])
aSp <- data.frame(sumA[grep("alpha_group\\[", rownames(sumA)), c("mean","2.5%", "97.5%")])
quant <- t(apply(data.frame(post), 2, quantile2575) )

mdlOut <- data.frame(species = sort(unique(subAbund$taxonL)))
mdlOut$bSp <- bSp[,"mean"]
mdlOut$aSp <- aSp[,"mean"]
ggplot(subAbund, aes(x = year, y = abundance)) + geom_point() + facet_wrap(~taxonL)

########################
# how many rare spp are there?

head(abund)

avgAbund <- aggregate(abund["abundance"],
                  abund[c("site", "taxonL")],
                  FUN = mean)


rare <- aggregate(avgAbund["abundance"],
                  avgAbund[c("site", "taxonL")],
                           FUN = sum)

rTab <- data.frame(unique(rare$site))
rTab$iter <- seq(1:length(unique(rare$site)))
rTab$noSp <- NA
rTab$perRare <- NA
names(rTab) <- c("Site","iter", "noSp", "perRare")

siteNam <- sort(unique(rare$site))

for(i in 1:length(unique(rare$site))){
rareTab <- unique(abund[,c("site")])
temp <- subset(rare, site == siteNam[i])
noSp <- length(unique(temp$taxonL))
tempRare <- subset(temp, abundance > 5)
noRare <- nrow(tempRare)
perRare <- round((noRare/noSp)*100,0)
rTab$noSp[rTab$iter == i] <- noSp
rTab$perRare[rTab$iter == i] <- perRare
}
rTab


