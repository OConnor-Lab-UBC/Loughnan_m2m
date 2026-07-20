# Started Feb 5, 2026 by D. Loughnan

# aim of this code is to model changes in the abundance of the 
# 13 most abundant and ecoligcally importnat spp

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)
require(shinystan)
require(bayesplot)
require(brms)
require(ggplot2)

setwd("~/Documents/github/Loughnan_m2m")

source("rcode/alphaDiv/prepData.R")

domSp <- c("Amphissa","Ampithoidae","Aoroides","Caprella","Exogone","Harpacticoida","Idoteidae",
            "Lacuna","Leptochelia","Mytilidae","Paguridae","Platynereis","Pugettia")

subAbund <- abund[abund$taxonL %in% domSp,]

subAbund <- subset(subAbund, site != "goose_north")

abundData <- list(y = subAbund$abundance,
                    N = nrow(subAbund),
                    Nspp = length(unique(subAbund$taxonL)),
                    Nsite = length(unique(subAbund$site)),
                    species = as.integer(as.factor(subAbund$taxonL)),
                    site = as.integer(as.factor(subAbund$site)),
                    year = subAbund$year-2014
                  )

mdlA <- stan("stan/spAbund.stan",
                   data = abundData,
                   iter = 4000, warmup = 3000, chains=4,
                   include = FALSE, pars = c("y_hat"))
            #, control = list(max_treedepth = 12, adapt_delta = 0.99))
save(mdlA, file = "output/mdlAbundChangeSite.Rdata")

sumA <- summary(mdlA)$summary
y <- subAbund$abundance

post <- rstan::extract(mdlA)
y.ext <- post$y_rep

# pdf(file = "figures/mdl_densityplot_temp.pdf", width = 4, height = 4)
# par(mfrow = c(1,2))
ppc_dens_overlay(y, y.ext[1:1000, ])

# dev.off()

#### plotting:
quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.25, 0.75,0.95))
  return(returnQuanilte)
}


bSp <- data.frame(sumA[grep("beta_group\\[", rownames(sumA)), c("mean")])
aSp <- data.frame(sumA[grep("alpha_group\\[", rownames(sumA)), c("mean","2.5%", "97.5%")])
aSite <- data.frame(sumA[grep("alpha_groupS\\[", rownames(sumA)), c("mean","2.5%", "97.5%")])

quant <- t(apply(data.frame(post), 2, quantile2575) )

mdlOut <- data.frame(species = sort(unique(subAbund$taxonL)))
mdlOut$bSp <- bSp[,"mean"]
mdlOut$aSp <- aSp[,"mean"]
ggplot(subAbund, aes(x = year, y = abundance)) + geom_point() + 
  facet_wrap(~taxonL, scales = "free")

# mu plot of species:
fitSp <- data.frame(post$beta_group)

bSpMean <- (data.frame(colMeans(fitSp)))
quantSpM <- t(apply(fitSp, 2, quantile2575) )
bSlopes <- data.frame(cbind(bSpMean, quantSpM)); names(bSlopes) <- c("mean","X5", "X25", "X75", "X95")
bSlopes$spId <- domSp

spAbund <- ggplot()  +
  geom_pointrange(bSlopes, mapping = aes(x = factor(spId, level = domSp), y = mean, ymin=X5, ymax= X95),
                  size =0.55, position=position_dodge(width=0.5), color = "maroon") +
  geom_linerange( bSlopes, mapping = aes(x = factor(spId, level = domSp), ymin = X25, ymax= X75), 
                  size =1.25, position=position_dodge(width=0.5), color = "maroon") +
  xlab("Dominant epifauna taxa") + ylab("Change in abundance per quadrat") +
 # annotate("text", 1 ,1, label = "a)", size = 5.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 12),axis.text.x = element_text( size= 10,angle = 55, hjust=1), legend.key=element_rect(fill="white"),legend.text=element_text(size=25)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")


pdf("figures/changeAbunda.pdf", width = 6, height = 4)
spAbund
dev.off()
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


