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
                   iter = 2000, warmup = 1000, chains=4,
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