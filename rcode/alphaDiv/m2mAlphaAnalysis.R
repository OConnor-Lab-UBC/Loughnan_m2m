# code started in May 2025 by D. Loughnan

# aim of this code is to run the joint bayesian analysis on the microbial and macrbiota data
# based on previous analyses and the availability of data the three primary parameters are: temperature, biomass, and macroAlgae

# June 29: met with Mary---suggested using LAI instead of biomass, including intercept for person who id'ed things
# add abundance, remove bed area, can we include epiphytes?
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)
require(shinystan)

setwd("~/Documents/github/Loughnan_m2m/rcode/alphaDiv")

source("prepData.R")

# try converting year from 1-4 first

## prok16
# the number of sites differ---does it matter? for now subsetting to be the same.

subSite <- c( "pruth_bay", "goose_north")
alpha16 <- alpha16[!alpha16$site %in% subSite,]
alpha18 <- alpha18[!alpha18$site %in% subSite,]

alphaMacro <- alphaFive
alphaMacro <- alphaMacro[!alphaMacro$site %in% subSite,]

alphaMacro$siteFactor <- as.numeric(as.factor(alphaMacro$site))
#alphaMacro$macroAm2 <- alphaMacro$quadrat_macroalgae_g * 16

#### Remove the first year and see if parameter esti change: only temp and LAI esti changed
# alphaMacro_3 <- subset(alphaMacro, year != "2014")


### Five parameter model 
# using LAI, temp, depth and macroAlgae instead of biomass or bed area
alphaMacro$actual_depth <- as.numeric(alphaMacro$actual_depth)

prok16.data <- list(yMicroi = alpha16$shannon,
                    N = nrow(alpha16),
                    n_site = length(unique(alpha16$site)),
                    micro_site = as.numeric(as.factor(alpha16$site)),
                    yeari = alpha16$year-2014,
                    Nm = nrow(alphaMacro),
                    n_siteM = length(unique(alphaMacro$site)),
                    macro_site = as.numeric(as.factor(alphaMacro$site)),
                    yMacroi = alphaMacro$shannon,
                    year2i = alphaMacro$year - 2014,
                    depthi = (alphaMacro$actual_depth-mean(alphaMacro$actual_depth))/(sd(alphaMacro$actual_depth)),
                    #  bedAi = (alphaMacro$bed_area_m2-mean(alphaMacro$bed_area_m2))/(sd(alphaMacro$bed_area_m2)),
                    LAIi = ((alphaMacro$quadrat_lai) - mean((alphaMacro$quadrat_lai)))/(sd((alphaMacro$quadrat_lai))),
                    tempi = (alphaMacro$comboTemp - mean(alphaMacro$comboTemp))/(sd(alphaMacro$comboTemp)),
                    macroAi = (alphaMacro$quadrat_macroalgae_g*16 - mean(alphaMacro$quadrat_macroalgae_g*16))/(sd(alphaMacro$quadrat_macroalgae_g*16)))



mdl4Shan16 <- stan("..//..//stan/jointTLDParam.stan",
                   data = prok16.data,
                   iter = 4000, warmup = 3000, chains=4,
                   include = FALSE, pars = c("y_hat"), control = list(max_treedepth = 12, adapt_delta = 0.99))
save(mdl4Shan16, file = "..//output/mdlShan16FourParaLAI.Rdata")


ssm <-  as.shinystan(macroAShan16)
launch_shinystan(ssm)

temp <- summary(macroAShan16)$summary

sumer16 <- round(summary(mdl4Shan16_3)$summary[c("mu_grand","muyear", "sigmayear", "sigma_microy","sigmaSt", 
                                              "mu_grandM", "muyear2", "sigmayear2","muDepth",
                                              "mutempSt","sigmatempSt","muLAISt","sigmaLAISt",
                                              # "mumacroASt","sigmamacroASt",
                                              "sigmaMacroSt","sigmaMacro_y",
                                              "betaMicroxtemp","betaMicroxLAI"
                                              # ,"betaMicroxmacroA"
                                              ), c("mean","2.5%","25%", "75%","97.5%")],2); sumer16


post <- rstan::extract(mdl4Shan16_3)

hist(rnorm(1000, 0, 120), col=rgb(0,0,1,1/4), xlim = c(-1000, 1000))
hist(post$mu_grand, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$muyear, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigmayear, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigma_microy, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-200,200))
hist(post$sigmaSt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$mu_grandM, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0, 20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$muyear2, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0, 20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$mutempSt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$muLAISt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$muDepth, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigmayear2, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigmatempSt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigmaLAISt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$sigmamacroASt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$sigmaMacroSt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,9), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$betaMicroxtemp, col=rgb(1,0,1,1/4), add = T)

# hist(rnorm(1000, 0, 20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
# hist(post$betaMicroxmacroA, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,9), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$betaMicroxLAI, col=rgb(1,0,1,1/4), add = T)


############ Micro 18 ######

micro18.data <- list(yMicroi = alpha18$shannon,
                     N = nrow(alpha18),
                     n_site = length(unique(alpha18$site)),
                     micro_site = as.numeric(as.factor(alpha18$site)),
                     yeari = alpha18$year-2014,
                     Nm = nrow(alphaMacro),
                     n_siteM = length(unique(alphaMacro$site)),
                     macro_site = as.numeric(as.factor(alphaMacro$site)),
                     yMacroi = alphaMacro$shannon,
                     year2i = alphaMacro$year-2014,
                     depthi = (alphaMacro$actual_depth-mean(alphaMacro$actual_depth))/(sd(alphaMacro$actual_depth)),
                     bedAi = (alphaMacro$bed_area_m2-mean(alphaMacro$bed_area_m2))/(sd(alphaMacro$bed_area_m2)),
                     LAIi = ((alphaMacro$quadrat_lai) - mean((alphaMacro$quadrat_lai)))/(sd((alphaMacro$quadrat_lai))),
                     tempi = (alphaMacro$comboTemp - mean(alphaMacro$comboTemp))/(sd(alphaMacro$comboTemp)),
                     macroAi = (alphaMacro$quadrat_macroalgae_g - mean(alphaMacro$quadrat_macroalgae_g))/(sd(alphaMacro$quadrat_macroalgae_g)))

mdl4Shan18 <- stan("..//stan/jointTLDParam.stan",
                   data = micro18.data,
                   iter = 6000, warmup = 3000, chains=4,
                   include = FALSE, pars = c("y_hat"), control = list(max_treedepth = 12, adapt_delta = 0.99))

save(mdl4Shan18, file = "..//output/mdlShan18ParamTLD_temp.Rdata")

ssm <-  as.shinystan(mdl4Shan18)
launch_shinystan(ssm)

temp <- summary(mdl4Shan18)$summary
sumer18 <- round(summary(mdl4Shan18)$summary[c("mu_grand","muyear", "sigmayear", "sigma_microy","sigmaSt", 
                                               "mu_grandM", "muyear2", "sigmayear2","muDepth",
                                               "mutempSt","sigmatempSt","mubiomassSt","sigmabiomassSt",
                                               # "mumacroASt","sigmamacroASt",
                                               "sigmaMacroSt","sigmaMacro_y",
                                               "betaMicroxtemp","betaMicroxbiomass"
),c("mean","2.5%","25%", "75%","97.5%","Rhat","n_eff")],2); sumer18


# load("..//output/mdlShan18ParamTLD.Rdata")
post <- rstan::extract(mdl4Shan18)

hist(rnorm(1000, 0, 20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$mu_grand, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$muyear, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigmayear, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigma_microy, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-200,200))
hist(post$sigmaSt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$mu_grandM, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0, 20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$muyear2, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0, 20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$mutempSt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$muLAISt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$muDepth, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigmayear2, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigmatempSt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-100, 100))
hist(post$sigmaLAISt, col=rgb(1,0,1,1/4), add = T)

# hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
# hist(post$sigmamacroASt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$sigmaMacroSt, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$betaMicroxtemp, col=rgb(1,0,1,1/4), add = T)
# 
# hist(rnorm(1000, 0, 20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
# hist(post$betaMicroxmacroA, col=rgb(1,0,1,1/4), add = T)

hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4), xlim = c(-80, 80))
hist(post$betaMicroxlai, col=rgb(1,0,1,1/4), add = T)

### Other diversity indixes #############################################

# eveness/ Pielou
prok16R.data <- list(yMicroi = alpha16$richness,
                    N = nrow(alpha16),
                    n_site = length(unique(alpha16$site)),
                    micro_site = as.numeric(as.factor(alpha16$site)),
                    yeari = alpha16$year-2014,
                    Nm = nrow(alphaMacro),
                    n_siteM = length(unique(alphaMacro$site)),
                    macro_site = as.numeric(as.factor(alphaMacro$site)),
                    yMacroi = alphaMacro$richness,
                    year2i = alphaMacro$year - 2014,
                    depthi = (alphaMacro$actual_depth-mean(alphaMacro$actual_depth))/(sd(alphaMacro$actual_depth)),
                    #  bedAi = (alphaMacro$bed_area_m2-mean(alphaMacro$bed_area_m2))/(sd(alphaMacro$bed_area_m2)),
                    LAIi = ((alphaMacro$quadrat_lai) - mean((alphaMacro$quadrat_lai)))/(sd((alphaMacro$quadrat_lai))),
                    tempi = (alphaMacro$comboTemp - mean(alphaMacro$comboTemp))/(sd(alphaMacro$comboTemp)),
                    macroAi = (alphaMacro$quadrat_macroalgae_g*16 - mean(alphaMacro$quadrat_macroalgae_g*16))/(sd(alphaMacro$quadrat_macroalgae_g*16)))


mdl4Rich16 <- stan("..//stan/jointTLDParam.stan",
                   data = prok16R.data,
                   iter = 5000, warmup = 3000, chains=4,
                   include = FALSE, pars = c("y_hat"), control = list(max_treedepth = 12, adapt_delta = 0.99))

save(mdl4Rich16, file = "..//output/mdlRich16ParamTLD.Rdata")

ssm <-  as.shinystan(mdl4Shan16)
launch_shinystan(ssm)

sumer16Rich <- round(summary(mdl4Rich16)$summary[c("mu_grand","muyear", "sigmayear", "sigma_microy","sigmaSt", 
                                               "mu_grandM", "muyear2", "sigmayear2","muDepth",
                                               "mutempSt","sigmatempSt","muLAISt","sigmaLAISt",
                                               "sigmaMacroSt","sigmaMacro_y",
                                               "betaMicroxtemp","betaMicroxLAI"),c("mean","2.5%","25%", "75%","97.5%","Rhat","n_eff")],2); sumer16Rich
                                               
micro18R.data <- list(yMicroi = alpha18$richness,
                     N = nrow(alpha18),
                     n_site = length(unique(alpha18$site)),
                     micro_site = as.numeric(as.factor(alpha18$site)),
                     yeari = alpha18$year-2014,
                     Nm = nrow(alphaMacro),
                     n_siteM = length(unique(alphaMacro$site)),
                     macro_site = as.numeric(as.factor(alphaMacro$site)),
                     yMacroi = alphaMacro$richness,
                     year2i = alphaMacro$year-2014,
                     depthi = (alphaMacro$actual_depth-mean(alphaMacro$actual_depth))/(sd(alphaMacro$actual_depth)),
                     # bedAi = (alphaMacro$bed_area_m2-mean(alphaMacro$bed_area_m2))/(sd(alphaMacro$bed_area_m2)),
                     LAIi = ((alphaMacro$quadrat_lai) - mean((alphaMacro$quadrat_lai)))/(sd((alphaMacro$quadrat_lai))),
                     tempi = (alphaMacro$comboTemp - mean(alphaMacro$comboTemp))/(sd(alphaMacro$comboTemp)))
#macroAi = (alphaMacro$quadrat_macroalgae_g - mean(alphaMacro$quadrat_macroalgae_g))/(sd(alphaMacro$quadrat_macroalgae_g)))

mdl4Rich18 <- stan("..//stan/jointTLDParam.stan",
                   data = micro18R.data,
                   iter = 4000, warmup = 3000, chains=4,
                   include = FALSE, pars = c("y_hat"), control = list(max_treedepth = 12, adapt_delta = 0.99))

save(mdl4Rich18, file = "..//output/mdlRich18ParamTLD.Rdata")

sumer18Rich <- round(summary(mdl4Rich18)$summary[c("mu_grand","muyear", "sigmayear", "sigma_microy","sigmaSt", 
                                                   "mu_grandM", "muyear2", "sigmayear2","muDepth",
                                                   "mutempSt","sigmatempSt","muLAISt","sigmaLAISt",
                                                   "sigmaMacroSt","sigmaMacro_y",
                                                   "betaMicroxtemp","betaMicroxLAI"),c("mean","2.5%","25%", "75%","97.5%","Rhat","n_eff")],2); sumer18Rich


ssm <-  as.shinystan(mdl4Shan18)
launch_shinystan(ssm)


# remove the 2015 macro data and see if it changes any of our parameter estimates:
alphaMacro_3$actual_depth <- as.numeric(alphaMacro_3$actual_depth)

prok16.data.3yr <- list(yMicroi = alpha16$shannon,
                    N = nrow(alpha16),
                    n_site = length(unique(alpha16$site)),
                    micro_site = as.numeric(as.factor(alpha16$site)),
                    yeari = alpha16$year-2014,
                    Nm = nrow(alphaMacro_3),
                    n_siteM = length(unique(alphaMacro_3$site)),
                    macro_site = as.numeric(as.factor(alphaMacro_3$site)),
                    yMacroi = alphaMacro_3$shannon,
                    year2i = alphaMacro_3$year - 2014,
                    depthi = (alphaMacro_3$actual_depth-mean(alphaMacro_3$actual_depth))/(sd(alphaMacro_3$actual_depth)),
                    #  bedAi = (alphaMacro$bed_area_m2-mean(alphaMacro$bed_area_m2))/(sd(alphaMacro$bed_area_m2)),
                    LAIi = ((alphaMacro_3$quadrat_lai) - mean((alphaMacro_3$quadrat_lai)))/(sd((alphaMacro_3$quadrat_lai))),
                    tempi = (alphaMacro_3$comboTemp - mean(alphaMacro_3$comboTemp))/(sd(alphaMacro_3$comboTemp)),
                    macroAi = (alphaMacro_3$quadrat_macroalgae_g*16 - mean(alphaMacro_3$quadrat_macroalgae_g*16))/(sd(alphaMacro_3$quadrat_macroalgae_g*16)))


mdl4Shan16_3yr <- stan("m2m/stan/jointTLDParam.stan",
                   data = prok16.data.3yr,
                   iter = 4000, warmup = 3000, chains=4,
                   include = FALSE, pars = c("y_hat"), control = list(max_treedepth = 12, adapt_delta = 0.99))

save(mdl4Shan16_3yr, file = "m2m/output/mdlShan16ParamTLD_3yrMacro.Rdata")

###############################################################################################
micro18.data_3yr <- list(yMicroi = alpha18$shannon,
                     N = nrow(alpha18),
                     n_site = length(unique(alpha18$site)),
                     micro_site = as.numeric(as.factor(alpha18$site)),
                     yeari = alpha18$year-2014,
                     Nm = nrow(alphaMacro_3),
                     n_siteM = length(unique(alphaMacro_3$site)),
                     macro_site = as.numeric(as.factor(alphaMacro_3$site)),
                     yMacroi = alphaMacro_3$shannon,
                     year2i = alphaMacro_3$year - 2014,
                     depthi = (alphaMacro_3$actual_depth-mean(alphaMacro_3$actual_depth))/(sd(alphaMacro_3$actual_depth)),
                     #  bedAi = (alphaMacro$bed_area_m2-mean(alphaMacro$bed_area_m2))/(sd(alphaMacro$bed_area_m2)),
                     LAIi = ((alphaMacro_3$quadrat_lai) - mean((alphaMacro_3$quadrat_lai)))/(sd((alphaMacro_3$quadrat_lai))),
                     tempi = (alphaMacro_3$comboTemp - mean(alphaMacro_3$comboTemp))/(sd(alphaMacro_3$comboTemp)),
                     macroAi = (alphaMacro_3$quadrat_macroalgae_g*16 - mean(alphaMacro_3$quadrat_macroalgae_g*16))/(sd(alphaMacro_3$quadrat_macroalgae_g*16)))

mdl4Shan18_3yr <- stan("..//stan/jointTLDParam.stan",
                   data = micro18.data_3yr,
                   iter = 8000, warmup = 7000, chains=4,
                   include = FALSE, pars = c("y_hat"), control = list(max_treedepth = 12, adapt_delta = 0.99))

save(mdl4Shan18_3yr, file = "..//output/mdlShan18ParamTLD_3yrMacro.Rdata")
