# aim of this code is to run a Bukley-terry model for the beta diversity model
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)
require(shinystan)
require(tidyverse)
require(vegan)
require(otuSummary)
require(bayesplot)

setwd("~/Documents/github/Calvert_O-Connor_eelgrass/m2m")

# get abiotic factors:
abiotic <- read.csv("output/alphaMacroOct212025.csv")
# abiotic$comm_id1 <- abiotic$site
# 
# abiotic$comm_id1[which(abiotic$comm_id1 == "choked_inner")] <- "chkd_"
# abiotic$comm_id1[which(abiotic$comm_id1 == "choked_sandspit")] <- "chkd_s"
# abiotic$comm_id1[which(abiotic$comm_id1 == "goose_south_east")] <- "gs_st_"
# abiotic$comm_id1[which(abiotic$comm_id1 == "goose_south_east")] <- "gs_s_"
# abiotic$comm_id1[which(abiotic$comm_id1 == "goose_south_west")] <- "gs_sth_w"
# abiotic$comm_id1[which(abiotic$comm_id1 == "mcmullin_north")] <- "mcmllns_n"
# abiotic$comm_id1[which(abiotic$comm_id1 == "mcmullin_south")] <- "mcmllns_s"
# abiotic$comm_id1[which(abiotic$comm_id1 == "choked_sandspit")] <- "cmcmllns_s"
# abiotic$comm_id1[which(abiotic$comm_id1 == "triquet_north")] <- "trqt_n"
# abiotic$comm_id1[which(abiotic$comm_id1 == "triquet_south")] <- "trqt_s"

abiotic$comm_id1 <- paste(abiotic$site, abiotic$quadrat_id, sep = "_")
abiotic <- abiotic[,c("region", "year", "site","quadrat_id","quadrat_macroalgae_g", "quadrat_lai","comboTemp")]
abiotic$comm_id1 <- paste(abiotic$site, abiotic$quadrat_id, sep = "_")
######## 2014 ########################################
ab2014 <- unique(subset(abiotic, year == "2014"))

bray2014 <- read.csv("output/betaDiv/betaDiv2014_macroeuk_braycurtis_finest.csv")
dim(bray2014)

rownames(bray2014) <- colnames(bray2014)

spec2014 <- matrixConvert(bray2014, colname = c("comm_id1", "comm_id2", "bray"))
spec2014$comm1Fact <- as.numeric(spec2014$comm_id1)
spec2014$comm2Fact <- as.numeric(spec2014$comm_id2)
spec2014 <- subset(spec2014, bray < 1)

lai2014 <- lai2014[complete.cases(lai2014$leaf_area_index),]
lai2014$comm_id1 <- paste(lai2014$site, lai2014$sample_no, sep = "_")
  
df<-data.frame(x=lai2014$leaf_area_index)
rownames(df)<- lai2014$comm_id
dist_lai <- dist(df)
lai_dist<-matrixConvert(dist_lai, colname = c("comm_id1", "comm_id2", "lai"))

data_model<-list(spec2014,lai_dist) %>% reduce(inner_join, by=c("comm_id1","comm_id2"))
# 
# spec2014A <- merge(spec2014, lai2014, by = "comm_id1")

temp <- str_split_fixed(data_model$comm_id1, "_", 3)
data_model$site <- paste(temp[,1],temp[,2], sep = "_")

dat2014<-list(N = nrow(data_model),
          N_sample = nrow(data_model),
          s = data_model$bray,
          regions = as.numeric(as.factor(data_model$site)),
          n_reg = length(unique(data_model$site)),
          x1= data_model$lai,
          # x2=x2,
          idx1= as.numeric(data_model$comm_id1),
          idx2 = as.numeric(data_model$comm_id2)
          )

# 
# mdlB2014 <- stan("stan/m2mBetaDivSite.stan",
#              data = dat2014,
#              iter = 4000, warmup = 3000, chains=4)

mdlB2014Simp <- stan("m2m/stan/m2mBetaDivSite.stan",
                 data = dat2014,
                 iter = 4000, warmup = 3000, chains=4)
save(mdlB, file = "output/mdlBeta2014.Rdata")

post <- rstan::extract(mdlB2014Simp)
pred<-post$pred
ppd_dens_overlay(pred[1:300,])+theme_bw()+legend_none()+labs(y="Density",x="Bray-Curtis indices")+
  labs(title="Prior samples")

s <- data_model$bray
# x1<-dat$x1
# x2<-dat$x2
mean_pred<-apply(pred,2,mean)
res<-post$res
res_mean<-apply(res,2,mean)

library(ggpubr)
library(bayesplot)
ggarrange(
  ppc_dens_overlay(s,pred[1:300,])+theme_bw()+legend_none()+
    labs(y="Density",x="Bray-Curtis indices")+
    theme(axis.title = element_text(size = 6)),
  
  ppc_intervals(res_mean,res,s)+
    labs(y="Residuals (Predicted - Observed)",
         x="Bray-Curtis indices") +
    theme_bw()+
    theme(axis.title = element_text(size = 6)) +
    hline_0()+legend_none(),
  
  ppc_intervals(mean_pred,pred,s)+labs(y="Predicted indices",
                                       x="Bray-Curtis indices")+
    theme_bw()+
    theme(axis.title = element_text(size = 6))+
    abline_01()+legend_none(),
  
  ppc_intervals(res_mean,res,x1)+labs(y="Residuals (Predicted - Observed)",
                                      x="Euc distance")+
    theme_bw()+
    theme(axis.title = element_text(size = 6))+
    hline_0()+legend_none(),
  
  ppc_intervals(res_mean,res,x2)+labs(y="Residuals (Predicted - Observed)",
                                      x="Precipitation difference")+
    theme_bw()+
    theme(axis.title = element_text(size = 6))+
    hline_0()+legend_none(),
  
  ncol=2,nrow=3)

ggarrange(
  mcmc_hist(fit$draws(c("beta1")))+theme_bw()+
    theme(axis.title = element_text(size = 6))+
    labs(x="Geographical distance slope",y="Frequency"),
  mcmc_hist(fit$draws(c("beta2")))+theme_bw()+
    theme(axis.title = element_text(size = 6))+
    labs(x="Elevation slope",y="Frequency"),
  ncol=1,nrow=2)


sum <- summary(mdlB2014Simp)$summary

# (c("alpha","kappa"),~quantile(., probs = c(0.05, 0.95)))

########### 2015 ##############################################
bray2015 <- read.csv("..//R_Code_and_Analysis/betadiversity/Bray-Curtis/2015_macroeuk_braycurtis_finest.csv")
dim(bray2015)

rownames(bray2015) <- colnames(bray2015)

spec2015 <- matrixConvert(bray2015, colname = c("comm_id1", "comm_id2", "bray"))
spec2015$comm1Fact <- as.numeric(spec2015$comm_id1)
spec2015$comm2Fact <- as.numeric(as.character(spec2015$comm_id2))

spec2015 <- subset(spec2015, bray < 1)
spec2015 <- subset(spec2015, bray > 0)

temp <- str_split_fixed(spec2015$comm_id1, "_", 3)
spec2015$site <- temp[,1]

dat2015<-list(N = nrow(spec2015),
              N_sample = nrow(spec2015),
              s = spec2015$bray,
              regions = as.numeric(as.factor(spec2015$site)),
              n_reg = length(unique(spec2015$site)),
              # x1=x1,
              # x2=x2,
              idx1= as.numeric(spec2015$comm_id1),
              idx2 = as.numeric(spec2015$comm_id2)
)

# 
# mdlB2015 <- stan("stan/m2mBetaDivSite.stan",
#              data = dat2015,
#              iter = 4000, warmup = 3000, chains=4)

mdlB2015Simp <- stan("stan/m2mBetaDiv.stan",
                     data = dat2015,
                     iter = 4000, warmup = 3000, chains=4)
save(mdlB2015Simp, file = "output/mdlBeta2015.Rdata")

post <- rstan::extract(mdlB2015Simp)
pred<-post$pred
ppd_dens_overlay(pred[1:300,])+theme_bw()+legend_none()+labs(y="Density",x="Bray-Curtis indices")


########### 2016 ##############################################

bray2016 <- read.csv("..//R_Code_and_Analysis/betadiversity/Bray-Curtis/2016_macroeuk_braycurtis_finest.csv")
dim(bray2016)

rownames(bray2016) <- colnames(bray2016)

spec2016 <- matrixConvert(bray2016, colname = c("comm_id1", "comm_id2", "bray"))
spec2016$comm1Fact <- as.numeric(spec2016$comm_id1)
spec2016$comm2Fact <- as.numeric(as.character(spec2016$comm_id2))

spec2016 <- subset(spec2016, bray < 1)
spec2016 <- subset(spec2016, bray > 0)
temp <- str_split_fixed(spec2016$comm_id1, "_", 3)
spec2016$site <- temp[,1]

dat2016<-list(N = nrow(spec2016),
              N_sample = nrow(spec2016),
              s = spec2016$bray,
              regions = as.numeric(as.factor(spec2016$site)),
              n_reg = length(unique(spec2016$site)),
              # x1=x1,
              # x2=x2,
              idx1= as.numeric(spec2016$comm_id1),
              idx2 = as.numeric(spec2016$comm_id2)
)

# 
# mdlB2016 <- stan("stan/m2mBetaDivSite.stan",
#              data = dat2016,
#              iter = 4000, warmup = 3000, chains=4)

mdlB2016Simp <- stan("stan/m2mBetaDiv.stan",
                     data = dat2016,
                     iter = 4000, warmup = 3000, chains=4)
save(mdlB2016Simp, file = "output/mdlBeta2016.Rdata")

post <- rstan::extract(mdlB2016Simp)
pred<-post$pred
ppd_dens_overlay(pred[1:300,])+theme_bw()+legend_none()+labs(y="Density",x="Bray-Curtis indices")
########### 2017 ##############################################

bray2017 <- read.csv("..//R_Code_and_Analysis/betadiversity/Bray-Curtis/2017_macroeuk_braycurtis_finest.csv")
dim(bray2017)

rownames(bray2017) <- colnames(bray2017)

spec2017 <- matrixConvert(bray2017, colname = c("comm_id1", "comm_id2", "bray"))
spec2017$comm1Fact <- as.numeric(spec2017$comm_id1)
spec2017$comm2Fact <- as.numeric(as.character(spec2017$comm_id2))

spec2017 <- subset(spec2017, bray < 1)
spec2017 <- subset(spec2017, bray > 0)

temp <- str_split_fixed(spec2017$comm_id1, "_", 3)
spec2017$site <- temp[,1]

dat2017<-list(N = nrow(spec2017),
              N_sample = nrow(spec2017),
              s = spec2017$bray,
              regions = as.numeric(as.factor(spec2017$site)),
              n_reg = length(unique(spec2017$site)),
              # x1=x1,
              # x2=x2,
              idx1= as.numeric(spec2017$comm_id1),
              idx2 = as.numeric(spec2017$comm_id2)
)

# 
# mdlB2017 <- stan("stan/m2mBetaDivSite.stan",
#              data = dat2017,
#              iter = 4000, warmup = 3000, chains=4)

mdlB2017Simp <- stan("stan/m2mBetaDiv.stan",
                     data = dat2017,
                     iter = 4000, warmup = 3000, chains=4)
save(mdlB2017Simp, file = "output/mdlBeta2017.Rdata")

post <- rstan::extract(mdlB2017Simp)
pred<-post$pred
ppd_dens_overlay(pred[1:300,])+theme_bw()+legend_none()+labs(y="Density",x="Bray-Curtis indices")

# data vis
col1.sp <-c( rgb(204 / 255, 105 / 255, 112 / 255, alpha = 0.8))
col2.sp <- c( rgb(205 / 255, 122 / 255, 0 / 255, alpha = 0.5))
# col2.sp <-c( rgb(0 / 255, 166 / 255, 167 / 255, alpha = 0.05))
col4.sp <- c( rgb(34 / 255, 166 / 255, 167 / 255, alpha = 0.5))
col5.sp <- c( rgb(141 / 255, 34 / 255, 171 / 255, alpha = 0.5))



par(mfrow = c(2,2))
hist(spec2014$bray, xlab = "2014 Bray-Curtis dissimilarity", ylab = "",main = NULL, col= col1.sp)
hist(spec2015$bray, xlab = "2015 Bray-Curtis dissimilarity", ylab = "",main = NULL, col= col2.sp)
hist(spec2016$bray, xlab = "2016 Bray-Curtis dissimilarity", ylab = "",main = NULL, col= col4.sp)
hist(spec2017$bray, xlab = "2017 Bray-Curtis dissimilarity", ylab = "",main = NULL, col= col5.sp)


library(ggpubr)
library(bayesplot)

post <- rstan::extract(mdlB2014Simp)
pred<-post$pred

s <- spec2014$bray
# x1<-dat$x1
# x2<-dat$x2
mean_pred<-apply(pred,2,mean)
res<-post$res
res_mean<-apply(res,2,mean)

plot2014 <- ppc_dens_overlay(s,pred[1:300,])+theme_bw()+legend_none()+
  labs(y="Density",x="Bray-Curtis indices")+
  theme(axis.text= element_text(size = 6))

post <- rstan::extract(mdlB2015Simp)
pred<-post$pred
s <- spec2015$bray
# x1<-dat$x1
# x2<-dat$x2
mean_pred<-apply(pred,2,mean)
res<-post$res
res_mean<-apply(res,2,mean)

plot2015 <- ppc_dens_overlay(s,pred[1:300,])+theme_bw()+legend_none()+
  labs(y="Density",x="Bray-Curtis indices")+
  theme(axis.text= element_text(size = 6))

post <- rstan::extract(mdlB2016Simp)
pred<-post$pred
s <- spec2016$bray
# x1<-dat$x1
# x2<-dat$x2
mean_pred<-apply(pred,2,mean)
res<-post$res
res_mean<-apply(res,2,mean)

plot2016 <- ppc_dens_overlay(s,pred[1:300,])+theme_bw()+legend_none()+
  labs(y="Density",x="Bray-Curtis indices")+
  theme(axis.text= element_text(size = 6))

post <- rstan::extract(mdlB2017Simp)
pred<-post$pred
s <- spec2017$bray
# x1<-dat$x1
# x2<-dat$x2
mean_pred<-apply(pred,2,mean)
res<-post$res
res_mean<-apply(res,2,mean)

plot2017 <- ppc_dens_overlay(s,pred[1:300,])+theme_bw()+legend_none()+
  labs(y="Density",x="Bray-Curtis indices")+
  theme(axis.text= element_text(size = 6))

pdf("figures/brayCurtisYpred.pdf", height = 5, width = 5)
plot_grid(plot2014, plot2015, plot2016, plot2017, labels = c('A', 'B', 'C', 'D'), label_size = 12)
dev.off()
ggarrange(
  ppc_dens_overlay(dat2014$s,pred[1:300,])+theme_bw()+legend_none()+
    labs(y="Density",x="Bray-Curtis indices")+
    theme(axis.title = element_text(size = 6)),
  
  ppc_dens_overlay(dat2014$bray,pred[1:300,])+theme_bw()+legend_none()+
    labs(y="Density",x="Bray-Curtis indices")+
    theme(axis.title = element_text(size = 6)),
  
  ppc_dens_overlay(dat2014$bray,pred[1:300,])+theme_bw()+legend_none()+
    labs(y="Density",x="Bray-Curtis indices")+
    theme(axis.title = element_text(size = 6)),
  
  ppc_dens_overlay(dat2014$bray,pred[1:300,])+theme_bw()+legend_none()+
    labs(y="Density",x="Bray-Curtis indices")+
    theme(axis.title = element_text(size = 6))+
    hline_0()+legend_none(),
  ncol=2, nrow=3)
