##################################################### Detrital food webs ############################
###################################################### Imbalances ##################################
#################################################### Rebecca Oester #################################
#################################################### September, 2023 ################################

rm(list=ls())
#### 1. Packages and Data ####
library(rstan)
library(brms)
library(bayesplot)
library(grid)
library(plyr)
library(dplyr)
library(outliers)
library(stringr)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(plyr)
library(vegan)
library(reshape2)
library(ggrepel)
library(ggExtra)
library(lme4)
library(sjPlot)
library(ggpubr)
library(ggExtra)

## color vectors
Leaves <- c( "firebrick3", "orange", "darkorange3") ## Alder, Ash, Mix
Shreds <- c( "#762a83", "#5ab4ac",  "#addd8e", "darkorange", "red") ## Crustacea, Plecoptera, Trichoptera
FNF <- c("springgreen4","burlywood3") ## Forest, Non-Forest
bw=c("white", "black")


## scaling function to standardize all numeric variables
scale2 <- function(x, na.rm = TRUE) ((x - mean(x, na.rm = na.rm)) / sd(x, na.rm))/2 ### na.rm = TRUE leaves NAs in the scaling and just considers values


## vaguely informative prior with normal distribution around 0
prior1 <- c(prior(normal(0, 1), class = Intercept),
            prior(normal(0, 0.1), class = "b"))
### load data ####
dat.all<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.all.txt")

### 2. d13C ####
###### 2.1 Crustacea ####
dat.all.crust<- subset(dat.all, Type=="Leaf" | Order=="Crustacea")
dat.all.crust <- subset(dat.all.crust, T=="sample")
dat.all.crust.C <- dat.all.crust %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Median_dC = median(Delta13_C_VPDB_corr), na.rm=TRUE)
## only get MZB and Leaves
dat.all.crust.C <- subset(dat.all.crust.C, (dat.all.crust.C$Type=="MZB" | dat.all.crust.C$Type=="Leaf") & dat.all.crust.C$T=="sample" & !is.na(dat.all.crust.C$Sample_ID))
## get only 1-2 values per bag
dat.all.crust.C <- dat.all.crust.C %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.crust.C$Leaf1 <- dat.all.crust.C$Leaf
dat.all.crust.C.dcast <- dcast(dat.all.crust.C, Sample_ID ~Type, value.var="Median_dC")
dat.all.crust.C.dcast$Diff_C_abs <- dat.all.crust.C.dcast$MZB-dat.all.crust.C.dcast$Leaf
dat.all.crust.C <- plyr::join(dat.all.crust.C.dcast, dat.all.crust.C, by="Sample_ID", type="inner")

dat.all.crust.C.dcast.1 <- dcast(dat.all.crust.C, Sample_ID~Type, value.var = "Median_dC")
dat.all.crust.C.dcast.1 <- dat.all.crust.C.dcast.1 %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.crust.C.dcast.1$Vegetation<- as.factor(dat.all.crust.C.dcast.1$Vegetation)
levels(dat.all.crust.C.dcast.1$Vegetation)
dat.all.crust.C.dcast.1 <- dat.all.crust.C.dcast.1 %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.crust.C.dcast.1 <- dat.all.crust.C.dcast.1 %>%
  mutate(Stream = substring(Sample_ID, 1, 3))


crustacead13c<- ggplot(data=subset(dat.all.crust.C.dcast.1))+
  geom_point(aes(x=Leaf, y=MZB,  alpha=Vegetation), size=3, color="#762a83")+
  #geom_smooth(aes(x=Leaf, y=MZB,color="#762a83"), method="lm",  se=FALSE, linetype = "dashed")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(delta^{13}, C," " , consumer, " [\u2030]")))+
  xlab(expression(paste(delta^{13}, C," " , resource, " [\u2030]")))+
  scale_y_continuous(limits=c(-32, -24))+
  scale_x_continuous(limits=c(-32, -28))
crustacead13c


aggregate(sqrt(Diff_C_abs) ~ Vegetation, data=subset(dat.all.crust.C, T=="sample"& Type=="MZB"), mean)

((1.914836-1.413293)/1.914836)*100## sqrt
((3.741591-2.073913)/3.741591)*100## not sqrt



##BRMS
dat.all.CmeanCrust <- subset(dat.all.crust.C, T=="sample"& Type=="MZB")

dat.all.CmeanCrustSTAN <-  dat.all.CmeanCrust %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Diff_C_abs) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Diff_C_abs), scale2) # scale numberic variables to 0


CmeanCrust <- brm(Diff_C_abs ~ Vegetation+(1|Stream),
                  data = dat.all.CmeanCrustSTAN,
                  prior = prior1)
summary(CmeanCrust)




###### 2.2 Plecoptera ####
dat.all.plec.C<- subset(dat.all, Type=="Leaf" | Order=="Plecoptera")
dat.all.plec.C <- subset(dat.all.plec.C, T=="sample")
dat.all.plec.C <- dat.all.plec.C %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Median_dC = median(Delta13_C_VPDB_corr), na.rm=TRUE)
## only get MZB and Leaves
dat.all.plec.C <- subset(dat.all.plec.C, (dat.all.plec.C$Type=="MZB" | dat.all.plec.C$Type=="Leaf") & dat.all.plec.C$T=="sample" & !is.na(dat.all.plec.C$Sample_ID))
## get only 1-2 values per bag
dat.all.plec.C <- dat.all.plec.C %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.plec.C$Leaf1 <- dat.all.plec.C$Leaf
dat.all.plec.C.dcast <- dcast(dat.all.plec.C, Sample_ID ~Type, value.var="Median_dC")
dat.all.plec.C.dcast$Diff_C_abs <- dat.all.plec.C.dcast$MZB-dat.all.plec.C.dcast$Leaf
dat.all.plec.C <- plyr::join(dat.all.plec.C.dcast, dat.all.plec.C, by="Sample_ID", type="inner")


dat.all.plec.C.dcast.1 <- dcast(dat.all.plec.C, Sample_ID~Type, value.var = "Median_dC")
dat.all.plec.C.dcast.1 <- dat.all.plec.C.dcast.1 %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.plec.C.dcast.1$Vegetation<- as.factor(dat.all.plec.C.dcast.1$Vegetation)
levels(dat.all.plec.C.dcast.1$Vegetation)
dat.all.plec.C.dcast.1 <- dat.all.plec.C.dcast.1 %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.plec.C.dcast.1 <- dat.all.plec.C.dcast.1 %>%
  mutate(Stream = substring(Sample_ID, 1, 3))


plecopterad13c<- ggplot(data=subset(dat.all.plec.C.dcast.1))+
  geom_point(aes(x=Leaf, y=MZB,  alpha=Vegetation), size=3, color="#5ab4ac")+
  #geom_smooth(aes(x=Leaf, y=MZB,color="#762a83"), method="lm",  se=FALSE, linetype = "dashed")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_text(size=18),
        axis.text=element_text(size=18),
        axis.title.y =element_blank(),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(delta^{13}, C," " , consumer, " [\u2030]")))+
  xlab(expression(paste(delta^{13}, C," " , resource, " [\u2030]")))+
  scale_y_continuous(limits=c(-32, -24))+
  scale_x_continuous(limits=c(-32, -28))
plecopterad13c


aggregate(sqrt(Diff_C_abs) ~ Vegetation, data=subset(dat.all.plec.C, T=="sample"& Type=="MZB"), mean)


((1.555988-1.497569)/1.555988)*100## not sqrt

## BRMS
dat.all.CmeanPlec <- subset(dat.all.plec.C, T=="sample"& Type=="MZB")

dat.all.CmeanPlecSTAN <-  dat.all.CmeanPlec %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Diff_C_abs) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Diff_C_abs), scale2) # scale numberic variables to 0


CmeanPlec <- brm(Diff_C_abs ~ Vegetation+(1|Region/Stream),
                 data = dat.all.CmeanPlecSTAN,
                 prior = prior1)
summary(CmeanPlec)

###### 2.3 Trichoptera ####
dat.all.tric.C<- subset(dat.all, Type=="Leaf" | Order=="Trichoptera")
dat.all.tric.C <- subset(dat.all.tric.C, T=="sample")
dat.all.tric.C <- dat.all.tric.C %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Median_dC = median(Delta13_C_VPDB_corr), na.rm=TRUE)
## only get MZB and Leaves
dat.all.tric.C <- subset(dat.all.tric.C, (dat.all.tric.C$Type=="MZB" | dat.all.tric.C$Type=="Leaf") & dat.all.tric.C$T=="sample" & !is.na(dat.all.tric.C$Sample_ID))
## get only 1-2 values per bag
dat.all.tric.C <- dat.all.tric.C %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.tric.C$Leaf1 <- dat.all.tric.C$Leaf
dat.all.tric.C.dcast <- dcast(dat.all.tric.C, Sample_ID ~Type, value.var="Median_dC")
dat.all.tric.C.dcast$Diff_C_abs <- dat.all.tric.C.dcast$MZB-dat.all.tric.C.dcast$Leaf
dat.all.tric.C <- plyr::join(dat.all.tric.C.dcast, dat.all.tric.C, by="Sample_ID", type="inner")


dat.all.tric.C.dcast.1 <- dcast(dat.all.tric.C, Sample_ID~Type, value.var = "Median_dC")
dat.all.tric.C.dcast.1 <- dat.all.tric.C.dcast.1 %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.tric.C.dcast.1$Vegetation<- as.factor(dat.all.tric.C.dcast.1$Vegetation)
levels(dat.all.tric.C.dcast.1$Vegetation)
dat.all.tric.C.dcast.1 <- dat.all.tric.C.dcast.1 %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.tric.C.dcast.1 <- dat.all.tric.C.dcast.1 %>%
  mutate(Stream = substring(Sample_ID, 1, 3))


trichopterad13c<- ggplot(data=subset(dat.all.tric.C.dcast.1))+
  geom_point(aes(x=Leaf, y=MZB,  alpha=Vegetation), size=3, color="#addd8e")+
  #geom_smooth(aes(x=Leaf, y=MZB,color="#762a83"), method="lm",  se=FALSE, linetype = "dashed")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title.y = element_blank(),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(delta^{13}, C," " , consumer, " [\u2030]")))+
  xlab(expression(paste(delta^{13}, C," " , resource, " [\u2030]")))+
  scale_y_continuous(limits=c(-32, -24))+
  scale_x_continuous(limits=c(-32, -28))
trichopterad13c


## BRMS
dat.all.CmeanTric <- subset(dat.all.tric.C, T=="sample"& Type=="MZB")

dat.all.CmeanTricSTAN <-  dat.all.CmeanTric %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Diff_C_abs) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Diff_C_abs), scale2) # scale numberic variables to 0


CmeanTric <- brm(Diff_C_abs ~ Vegetation+(1|Region/Stream),
                 data = dat.all.CmeanTricSTAN,
                 prior = prior1)
summary(CmeanTric)
### 3. C:N ####
###### 3.1 Crustacea ####
dat.all.crust<- subset(dat.all, Type=="Leaf" | Order=="Crustacea")
dat.all.crust <- subset(dat.all.crust, T=="sample")
dat.all.crust.CN <- dat.all.crust %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.crust.CN <- subset(dat.all.crust.CN, (dat.all.crust.CN$Type=="MZB" | dat.all.crust.CN$Type=="Leaf") & dat.all.crust.CN$T=="sample" & !is.na(dat.all.crust.CN$Sample_ID))
## get only 1-2 values per bag
dat.all.crust.CN <- dat.all.crust.CN %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.crust.CN$Leaf1 <- dat.all.crust.CN$Leaf
dat.all.crust.CN.dcast <- dcast(dat.all.crust.CN, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
dat.all.crust.CN.dcast$Ratio_molar_CN <- dat.all.crust.CN.dcast$Leaf/dat.all.crust.CN.dcast$MZB
#dat.all.crust.CN.dcast$LRR <- log(dat.all.crust.CN.dcast$Leaf/dat.all.crust.CN.dcast$MZB)
#dat.all.crust.CN.dcast$Diff_molar_CN <- dat.all.crust.CN.dcast$Leaf - dat.all.crust.CN.dcast$MZB
dat.all.crust.CN.dcast$Ratio_molar_CN_abs <- abs(dat.all.crust.CN.dcast$Ratio_molar_CN)
dat.all.crust.CN <- plyr::join(dat.all.crust.CN.dcast, dat.all.crust.CN, by="Sample_ID", type="inner")
dat.all.crust.CN.dcast.1 <- dcast(dat.all.crust.CN, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.crust.CN.dcast.1 <- dat.all.crust.CN.dcast.1 %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.crust.CN.dcast.1$Vegetation<- as.factor(dat.all.crust.CN.dcast.1$Vegetation)
levels(dat.all.crust.CN.dcast.1$Vegetation)
dat.all.crust.CN.dcast.1 <- dat.all.crust.CN.dcast.1 %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.crust.CN.dcast.1 <- dat.all.crust.CN.dcast.1 %>%
  mutate(Stream = substring(Sample_ID, 1, 3))


crustaceaCN<- ggplot(data=subset(dat.all.crust.CN.dcast.1))+
  geom_point(aes(x=(Leaf), y=(MZB),  alpha=Vegetation), size=3, color="#762a83")+
  #geom_smooth(aes(x=Leaf, y=MZB,color="#762a83"), method="lm",  se=FALSE, linetype = "dashed")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(C:N, " ", consumer)))+
  xlab(expression(paste(C:N, " ",resource)))+
  scale_y_continuous(limits=c(1, 2.5))+
  scale_x_continuous(limits=c(2, 3.5))
crustaceaCN

inletC <- ggplot(data=subset(dat.all.crust.CN.dcast.1))+
  geom_point(aes(x=(Leaf), y=(MZB),  alpha=Vegetation), size=3, color="#762a83")+
  geom_smooth(aes(x=Leaf, y=MZB,), method="lm",  se=FALSE, linetype = "dashed",color="#762a83")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=18),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(C:N, " ", consumer)))+
  xlab(expression(paste(C:N, " ",resource)))
inletC

ggsave("c_inlet.png", width = 4, height = 3)



aggregate(Ratio_molar_CN ~ Vegetation, data=subset(dat.all.crust.CN, T=="sample"& Type=="MZB"), mean)


((1.942700- 1.863295)/1.942700)*100## %log C:N mismatch difference


## BRMS
dat.all.CNCrust <- subset(dat.all.crust.CN, T=="sample"& Type=="MZB")

dat.all.CNCrustSTAN <-  dat.all.CNCrust %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Ratio_molar_CN) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Ratio_molar_CN), scale2) # scale numberic variables to 0


CNCrust <- brm(Ratio_molar_CN ~ Vegetation+(1|Stream),
               data = dat.all.CNCrustSTAN,
               prior = prior1)
summary(CNCrust)

###### 3.2 Plecoptera ####
dat.all.plec<- subset(dat.all, Type=="Leaf" | Order=="Plecoptera")
dat.all.plec <- subset(dat.all.plec, T=="sample")
dat.all.plec.CN <- dat.all.plec %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.plec.CN <- subset(dat.all.plec.CN, (dat.all.plec.CN$Type=="MZB" | dat.all.plec.CN$Type=="Leaf") & dat.all.plec.CN$T=="sample" & !is.na(dat.all.plec.CN$Sample_ID))
## get only 1-2 values per bag
dat.all.plec.CN <- dat.all.plec.CN %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.plec.CN$Leaf1 <- dat.all.plec.CN$Leaf
dat.all.plec.CN.dcast <- dcast(dat.all.plec.CN, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
dat.all.plec.CN.dcast$Ratio_molar_CN <- dat.all.plec.CN.dcast$Leaf/dat.all.plec.CN.dcast$MZB
#dat.all.plec.CN.dcast$LRR <- log(dat.all.plec.CN.dcast$Leaf/dat.all.plec.CN.dcast$MZB)
#dat.all.plec.CN.dcast$Diff_molar_CN <- dat.all.plec.CN.dcast$Leaf - dat.all.plec.CN.dcast$MZB
dat.all.plec.CN.dcast$Ratio_molar_CN_abs <- abs(dat.all.plec.CN.dcast$Ratio_molar_CN)
dat.all.plec.CN <- plyr::join(dat.all.plec.CN.dcast, dat.all.plec.CN, by="Sample_ID", type="inner")
dat.all.plec.CN.dcast.1 <- dcast(dat.all.plec.CN, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.plec.CN.dcast.1 <- dat.all.plec.CN.dcast.1 %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.plec.CN.dcast.1$Vegetation<- as.factor(dat.all.plec.CN.dcast.1$Vegetation)
levels(dat.all.plec.CN.dcast.1$Vegetation)
dat.all.plec.CN.dcast.1 <- dat.all.plec.CN.dcast.1 %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.plec.CN.dcast.1 <- dat.all.plec.CN.dcast.1 %>%
  mutate(Stream = substring(Sample_ID, 1, 3))



plecopteraCN<- ggplot(data=subset(dat.all.plec.CN.dcast.1))+
  geom_point(aes(x=(Leaf), y=(MZB),  alpha=Vegetation), size=3, color="#5ab4ac")+
  #geom_smooth(aes(x=Leaf, y=MZB,color="#762a83"), method="lm",  se=FALSE, linetype = "dashed")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_text(size=18),
        axis.text=element_text(size=18),
        axis.title.y = element_blank(),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(C:N, " ", consumer)))+
  xlab(expression(paste(C:N, " ",resource)))+
  scale_y_continuous(limits=c(1, 2.5))+
  scale_x_continuous(limits=c(2, 3.5))
plecopteraCN



inletP <- ggplot(data=subset(dat.all.plec.CN.dcast.1))+
  geom_point(aes(x=(Leaf), y=(MZB),  alpha=Vegetation), size=3, color="#5ab4ac")+
  geom_smooth(aes(x=Leaf, y=MZB,), method="lm",  se=FALSE, linetype = "dashed",color="#5ab4ac")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=18),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(C:N, " ", consumer)))+
  xlab(expression(paste(C:N, " ",resource)))
inletP

ggsave("p_inlet.png", width = 4, height = 3)


aggregate(Ratio_molar_CN ~ Vegetation, data=subset(dat.all.plec.CN, T=="sample"& Type=="MZB"), mean)


((2.027773- 1.990939)/2.027773)*100## %log C:N mismatch difference

## BRMS
dat.all.CNPlec <- subset(dat.all.plec.CN, T=="sample"& Type=="MZB")

dat.all.CNPlecSTAN <-  dat.all.CNPlec %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Ratio_molar_CN) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Ratio_molar_CN), scale2) # scale numberic variables to 0


CNPlec <- brm(Ratio_molar_CN ~ Vegetation+(1|Region/Stream),
              data = dat.all.CNPlecSTAN,
              prior = prior1)
summary(CNPlec)

###### 3.3 Trichoptera ####
dat.all.tric<- subset(dat.all, Type=="Leaf" | Order=="Trichoptera")
dat.all.tric <- subset(dat.all.tric, T=="sample")
dat.all.tric.CN <- dat.all.tric %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.tric.CN <- subset(dat.all.tric.CN, (dat.all.tric.CN$Type=="MZB" | dat.all.tric.CN$Type=="Leaf") & dat.all.tric.CN$T=="sample" & !is.na(dat.all.tric.CN$Sample_ID))
## get only 1-2 values per bag
dat.all.tric.CN <- dat.all.tric.CN %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.tric.CN$Leaf1 <- dat.all.tric.CN$Leaf
dat.all.tric.CN.dcast <- dcast(dat.all.tric.CN, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
dat.all.tric.CN.dcast$Ratio_molar_CN <- dat.all.tric.CN.dcast$Leaf/dat.all.tric.CN.dcast$MZB
#dat.all.tric.CN.dcast$LRR <- log(dat.all.tric.CN.dcast$Leaf/dat.all.tric.CN.dcast$MZB)
#dat.all.tric.CN.dcast$Diff_molar_CN <- dat.all.tric.CN.dcast$Leaf - dat.all.tric.CN.dcast$MZB
dat.all.tric.CN.dcast$Ratio_molar_CN_abs <- abs(dat.all.tric.CN.dcast$Ratio_molar_CN)
dat.all.tric.CN <- plyr::join(dat.all.tric.CN.dcast, dat.all.tric.CN, by="Sample_ID", type="inner")
dat.all.tric.CN.dcast.1 <- dcast(dat.all.tric.CN, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.tric.CN.dcast.1 <- dat.all.tric.CN.dcast.1 %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.tric.CN.dcast.1$Vegetation<- as.factor(dat.all.tric.CN.dcast.1$Vegetation)
levels(dat.all.tric.CN.dcast.1$Vegetation)
dat.all.tric.CN.dcast.1 <- dat.all.tric.CN.dcast.1 %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.tric.CN.dcast.1 <- dat.all.tric.CN.dcast.1 %>%
  mutate(Stream = substring(Sample_ID, 1, 3))


trichopteraCN<- ggplot(data=subset(dat.all.tric.CN.dcast.1))+
  geom_point(aes(x=(Leaf), y=(MZB),  alpha=Vegetation), size=3, color="#addd8e")+
  #geom_smooth(aes(x=Leaf, y=MZB,color="#762a83"), method="lm",  se=FALSE, linetype = "dashed")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title.y = element_blank(),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(C:N, " ", consumer)))+
  xlab(expression(paste(C:N, " ",resource)))+
  scale_y_continuous(limits=c(1, 2.5))+
  scale_x_continuous(limits=c(2, 3.5))
trichopteraCN

inletT <- ggplot(data=subset(dat.all.tric.CN.dcast.1))+
  geom_point(aes(x=(Leaf), y=(MZB),  alpha=Vegetation), size=3, color="#addd8e")+
  geom_smooth(aes(x=Leaf, y=MZB),color="#addd8e", method="lm",  se=FALSE, linetype = "dashed")+
  theme_bw()+
  scale_alpha_manual(values=c(0.4,1))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title.y = element_blank(),
        plot.title = element_text(size=18))+
  geom_abline(intercept = 0, slope=1)+
  ylab(expression(paste(C:N, " ", consumer)))+
  xlab(expression(paste(C:N, " ",resource)))
inletT
ggsave("t_inlet.png", width = 4, height = 3)


aggregate(Ratio_molar_CN ~ Vegetation, data=subset(dat.all.tric.CN, T=="sample"& Type=="MZB"), mean)


((1.731271- 1.675971)/1.731271)*100## %log C:N mismatch difference


## BRMS
dat.all.CNTric <- subset(dat.all.tric.CN, T=="sample"& Type=="MZB")

dat.all.CNTricSTAN <-  dat.all.CNTric %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Ratio_molar_CN) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Ratio_molar_CN), scale2) # scale numberic variables to 0


CNTric <- brm(Ratio_molar_CN ~ Vegetation+(1|Region/Stream),
              data = dat.all.CNTricSTAN,
              prior = prior1)
summary(CNTric)



#### 4. Plot ####

cc<-subset(dat.all.crust.C, T=="sample"& Type=="MZB")
cc <- select(cc, Sample_ID, Diff_C_abs, Vegetation, Order )
pc<-subset(dat.all.plec.C, T=="sample"& Type=="MZB")
pc <- select(pc, Sample_ID, Diff_C_abs, Vegetation, Order )
tc<-subset(dat.all.tric.C, T=="sample"& Type=="MZB")
tc <- select(tc, Sample_ID, Diff_C_abs, Vegetation , Order )
crust<- subset(dat.all.crust.CN, T=="sample"& Type=="MZB")
crust <- select(crust, Sample_ID, Ratio_molar_CN, Vegetation , Order )
plec<-subset(dat.all.plec.CN, T=="sample"& Type=="MZB")
plec <- select(plec, Sample_ID, Ratio_molar_CN, Vegetation , Order )
tric<- subset(dat.all.tric.CN, T=="sample"& Type=="MZB")
tric <- select(tric, Sample_ID, Ratio_molar_CN, Vegetation , Order )

a<-rbind( cc, pc, tc)
b<- rbind(crust, plec, tric)

ab<- merge(a,b)

bw<- c("white", "black")


c<-ggplot(ab)+theme_bw()+
  geom_boxplot(aes(x=Order, y=sqrt(Diff_C_abs), fill=Order,color=Order, alpha=Vegetation))+
  scale_fill_manual(values=Shreds)+
  scale_color_manual(values=Shreds)+
  scale_alpha_manual(values=c(0.4,0.7))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        plot.title = element_text(size=18))+
  ylab(expression(paste(delta^{13}, "C difference")))+
  annotate("text", y = 2.5, x= 1, label = "*", size=12)+
  annotate("text", y = 2.5, x= 2, label = "*", size=12)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
c

cn<-ggplot(ab)+theme_bw()+
  geom_boxplot(aes(x=Order, y=Ratio_molar_CN, fill=Order, color=Order, alpha=Vegetation))+
  scale_fill_manual(values=Shreds)+
  scale_color_manual(values=Shreds)+
  scale_alpha_manual(values=c(0.4,0.7))+
  guides(alpha= guide_legend(override.aes = list(fill= bw )))+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.3, 0.2), legend.title = element_text(size=16), legend.text = element_text(size=16), legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title.y = element_text(size=18),
        plot.title = element_text(size=18))+
  ylab(expression(paste("C:N difference")))+
  annotate("text", y = 2.2, x= 1, label = "*", size=12)+
  theme(axis.text.x = element_blank())
cn

ggarrange(cn, crustaceaCN, plecopteraCN, trichopteraCN,
          c, crustacead13c,plecopterad13c, trichopterad13c, ncol=4, nrow=2,
          align="h",
          widths= c(0.4, 0.3,0.3,0.3,0.4,0.3,0.3,0.3),
          labels=c("A", "B", "C", "D", "E", "F", "G", "H"))



ggsave("fig4_1.png", width = 16, height = 10)
