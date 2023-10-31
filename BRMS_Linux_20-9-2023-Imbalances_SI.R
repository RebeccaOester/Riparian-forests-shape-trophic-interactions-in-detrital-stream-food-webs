##################################################### Detrital food webs ############################
###################################################### Imbalances SI ##################################
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



### same for Subtraction####
### 4. Shred C:N ####
## min: 1.39 max: 2.08 log_molar_CN
###### 4.1 Crustacea ####
dat.all.crustS<- subset(dat.all, Type=="Leaf" | Order=="Crustacea")
dat.all.crustS <- subset(dat.all.crustS, T=="sample")
dat.all.crust.CNS <- dat.all.crustS %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.crust.CNS <- subset(dat.all.crust.CNS, (dat.all.crust.CNS$Type=="MZB" | dat.all.crust.CNS$Type=="Leaf") & dat.all.crust.CNS$T=="sample" & !is.na(dat.all.crust.CNS$Sample_ID))
## get only 1-2 values per bag
dat.all.crust.CNS <- dat.all.crust.CNS %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.crust.CNS$Leaf1 <- dat.all.crust.CNS$Leaf
dat.all.crust.CN.dcastS <- dcast(dat.all.crust.CNS, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
# dat.all.crust.CN.dcast$Ratio_molar_CN <- dat.all.crust.CN.dcast$Leaf/dat.all.crust.CN.dcast$MZB
# dat.all.crust.CN.dcast$LRR <- log(dat.all.crust.CN.dcast$Leaf/dat.all.crust.CN.dcast$MZB)
dat.all.crust.CN.dcastS$Diff_molar_CN <- dat.all.crust.CN.dcastS$Leaf - dat.all.crust.CN.dcastS$MZB
dat.all.crust.CN.dcastS$Diff_molar_CN_abs <- abs(dat.all.crust.CN.dcastS$Diff_molar_CN)
dat.all.crust.CNS <- plyr::join(dat.all.crust.CN.dcastS, dat.all.crust.CNS, by="Sample_ID", type="inner")
dat.all.crust.CN.dcast.1S <- dcast(dat.all.crust.CNS, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.crust.CN.dcast.1S <- dat.all.crust.CN.dcast.1S %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.crust.CN.dcast.1S$Vegetation<- as.factor(dat.all.crust.CN.dcast.1S$Vegetation)
levels(dat.all.crust.CN.dcast.1S$Vegetation)
dat.all.crust.CN.dcast.1S <- dat.all.crust.CN.dcast.1S %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.crust.CN.dcast.1S <- dat.all.crust.CN.dcast.1S %>%
  mutate(Stream = substring(Sample_ID, 1, 3))

aggregate(Diff_molar_CN ~ Vegetation, data=subset(dat.all.crust.CNS, T=="sample"& Type=="MZB"), mean)


((1.456143- 1.348445)/1.456143)*100## %log C:N mismatch difference


### BRMS ####
##
dat.all.CNCrust <- subset(dat.all.crust.CNS, T=="sample"& Type=="MZB")

dat.all.CNCrustSTAN <-  dat.all.CNCrust %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Diff_molar_CN) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Diff_molar_CN), scale2) # scale numberic variables to 0


CNCrust <- brm(Diff_molar_CN ~ Vegetation+(1|Stream),
               data = dat.all.CNCrustSTAN,
               prior = prior1)
summary(CNCrust)

###### 4.2 Plecoptera ####
dat.all.plecS<- subset(dat.all, Type=="Leaf" | Order=="Plecoptera")
dat.all.plecS <- subset(dat.all.plecS, T=="sample")
dat.all.plec.CNS <- dat.all.plecS %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.plec.CNS <- subset(dat.all.plec.CNS, (dat.all.plec.CNS$Type=="MZB" | dat.all.plec.CNS$Type=="Leaf") & dat.all.plec.CNS$T=="sample" & !is.na(dat.all.plec.CNS$Sample_ID))
## get only 1-2 values per bag
dat.all.plec.CNS <- dat.all.plec.CNS %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.plec.CNS$Leaf1 <- dat.all.plec.CNS$Leaf
dat.all.plec.CN.dcastS <- dcast(dat.all.plec.CNS, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
# dat.all.plec.CN.dcast$Ratio_molar_CN <- dat.all.plec.CN.dcast$Leaf/dat.all.plec.CN.dcast$MZB
# dat.all.plec.CN.dcast$LRR <- log(dat.all.plec.CN.dcast$Leaf/dat.all.plec.CN.dcast$MZB)
dat.all.plec.CN.dcastS$Diff_molar_CN <- dat.all.plec.CN.dcastS$Leaf - dat.all.plec.CN.dcastS$MZB
dat.all.plec.CN.dcastS$Diff_molar_CN_abs <- abs(dat.all.plec.CN.dcastS$Diff_molar_CN)
dat.all.plec.CNS <- plyr::join(dat.all.plec.CN.dcastS, dat.all.plec.CNS, by="Sample_ID", type="inner")
dat.all.plec.CN.dcast.1S <- dcast(dat.all.plec.CNS, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.plec.CN.dcast.1S <- dat.all.plec.CN.dcast.1S %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.plec.CN.dcast.1S$Vegetation<- as.factor(dat.all.plec.CN.dcast.1S$Vegetation)
levels(dat.all.plec.CN.dcast.1S$Vegetation)
dat.all.plec.CN.dcast.1S <- dat.all.plec.CN.dcast.1S %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.plec.CN.dcast.1S <- dat.all.plec.CN.dcast.1S %>%
  mutate(Stream = substring(Sample_ID, 1, 3))

aggregate(Diff_molar_CN ~ Vegetation, data=subset(dat.all.plec.CNS, T=="sample"& Type=="MZB"), mean)

((1.504219- 1.453524)/1.504219)*100## %log C:N mismatch difference


### BRMS ####
dat.all.CNPlec <- subset(dat.all.plec.CNS, T=="sample"& Type=="MZB")

dat.all.CNPlecSTAN <-  dat.all.CNPlec %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Diff_molar_CN) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Diff_molar_CN), scale2) # scale numberic variables to 0


CNPlec <- brm(Diff_molar_CN ~ Vegetation+(1|Stream),
              data = dat.all.CNPlecSTAN,
              prior = prior1)
summary(CNPlec)



###### 4.3 Trichoptera ####
dat.all.tricS<- subset(dat.all, Type=="Leaf" | Order=="Trichoptera")
dat.all.tricS <- subset(dat.all.tricS, T=="sample")
dat.all.tric.CNS <- dat.all.tricS %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.tric.CNS <- subset(dat.all.tric.CNS, (dat.all.tric.CNS$Type=="MZB" | dat.all.tric.CNS$Type=="Leaf") & dat.all.tric.CNS$T=="sample" & !is.na(dat.all.tric.CNS$Sample_ID))
## get only 1-2 values per bag
dat.all.tric.CNS <- dat.all.tric.CNS %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.tric.CNS$Leaf1 <- dat.all.tric.CNS$Leaf
dat.all.tric.CN.dcastS <- dcast(dat.all.tric.CNS, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
# dat.all.tric.CN.dcast$Ratio_molar_CN <- dat.all.tric.CN.dcast$Leaf/dat.all.tric.CN.dcast$MZB
# dat.all.tric.CN.dcast$LRR <- log(dat.all.tric.CN.dcast$Leaf/dat.all.tric.CN.dcast$MZB)
dat.all.tric.CN.dcastS$Diff_molar_CN <- dat.all.tric.CN.dcastS$Leaf - dat.all.tric.CN.dcastS$MZB
dat.all.tric.CN.dcastS$Diff_molar_CN_abs <- abs(dat.all.tric.CN.dcastS$Diff_molar_CN)
dat.all.tric.CNS <- plyr::join(dat.all.tric.CN.dcastS, dat.all.tric.CNS, by="Sample_ID", type="inner")
dat.all.tric.CN.dcast.1S <- dcast(dat.all.tric.CNS, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.tric.CN.dcast.1S <- dat.all.tric.CN.dcast.1S %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.tric.CN.dcast.1S$Vegetation<- as.factor(dat.all.tric.CN.dcast.1S$Vegetation)
levels(dat.all.tric.CN.dcast.1S$Vegetation)
dat.all.tric.CN.dcast.1S <- dat.all.tric.CN.dcast.1S %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.tric.CN.dcast.1S <- dat.all.tric.CN.dcast.1S %>%
  mutate(Stream = substring(Sample_ID, 1, 3))


### BRMS ####
dat.all.CNTric <- subset(dat.all.tric.CNS, T=="sample"& Type=="MZB")

dat.all.CNTricSTAN <-  dat.all.CNTric %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,Diff_molar_CN) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Diff_molar_CN), scale2) # scale numberic variables to 0


CNTric <- brm(Diff_molar_CN ~ Vegetation+(1|Stream),
              data = dat.all.CNTricSTAN,
              prior = prior1)
summary(CNTric)




### same for LRR ####
### 4. Shred C:N ####
## min: 1.39 max: 2.08 log_molar_CN

###### 4.1 Crustacea ####
dat.all.crustL<- subset(dat.all, Type=="Leaf" | Order=="Crustacea")
dat.all.crustL <- subset(dat.all.crustL, T=="sample")
dat.all.crust.CNL <- dat.all.crustL %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.crust.CNL <- subset(dat.all.crust.CNL, (dat.all.crust.CNL$Type=="MZB" | dat.all.crust.CNL$Type=="Leaf") & dat.all.crust.CNL$T=="sample" & !is.na(dat.all.crust.CNL$Sample_ID))
## get only 1-2 values per bag
dat.all.crust.CNL <- dat.all.crust.CNL %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.crust.CNL$Leaf1 <- dat.all.crust.CNL$Leaf
dat.all.crust.CN.dcastL <- dcast(dat.all.crust.CNL, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
#dat.all.crust.CN.dcast$Ratio_molar_CN <- dat.all.crust.CN.dcast$Leaf/dat.all.crust.CN.dcast$MZB
dat.all.crust.CN.dcastL$LRR <- log(dat.all.crust.CN.dcastL$Leaf/dat.all.crust.CN.dcastL$MZB)
#dat.all.crust.CN.dcast$Diff_molar_CN <- dat.all.crust.CN.dcast$Leaf - dat.all.crust.CN.dcast$MZB
dat.all.crust.CN.dcastL$LRR_abs <- abs(dat.all.crust.CN.dcastL$LRR)
dat.all.crust.CNL <- plyr::join(dat.all.crust.CN.dcastL, dat.all.crust.CNL, by="Sample_ID", type="inner")
dat.all.crust.CN.dcast.1L <- dcast(dat.all.crust.CNL, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.crust.CN.dcast.1L <- dat.all.crust.CN.dcast.1L %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.crust.CN.dcast.1L$Vegetation<- as.factor(dat.all.crust.CN.dcast.1L$Vegetation)
levels(dat.all.crust.CN.dcast.1L$Vegetation)
dat.all.crust.CN.dcast.1L <- dat.all.crust.CN.dcast.1L %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.crust.CN.dcast.1L <- dat.all.crust.CN.dcast.1L %>%
  mutate(Stream = substring(Sample_ID, 1, 3))



aggregate(LRR ~ Vegetation, data=subset(dat.all.crust.CNL, T=="sample"& Type=="MZB"), mean)

### BRMS ####

dat.all.CNCrust <- subset(dat.all.crust.CNL, T=="sample"& Type=="MZB")

dat.all.CNCrustSTANL <-  dat.all.CNCrust %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,LRR) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(LRR), scale2) # scale numberic variables to 0


CNCrustL <- brm(LRR ~ Vegetation+(1|Stream),
                data = dat.all.CNCrustSTANL,
                prior = prior1)
summary(CNCrustL)

###### 4.2 Plecoptera ####
dat.all.plecL<- subset(dat.all, Type=="Leaf" | Order=="Plecoptera")
dat.all.plecL <- subset(dat.all.plecL, T=="sample")
dat.all.plec.CNL <- dat.all.plecL %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.plec.CNL <- subset(dat.all.plec.CNL, (dat.all.plec.CNL$Type=="MZB" | dat.all.plec.CNL$Type=="Leaf") & dat.all.plec.CNL$T=="sample" & !is.na(dat.all.plec.CNL$Sample_ID))
## get only 1-2 values per bag
dat.all.plec.CNL <- dat.all.plec.CNL %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.plec.CNL$Leaf1 <- dat.all.plec.CNL$Leaf
dat.all.plec.CN.dcastL <- dcast(dat.all.plec.CNL, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
#dat.all.plec.CN.dcast$Ratio_molar_CN <- dat.all.plec.CN.dcast$Leaf/dat.all.plec.CN.dcast$MZB
dat.all.plec.CN.dcastL$LRR <- log(dat.all.plec.CN.dcastL$Leaf/dat.all.plec.CN.dcastL$MZB)
#dat.all.plec.CN.dcast$Diff_molar_CN <- dat.all.plec.CN.dcast$Leaf - dat.all.plec.CN.dcast$MZB
dat.all.plec.CN.dcastL$LRR_abs <- abs(dat.all.plec.CN.dcastL$LRR)
dat.all.plec.CNL <- plyr::join(dat.all.plec.CN.dcastL, dat.all.plec.CNL, by="Sample_ID", type="inner")
dat.all.plec.CN.dcast.1L <- dcast(dat.all.plec.CNL, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.plec.CN.dcast.1L <- dat.all.plec.CN.dcast.1L %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.plec.CN.dcast.1L$Vegetation<- as.factor(dat.all.plec.CN.dcast.1L$Vegetation)
levels(dat.all.plec.CN.dcast.1L$Vegetation)
dat.all.plec.CN.dcast.1L <- dat.all.plec.CN.dcast.1L %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.plec.CN.dcast.1L <- dat.all.plec.CN.dcast.1L %>%
  mutate(Stream = substring(Sample_ID, 1, 3))

aggregate(LRR ~ Vegetation, data=subset(dat.all.plec.CNL, T=="sample"& Type=="MZB"), mean)

### BRMS ####
dat.all.CNPlec <- subset(dat.all.plec.CNL, T=="sample"& Type=="MZB")

dat.all.CNPlecSTANL <-  dat.all.CNPlec %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,LRR) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(LRR), scale2) # scale numberic variables to 0


CNPlecL <- brm(LRR ~ Vegetation+(1|Stream),
               data = dat.all.CNPlecSTANL,
               prior = prior1)
summary(CNPlecL)


###### 4.3 Trichoptera ####
dat.all.tricL<- subset(dat.all, Type=="Leaf" | Order=="Trichoptera")
dat.all.tricL <- subset(dat.all.tricL, T=="sample")
dat.all.tric.CNL <- dat.all.tricL %>%
  dplyr::group_by(Sample_ID , Type) %>% ## so we have 1 value for the mixed bag leaves
  dplyr::mutate(Geom_Mean_molar_CN = median(log_molar_CN), na.rm=TRUE)
## only get MZB and Leaves
dat.all.tric.CNL <- subset(dat.all.tric.CNL, (dat.all.tric.CNL$Type=="MZB" | dat.all.tric.CNL$Type=="Leaf") & dat.all.tric.CNL$T=="sample" & !is.na(dat.all.tric.CNL$Sample_ID))
## get only 1-2 values per bag
dat.all.tric.CNL <- dat.all.tric.CNL %>% distinct(Type, Sample_ID, .keep_all = TRUE)
dat.all.tric.CNL$Leaf1 <- dat.all.tric.CNL$Leaf
dat.all.tric.CN.dcastL <- dcast(dat.all.tric.CNL, Sample_ID ~Type, value.var="Geom_Mean_molar_CN")
#dat.all.tric.CN.dcast$Ratio_molar_CN <- dat.all.tric.CN.dcast$Leaf/dat.all.tric.CN.dcast$MZB
dat.all.tric.CN.dcastL$LRR <- log(dat.all.tric.CN.dcastL$Leaf/dat.all.tric.CN.dcastL$MZB)
#dat.all.tric.CN.dcast$Diff_molar_CN <- dat.all.tric.CN.dcast$Leaf - dat.all.tric.CN.dcast$MZB
dat.all.tric.CN.dcastL$LRR_abs <- abs(dat.all.tric.CN.dcastL$LRR)
dat.all.tric.CNL <- plyr::join(dat.all.tric.CN.dcastL, dat.all.tric.CNL, by="Sample_ID", type="inner")
dat.all.tric.CN.dcast.1L <- dcast(dat.all.tric.CNL, Sample_ID~Type, value.var = "log_molar_CN")
dat.all.tric.CN.dcast.1L <- dat.all.tric.CN.dcast.1L %>%
  mutate(Vegetation = ifelse(str_detect(Sample_ID, "A"), "non-forested", "forested"))
dat.all.tric.CN.dcast.1L$Vegetation<- as.factor(dat.all.tric.CN.dcast.1L$Vegetation)
levels(dat.all.tric.CN.dcast.1L$Vegetation)
dat.all.tric.CN.dcast.1L <- dat.all.tric.CN.dcast.1L %>%
  mutate(Region = ifelse(str_detect(Sample_ID, "TG"), "North", "South"))
dat.all.tric.CN.dcast.1L <- dat.all.tric.CN.dcast.1L %>%
  mutate(Stream = substring(Sample_ID, 1, 3))

### BRMS ####
dat.all.CNTric <- subset(dat.all.tric.CNL, T=="sample"& Type=="MZB")

dat.all.CNTricSTANL <-  dat.all.CNTric %>%
  dplyr::select(Sample_ID,Stream, Region, Vegetation,LRR) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(LRR), scale2) # scale numberic variables to 0


CNTricL <- brm(LRR ~ Vegetation+(1|Stream),
               data = dat.all.CNTricSTANL,
               prior = prior1)
summary(CNTricL)




### 6. Plot SI C:N mismatch calculations####
all_sub <- rbind(dat.all.crust.CNS, dat.all.plec.CNS, dat.all.tric.CNS)

all_lrr <-rbind(dat.all.crust.CNL, dat.all.plec.CNL, dat.all.tric.CNL)

all_sub1 <- subset(all_sub, Order=="Crustacea" | Order=="Plecoptera" | Order=="Trichoptera")

all_lrr1 <- subset(all_lrr, Order=="Crustacea" | Order=="Plecoptera" | Order=="Trichoptera")


a<-ggplot(data=subset(all_sub1, Shredder=="Y"))+
  geom_boxplot(aes(y=Diff_molar_CN_abs, x=Order, fill=Order, color=Order, alpha=Vegetation))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_alpha_manual(values=c(0.4,0.7))+
  theme(axis.title.x=element_text(size=0),
        axis.text=element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15))+
  ylab("C:N difference (subtraction)")+
  scale_fill_manual(values=Shreds)+
  scale_color_manual(values=Shreds)+
  annotate("text", y = 1.7, x= 1, label = "*", size=8)
a


b<-ggplot(data=subset(all_lrr1, Shredder=="Y"))+
  geom_boxplot(aes(y=LRR, x=Order, fill=Order, color=Order, alpha=Vegetation))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_alpha_manual(values=c(0.4,0.7))+
  theme(axis.title.x=element_text(size=0),
        axis.text=element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15))+
  ylab("C:N difference (LRR)")+
  scale_fill_manual(values=Shreds)+
  scale_color_manual(values=Shreds)+
  annotate("text", y = 0.7, x= 1, label = "*", size=8)
b

library(patchwork)

a/b+ plot_annotation(tag_levels = 'A')
ggsave("SI_plot_differences.png", width = 4, height = 8)
