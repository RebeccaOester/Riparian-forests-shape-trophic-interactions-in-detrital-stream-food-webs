##################################################### Detrital food webs ############################
###################################################### SI additional analyses ######################
#################################################### Rebecca Oester #################################
#################################################### September, 2023 ################################

rm(list=ls())

###load packages ####
library(dplyr)
library(stringr)
library(tidyr)
library(plyr)
library(reshape2)
library(ggrepel)
library(lme4)
library(sjPlot)
library(rstan)
library(brms)
library(bayesplot)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(outliers)
library(stringr)
library(ggpubr)
library(patchwork)

Leaves <- c( "firebrick3", "orange", "darkorange3") ## Alder, Ash, Mix
Shreds <- c( "#762a83", "#5ab4ac",  "#addd8e", "darkorange", "red") ## Crustacea, Plecoptera, Trichoptera

#### scaling function
scale2 <- function(x, na.rm = TRUE) ((x - mean(x, na.rm = na.rm)) / sd(x, na.rm))/2 ### na.rm = TRUE leaves NAs in the scaling and just considers values

#### prior
prior1 <- c(prior(normal(0, 1), class = Intercept),
            prior(normal(0, 0.1), class = "b"))

### load data ####
dat.all<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.all.txt")
dat.tot<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.tot.txt")

### 1.Leaves ####
dat.all$Vegetation<- factor(dat.all$Vegetation, levels=(c("T0", "forested", "non-forested")))

### Alder BRMS
Alder <- subset(dat.all, LeafLeafLeaf=="Alder")
fungiSTAN_Alder <-  Alder %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region,T, Vegetation, EA_C_Percent_corr,EA_N_Percent_corr, Delta15_N_Air_corr) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr), scale2) # scale numberic variables to 0

## C
T0_Alder<- brm(EA_C_Percent_corr ~ T, data = fungiSTAN_Alder, prior=prior1)
summary(T0_Alder)
pp_check(T0_Alder)

##N
T0_AlderN<- brm(EA_N_Percent_corr ~ T, data = fungiSTAN_Alder, prior=prior1)
summary(T0_AlderN)
pp_check(T0_AlderN)


##d15N
T0_Alder15N<- brm(Delta15_N_Air_corr ~ T, data = fungiSTAN_Alder, prior=prior1)
summary(T0_Alder15N)
pp_check(T0_Alder15N)


## Ash
Ash <- subset(dat.all, LeafLeafLeaf=="Ash")
fungiSTAN_Ash <-  Ash %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region,T, Vegetation, EA_C_Percent_corr,EA_N_Percent_corr, Delta15_N_Air_corr) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(EA_C_Percent_corr,EA_N_Percent_corr, Delta15_N_Air_corr), scale2) # scale numberic variables to 0

## C
T0_Ash<- brm(EA_C_Percent_corr ~ T, data = fungiSTAN_Ash, prior=prior1)
summary(T0_Ash)
pp_check(T0_Ash)

##N
T0_AshN<- brm(EA_N_Percent_corr ~ T, data = fungiSTAN_Ash, prior=prior1)
summary(T0_AshN)
pp_check(T0_AshN)

##d15N
T0_Ash15N<- brm(Delta15_N_Air_corr ~ T, data = fungiSTAN_Ash, prior=prior1)
summary(T0_Ash15N)
pp_check(T0_Ash15N)

### BRMS
## Alder
AlderS <- subset(dat.all, LeafLeafLeaf=="Alder" & T=="sample")
fungiSTAN_Alder <-  AlderS %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region,T, Vegetation, EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr), scale2) # scale numberic variables to 0

V_Alder<- brm(EA_C_Percent_corr ~ Vegetation+(1|Region/Stream),
              data = fungiSTAN_Alder, prior = prior1)
summary(V_Alder)
pp_check(V_Alder)


V_AlderN<- brm(EA_N_Percent_corr ~ Vegetation+(1|Region/Stream),
               data = fungiSTAN_Alder, prior = prior1)
summary(V_AlderN)
pp_check(V_AlderN)


V_Alder15N<- brm(Delta15_N_Air_corr ~ Vegetation+(1|Region/Stream),
                 data = fungiSTAN_Alder, prior = prior1)
summary(V_Alder15N)
pp_check(V_Alder15N)

## Ash
AshS <- subset(dat.all, LeafLeafLeaf=="Ash" & T=="sample")
fungiSTAN_Ash <-  AshS %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region,T, Vegetation, EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr), scale2) # scale numberic variables to 0

V_Ash<- brm(EA_C_Percent_corr ~ Vegetation+(1|Region/Stream),
            data = fungiSTAN_Ash, prior = prior1)
summary(V_Ash)
pp_check(V_Ash)


V_AshN<- brm(EA_N_Percent_corr ~ Vegetation+(1|Region/Stream),
             data = fungiSTAN_Ash, prior = prior1)
summary(V_AshN)
pp_check(V_AshN)

V_Ash15N<- brm(Delta15_N_Air_corr ~ Vegetation+(1|Region/Stream),
               data = fungiSTAN_Ash, prior = prior1)
summary(V_Ash15N)
pp_check(V_Ash15N)


### 2 Shredders ####
#### BRMS
# C
Crusta <- subset(dat.all, Order=="Crustacea" & T=="sample")
STAN_Crusta <-  Crusta %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region, Vegetation, EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr, Delta13_C_VPDB_corr) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr,Delta13_C_VPDB_corr), scale2) # scale numberic variables to 0

V_Crusta<- brm(EA_C_Percent_corr ~ Vegetation+(1|Region/Stream),
               data = STAN_Crusta, prior = prior1)
summary(V_Crusta)
pp_check(V_Crusta)


V_CrustaN<- brm(EA_N_Percent_corr ~ Vegetation+(1|Region/Stream),
                data = STAN_Crusta, prior = prior1)
summary(V_CrustaN)
pp_check(V_CrustaN)

V_Crusta15N<- brm(Delta15_N_Air_corr ~ Vegetation+(1|Region/Stream),
                  data = STAN_Crusta, prior = prior1)
summary(V_Crusta15N)
pp_check(V_Crusta15N)

# P
Plec <- subset(dat.all, Order=="Plecoptera" & T=="sample")
Plec1 <- Plec %>%
  dplyr::group_by(Sample_ID)%>%
  dplyr::mutate(EA_C_Percent_corr_mean = mean(EA_C_Percent_corr),
                EA_N_Percent_corr_mean = mean(EA_N_Percent_corr))
Plec1 <- Plec1 %>%
  distinct(Sample_ID, .keep_all=T) #155 bags with Plecoptera

STAN_Plec <-  Plec1 %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region, Vegetation, EA_C_Percent_corr_mean,EA_N_Percent_corr_mean,Delta15_N_Air_corr, Delta13_C_VPDB_corr)

STAN_Plec$EA_C_Percent_corr_mean<- scale2(STAN_Plec$EA_C_Percent_corr_mean)
hist(STAN_Plec$EA_C_Percent_corr_mean)
STAN_Plec$EA_N_Percent_corr_mean<- scale2(STAN_Plec$EA_N_Percent_corr_mean)
hist(STAN_Plec$EA_N_Percent_corr_mean)

STAN_Plec$Delta15_N_Air_corr<- scale2(STAN_Plec$Delta15_N_Air_corr)
hist(STAN_Plec$Delta15_N_Air_corr)

STAN_Plec$Delta13_C_VPDB_corr<- scale2(STAN_Plec$Delta13_C_VPDB_corr)
hist(STAN_Plec$Delta13_C_VPDB_corr)


V_Plec<- brm(EA_C_Percent_corr_mean ~ Vegetation+(1|Region/Stream),
             data = STAN_Plec, prior = prior1)
summary(V_Plec)
pp_check(V_Plec)


V_PlecN<- brm(EA_N_Percent_corr_mean ~ Vegetation+(1|Region/Stream),
              data = STAN_Plec, prior = prior1)
summary(V_PlecN)
pp_check(V_PlecN)


V_Plec15N<- brm(Delta15_N_Air_corr ~ Vegetation+(1|Region/Stream),
                data = STAN_Plec, prior = prior1)
summary(V_Plec15N)
pp_check(V_Plec15N)


# T
Tric <- subset(dat.all, Order=="Trichoptera" & T=="sample")
Tric1 <- Tric %>%
  dplyr::group_by(Sample_ID)%>%
  dplyr::mutate(EA_C_Percent_corr_mean = mean(EA_C_Percent_corr),
                EA_N_Percent_corr_mean = mean(EA_N_Percent_corr))
Tric1 <- Tric1 %>%
  distinct(Sample_ID, .keep_all=T)

STAN_Tric <-  Tric1 %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region, Vegetation, EA_C_Percent_corr_mean,EA_N_Percent_corr_mean,Delta15_N_Air_corr, Delta13_C_VPDB_corr,Delta13_C_VPDB_corr_corr )

#_corr_corr is corrected for Ethanol preservation


STAN_Tric$EA_C_Percent_corr_mean<- scale2(STAN_Tric$EA_C_Percent_corr_mean)
hist(STAN_Tric$EA_C_Percent_corr_mean)
STAN_Tric$EA_N_Percent_corr_mean<- scale2(STAN_Tric$EA_N_Percent_corr_mean)
hist(STAN_Tric$EA_N_Percent_corr_mean)
STAN_Tric$Delta15_N_Air_corr<- scale2(STAN_Tric$Delta15_N_Air_corr)
hist(STAN_Tric$Delta15_N_Air_corr)
STAN_Tric$Delta13_C_VPDB_corr<- scale2(STAN_Tric$Delta13_C_VPDB_corr)
hist(STAN_Tric$Delta13_C_VPDB_corr)
STAN_Tric$Delta13_C_VPDB_corr_corr<- scale2(STAN_Tric$Delta13_C_VPDB_corr)
hist(STAN_Tric$Delta13_C_VPDB_corr_corr)

V_Tric<- brm(EA_C_Percent_corr_mean ~ Vegetation+(1|Region/Stream),
             data = STAN_Tric, prior = prior1)
summary(V_Tric)
pp_check(V_Tric)


V_TricN<- brm(EA_N_Percent_corr_mean ~ Vegetation+(1|Region/Stream),
              data = STAN_Tric, prior = prior1)
summary(V_TricN)
pp_check(V_TricN)



V_Tric15N<- brm(Delta15_N_Air_corr ~ Vegetation+(1|Region/Stream),
                data = STAN_Tric, prior = prior1)
summary(V_Tric15N)
pp_check(V_Tric15N)


#### Figure SI ####
Plec1$EA_C_Percent_corr<- Plec1$EA_C_Percent_corr_mean ## rewrite this column
Plec1$EA_C_Percent_corr_mean<- NULL ## delete this column
Plec1$EA_N_Percent_corr<- Plec1$EA_N_Percent_corr_mean ## rewrite this column
Plec1$EA_N_Percent_corr_mean<- NULL ## delete this column

Tric1$EA_C_Percent_corr<- Tric1$EA_C_Percent_corr_mean ## rewrite this column
Tric1$EA_C_Percent_corr_mean<- NULL ## delete this column
Tric1$EA_N_Percent_corr<- Tric1$EA_N_Percent_corr_mean ## rewrite this column
Tric1$EA_N_Percent_corr_mean<- NULL ## delete this column

all<- rbind(Alder, Ash, Crusta, Plec1, Tric1) # 136+137+47+68+88=476
all <- all%>%
  mutate(Order=ifelse(LeafLeafLeaf=="Alder", "Alder",
                      ifelse(LeafLeafLeaf=="Ash", "Ash", Order)))


all$Combination<- as.factor(paste(all$Order,all$Vegetation, sep=""))
levels(all$Combination)
all$Combination<- factor(all$Combination, levels=c("AlderT0",
                                                   "Alderforested",
                                                   "Aldernon-forested",
                                                   "AshT0",
                                                   "Ashforested",
                                                   "Ashnon-forested",
                                                   "Crustaceaforested",
                                                   "Crustaceanon-forested",
                                                   "Plecopteraforested",
                                                   "Plecopteranon-forested",
                                                   "Trichopteraforested",
                                                   "Trichopteranon-forested"))

Shred2 <- c("grey18", "red", "red", "grey55", "darkorange", "darkorange",
            "#762a83", "#762a83", "#5ab4ac", "#5ab4ac", "#addd8e", "#addd8e")

bw=c("grey18", "white", "black")

c_perc<-ggplot(all)+theme_bw()+
  geom_boxplot(aes(x=Order, y=EA_C_Percent_corr, fill=Combination, color=Combination, alpha=Vegetation))+
  scale_fill_manual(values=Shred2)+
  scale_color_manual(values=Shred2)+
  scale_alpha_manual(values=c(0.2,0.4,0.7))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.7, 0.7), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15))+
  ylab(expression(paste("%C")))+
  annotate("text", y = 50, x= 0.825, label = "*", size=8)+
  annotate("text", y = 50, x= 1.125, label = "*", size=8)+
  theme(axis.text.x = element_blank())

c_perc







n_perc<-ggplot(all)+theme_bw()+
  geom_boxplot(aes(x=Order, y=EA_N_Percent_corr, fill=Combination, color=Combination, alpha=Vegetation))+
  scale_fill_manual(values=Shred2)+
  scale_color_manual(values=Shred2)+
  scale_alpha_manual(values=c(0.2,0.4,0.7))+
  guides(alpha= guide_legend(override.aes = list(fill= bw )))+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.2, 0.7), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15))+
  ylab(expression(paste("%N")))+
  annotate("text", y = 3.35, x= 0.875, label = "*", size=8)+
  annotate("text", y = 3.35, x= 2.125, label = "*", size=8)+
  annotate("text", y = 3.35, x= 1.875, label = "*", size=8)+
  theme(axis.text.x = element_blank())
n_perc

n_15<-ggplot(all)+theme_bw()+
  geom_boxplot(aes(x=Order, y=Delta15_N_Air_corr, fill=Combination, color=Combination, alpha=Vegetation))+
  scale_fill_manual(values=Shred2)+
  scale_color_manual(values=Shred2)+
  scale_alpha_manual(values=c(0.2,0.4,0.7))+
  guides(alpha= FALSE, size = FALSE)+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.9), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15))+
  ylab(expression(paste(delta^{15}, "N [\u2030]")))+
  annotate("text", y = 4.5, x= 3, label = "*", size=8)+
  annotate("text", y = 4.5, x= 2.125, label = "*", size=8)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
n_15


(c_perc/ n_perc/n_15 )+ plot_annotation(tag_levels = 'A')
ggsave("SI_plot.png", width = 14, height = 8)


##3. Mixing effects on shredders ############################

#### 3.1 shredder abundance ~ leaf treatment ####
sdat<- subset(dat.tot, Genus=="Capnia" | Genus=="Leuctra" | Genus== "Nemoura" | Genus=="Protonemura" | Species=="Gruppe AURICOLLIS" | Genus=="Limnephilini" | Species=="fossarum")

sdat <- sdat %>%
  dplyr::group_by(Site, Treatment)%>%
  dplyr::summarise(mean = mean(Total), sd = sd(Total)) ## 3*16=48 but one treatment not enough datapoints

sdat <- sdat %>%
  mutate(Region = ifelse(str_detect(Site, "G"), "North", "South"))

sdatSTAN <-  sdat %>%
  dplyr::mutate_at(vars(mean), scale2) # scale numberic variables to 0

sa<- brm(mean ~ Treatment+(1|Region),
         data = sdatSTAN, prior = prior1)
summary(sa)
pp_check(sa)


#### 3.2 shredder dry weight ~ leaf treatment ####
ddat<- subset(dat.tot, Genus=="Capnia" | Genus=="Leuctra" | Genus== "Nemoura" | Genus=="Protonemura" | Species=="Gruppe AURICOLLIS" | Genus=="Limnephilini" | Species=="fossarum")

ddat <- ddat %>%
  dplyr::group_by(Site, Treatment)%>%
  dplyr::summarise(mean = mean(DryMass,na.rm=TRUE), sd = sd(DryMass, na.rm=TRUE)) ## 3*16=48 but one treatment not enough datapoints

ddat <- ddat %>%
  mutate(Region = ifelse(str_detect(Site, "G"), "North", "South"))

ddatSTAN <-  ddat %>%
  dplyr::mutate_at(vars(mean), scale2) # scale numberic variables to 0

da<- brm(mean ~ Treatment+(1|Region),
         data = ddatSTAN, prior = prior1)
summary(da)
pp_check(da)

#### 3.3 shredder biomass between leaf treatment ####
bdat<- subset(dat.tot, Genus=="Capnia" | Genus=="Leuctra" | Genus== "Nemoura" | Genus=="Protonemura" | Species=="Gruppe AURICOLLIS" | Genus=="Limnephilini" | Species=="fossarum")

bdat <- bdat %>%
  dplyr::group_by(Site, Treatment)%>%
  dplyr::summarise(sum = sum(TotalDryMass,na.rm=TRUE), sd = sd(TotalDryMass, na.rm=TRUE)) ## 3*16=48 but one treatment not enough datapoints

bdat <- bdat %>%
  mutate(Region = ifelse(str_detect(Site, "G"), "North", "South"))

bdatSTAN <-  bdat %>%
  dplyr::mutate_at(vars(sum), scale2) # scale numberic variables to 0

ba<- brm(sum ~ Treatment+(1|Region),
         data = bdatSTAN, prior = prior1)
summary(ba)
pp_check(ba)




## 4. Ethanol Preservation ##############
Shreds <- c( "#762a83", "#5ab4ac",  "#addd8e") ## Crustacea, Plecoptera, Trichoptera

#### 4.1. Effect  ####
### BRMS

# C
CrustaT0 <- subset(dat.all, Order=="Crustacea")
STAN_Crusta <-  CrustaT0 %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region, T, Vegetation, EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr, Delta13_C_VPDB_corr) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(EA_C_Percent_corr,EA_N_Percent_corr,Delta15_N_Air_corr,Delta13_C_VPDB_corr), scale2) # scale numberic variables to 0


T0_Crusta<- brm(Delta13_C_VPDB_corr ~ T,
                data = STAN_Crusta, prior = prior1)
summary(T0_Crusta)
pp_check(T0_Crusta)




##P
PlecT0 <- subset(dat.all, Order=="Plecoptera")

STAN_PlecT0 <-  PlecT0 %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, T,Region, Vegetation, Delta15_N_Air_corr, Delta13_C_VPDB_corr)


STAN_PlecT0$Delta15_N_Air_corr<- scale2(STAN_PlecT0$Delta15_N_Air_corr)
STAN_PlecT0$Delta13_C_VPDB_corr<- scale2(STAN_PlecT0$Delta13_C_VPDB_corr)


V_PlecT0<- brm(Delta13_C_VPDB_corr ~ T,
               data = STAN_PlecT0, prior = prior1)
summary(V_PlecT0)
pp_check(V_PlecT0)

##P
TricT0 <- subset(dat.all, Order=="Trichoptera")

STAN_TricT0 <-  TricT0 %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, T,Region, Vegetation, Delta15_N_Air_corr, Delta13_C_VPDB_corr)

STAN_TricT0$Delta15_N_Air_corr<- scale2(STAN_TricT0$Delta15_N_Air_corr)
STAN_TricT0$Delta13_C_VPDB_corr<- scale2(STAN_TricT0$Delta13_C_VPDB_corr)

V_TricT0<- brm(Delta13_C_VPDB_corr ~ T,
               data = STAN_TricT0, prior = prior1)
summary(V_TricT0)
pp_check(V_TricT0)


### All
all <- subset(dat.all, Order=="Crustacea" | Order=="Plecoptera" | Order=="Trichoptera")

STAN_all <-  all %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, T,Order,Region, Vegetation, Delta15_N_Air_corr, Delta13_C_VPDB_corr)

STAN_all$Delta15_N_Air_corr<- scale2(STAN_all$Delta15_N_Air_corr)
STAN_all$Delta13_C_VPDB_corr<- scale2(STAN_all$Delta13_C_VPDB_corr)

V_PlecT0<- brm(Delta13_C_VPDB_corr ~ T,
               data = STAN_all, prior = prior1)
summary(V_PlecT0)
pp_check(V_PlecT0)


V_PlecT0_1<- brm(Delta13_C_VPDB_corr ~ T*Order,
                 data = STAN_all, prior = prior1)
summary(V_PlecT0_1)
pp_check(V_PlecT0_1)


#### Figure SI ####
ggplot(data=subset(dat.all, Shredder=="Y"))+
  geom_boxplot(aes(y=Delta13_C_VPDB_corr, x=T, fill=Order, color=Order), alpha=0.8)+
  theme_bw()+
  theme(axis.title.x=element_text(size=0),
        axis.text=element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15))+
  theme(legend.position = "bottom")+
  ylab(expression(paste(delta^{13}, "C [\u2030]")))+
  scale_fill_manual(values=Shreds)+
  scale_color_manual(values=Shreds)

ggsave("SI_plot_Conservation.png", width = 6, height = 6)


