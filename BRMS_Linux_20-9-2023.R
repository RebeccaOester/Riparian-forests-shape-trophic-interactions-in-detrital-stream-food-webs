##################################################### Detrital food webs ############################
##################################################### Direct effects ################################
#################################################### Rebecca Oester #################################
#################################################### September, 2023 ################################

### Data available on Dryad-Link:
rm(list=ls())
#### 1. Packages ####
library(rstan)
library(brms)
library(bayesplot)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(outliers)
library(stringr)
library(patchwork)

#### scaling function to standardize all numeric variables
scale2 <- function(x, na.rm = TRUE) ((x - mean(x, na.rm = na.rm)) / sd(x, na.rm))/2 ### na.rm = TRUE leaves NAs in the scaling and just considers values

#### prior setting
#### vaguely informative prior with normal distribution around 0
prior1 <- c(prior(normal(0, 1), class = Intercept),
            prior(normal(0, 0.1), class = "b"))

#### 2. Fungi data ####
## set working directory
setwd("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes")

# Data
dat.erg<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.erg.txt")
dat.erg <- subset(dat.erg, T=="sample")

### prepare data for Resources ####
## Alder
Alder <- subset(dat.erg, LeafLeafLeaf=="Alder")
fungiSTAN_Alder <-  Alder %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region, LeafLeafLeaf, Vegetation,Delta13_C_VPDB,log_molar_CN,Fungal_Biom_Leaf_mgg_noRR) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Delta13_C_VPDB,log_molar_CN,Fungal_Biom_Leaf_mgg_noRR), scale2) # scale numberic variables to 0

hist(fungiSTAN_Alder$Delta13_C_VPDB)
hist(fungiSTAN_Alder$log_molar_CN)
hist(fungiSTAN_Alder$Fungal_Biom_Leaf_mgg_noRR)

## Ash
Ash <- subset(dat.erg, LeafLeafLeaf=="Ash")
fungiSTAN_Ash <-  Ash %>%
  dplyr::select(Sample_Name,Sample_ID,Stream, Region, LeafLeafLeaf, Vegetation,Delta13_C_VPDB,log_molar_CN,Fungal_Biom_Leaf_mgg_noRR) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Delta13_C_VPDB,log_molar_CN,Fungal_Biom_Leaf_mgg_noRR), scale2) # scale numberic variables to 0

hist(fungiSTAN_Ash$Delta13_C_VPDB)
hist(fungiSTAN_Ash$log_molar_CN)
hist(fungiSTAN_Ash$Fungal_Biom_Leaf_mgg_noRR)


## simple models Alder
funmod2<- brm(Fungal_Biom_Leaf_mgg_noRR ~ Vegetation+(1|Region/Stream),
              data = fungiSTAN_Alder, prior = prior1)
summary(funmod2)
pp_check(funmod2)
plot(funmod2)

funmod2P <- brm(log_molar_CN ~ Fungal_Biom_Leaf_mgg_noRR+Vegetation+(1|Region/Stream),
                data = fungiSTAN_Alder,
                prior = prior1)

summary(funmod2P)
pp_check(funmod2P)
plot(funmod2P)

funmod3 <- brm(Delta13_C_VPDB ~ Fungal_Biom_Leaf_mgg_noRR+Vegetation+(1|Region/Stream),
               data = fungiSTAN_Alder,
               prior = prior1)

summary(funmod3)
pp_check(funmod3)
plot(funmod3)

## Prepare SEM
funmod1.1 <- brmsformula(Fungal_Biom_Leaf_mgg_noRR ~ Vegetation+(1|Region/Stream))
funmod2.1 <- brmsformula(log_molar_CN ~ Fungal_Biom_Leaf_mgg_noRR+Vegetation+(1|Region/Stream))
funmod3.1 <- brmsformula(Delta13_C_VPDB ~ Fungal_Biom_Leaf_mgg_noRR+Vegetation+(1|Region/Stream))

## SEM
Stan_Alder_priors_4 <- brm(
  funmod1.1+
    funmod2.1+
    funmod3.1,
  chains=4, cores=4,iter=10000,
  family=c(gaussian,
           gaussian,
           gaussian),
  prior= prior1,
  inits=0,data = fungiSTAN_Alder,
  save_all_pars=T,control = list(adapt_delta = 0.99,max_treedepth=20))

Stan_Alder_priors_4$prior
summary(Stan_Alder_priors_4)

pp_check(Stan_Alder_priors_4, resp = 'FungalBiomLeafmggnoRR') ### looks ok
pp_check(Stan_Alder_priors_4, resp = 'logmolarCN') ### looks ok
pp_check(Stan_Alder_priors_4, resp = 'Delta13CVPDB') ### looks ok

## simple models Ash
funmod1_as <- brm(Fungal_Biom_Leaf_mgg_noRR ~ Vegetation+(1|Region/Stream),
                  data = fungiSTAN_Ash, prior=prior1)
summary(funmod1_as)
pp_check(funmod1_as)
plot(funmod1_as)

funmod2_as <- brm(log_molar_CN ~ Fungal_Biom_Leaf_mgg_noRR+Vegetation+(1|Region/Stream),
                  data = fungiSTAN_Ash, prior=prior1)
summary(funmod2_as)
pp_check(funmod2_as)
plot(funmod2_as)

funmod3_as <- brm(Delta13_C_VPDB ~ Fungal_Biom_Leaf_mgg_noRR+Vegetation+(1|Region/Stream),
                  data = fungiSTAN_Ash, prior=prior1)
summary(funmod3_as)
pp_check(funmod3_as)
plot(funmod3_as)


## prepare SEM
funmod1.1 <- brmsformula(Fungal_Biom_Leaf_mgg_noRR ~ Vegetation+(1|Region/Stream))
funmod2.1 <- brmsformula(log_molar_CN ~ Fungal_Biom_Leaf_mgg_noRR+Vegetation+(1|Region/Stream))
funmod3.1 <- brmsformula(Delta13_C_VPDB ~ Fungal_Biom_Leaf_mgg_noRR+Vegetation+(1|Region/Stream))


## SEM
Stan_Ash <- brm(
  funmod1.1+
    funmod2.1+
    funmod3.1,
  chains=4, cores=4,iter=10000,
  family=c(gaussian,
           gaussian,
           gaussian),
  prior= prior1,
  inits=0,data = subset(fungiSTAN, LeafLeafLeaf=="Ash"),
  save_all_pars=T,control = list(adapt_delta = 0.99,max_treedepth=20))


summary(Stan_Ash)
Stan_Ash$prior

pp_check(Stan_Ash, resp = 'FungalBiomLeafmggnoRR') ### looks ok
pp_check(Stan_Ash, resp = 'logmolarCN') ### looks ok
pp_check(Stan_Ash, resp = 'Delta13CVPDB') ### looks ok



### T0 check ####
dat.erg<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.erg.txt")


T0_Alder <- subset(dat.erg, Species=="Alder")
T0_Alder_T0 <- subset(T0_Alder, T=="T0")
mean(T0_Alder_T0$log_molar_CN, na.rm = T) # 3.271886
sd(T0_Alder_T0$log_molar_CN, na.rm = T) # 0.06394871
T0_Alder_Tx <- subset(T0_Alder, T!="T0")
mean(T0_Alder_Tx$log_molar_CN, na.rm = T) # 2.938827
sd(T0_Alder_Tx$log_molar_CN, na.rm = T) # 0.07844415

T0_Alder <- T0_Alder %>%
  dplyr::select(Sample_Name, T, log_molar_CN, Delta13_C_VPDB) %>%
  dplyr::mutate_at(vars(Delta13_C_VPDB,log_molar_CN), scale2)

T0_AlderCN <- brm(log_molar_CN ~ T, data = T0_Alder, prior=prior1)
summary(T0_AlderCN)
T0_AlderC <- brm(Delta13_C_VPDB ~ T, data = T0_Alder, prior=prior1)
summary(T0_AlderC)



T0_Ash <- subset(dat.erg, Species=="Ash")
T0_Ash_T0 <- subset(T0_Ash, T=="T0")
mean(T0_Ash_T0$log_molar_CN, na.rm = T) # 3.473431
sd(T0_Ash_T0$log_molar_CN, na.rm = T) # 0.06154893
T0_Ash_Tx <- subset(T0_Ash, T!="T0")
mean(T0_Ash_Tx$log_molar_CN, na.rm = T) # 2.983141
sd(T0_Ash_Tx$log_molar_CN, na.rm = T) # 0.08511103

T0_Ash <- T0_Ash %>%
  dplyr::select(Sample_Name, T, log_molar_CN, Delta13_C_VPDB) %>%
  dplyr::mutate_at(vars(Delta13_C_VPDB,log_molar_CN), scale2)

T0_AshCN <- brm(log_molar_CN ~ T, data = T0_Ash, prior=prior1)
summary(T0_AshCN)
T0_AshC <- brm(Delta13_C_VPDB ~ T, data = T0_Ash, prior=prior1)
summary(T0_AshC)



#### Shredder Data ####
dat.all<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.all.txt")
dat.tot<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.tot.txt")


dat.shred <- subset(dat.all, Shredder=="Y" & T=="sample")
mean(dat.shred$log_molar_CN, na.rm=T) #1.593865
aggregate(log_molar_CN ~ Order, mean, data=dat.shred)
# Order log_molar_CN
# 1   Crustacea     1.556847
# 2  Plecoptera     1.465633
# 3 Trichoptera     1.737009
aggregate(log_molar_CN ~ Order, sd, data=dat.shred)
# Order log_molar_CN
# 1   Crustacea   0.03864487
# 2  Plecoptera   0.04687921
# 3 Trichoptera   0.12281525

## quick comparison between T0-Tx and Shred C:N
3.271886/1.593865 #AlderT0/shredder
3.473431/1.593865 #AlderTx/shredder
2.938827/1.593865 #AshT0/shredder
2.0528/1.593865 #AshTx/shredder


##### 3. Crustacea ####
gam <- subset(dat.tot, Taxon=="Gammarus fossarum")
gam <- gam[!duplicated(gam$Sample_ID),] ## delete duplicate Sample_ID ==> 47

dat.all.gamm <- right_join(dat.all, gam, by="Sample_ID")
dat.all.gamm <- subset(dat.all.gamm, Order.x=="Crustacea")

### prepare data
crustaSTAN <-  dat.all.gamm %>%
  dplyr::select(Sample_Name,Sample_ID,Stream.x, Region.x, LeafLeafLeaf, Vegetation,Delta13_C_VPDB_corr,log_molar_CN, DryMass, TotalDryMass) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Delta13_C_VPDB_corr,log_molar_CN,DryMass, TotalDryMass), scale2) # scale numberic variables to 0


crustaSTAN <-  crustaSTAN[complete.cases(crustaSTAN), ]

### model for Crustacea
crumod1 <- brm(DryMass ~ Vegetation+(1|Stream.x), data = crustaSTAN, prior=prior1)
summary(crumod1)
crumod1B <- brm(TotalDryMass ~ Vegetation+(1|Stream.x), data = crustaSTAN, prior=prior1)
summary(crumod1B)

crumod2 <- brm(log_molar_CN ~ DryMass+Vegetation+(1|Stream.x), data = crustaSTAN, prior=prior1)
summary(crumod2)
crumod3 <- brm(Delta13_C_VPDB_corr ~ DryMass+Vegetation+(1|Stream.x), data = crustaSTAN,  prior=prior1)
summary(crumod3)

crumod1.1 <- brmsformula(DryMass ~ Vegetation+(1|Stream.x))
crumod2.1 <- brmsformula(log_molar_CN ~ DryMass+Vegetation+(1|Stream.x))
crumod3.1 <- brmsformula(Delta13_C_VPDB_corr ~ DryMass+Vegetation+(1|Stream.x))


Stan_Cru <- brm(
  crumod1.1+
    crumod2.1+
    crumod3.1,
  chains=4, cores=4,iter=10000,
  family=c(gaussian,
           gaussian,
           gaussian),
  prior= prior1,
  inits=0,data = crustaSTAN,
  save_all_pars=T,control = list(adapt_delta = 0.99,max_treedepth=20))


summary(Stan_Cru)
prior_summary(Stan_Cru)
pp_check(Stan_Cru, resp = 'DryMass') ### looks ok
pp_check(Stan_Cru, resp = 'logmolarCN') ### looks ok
pp_check(Stan_Cru, resp = 'Delta13_C_VPDB_corr') ### looks ok


##### 4. Plecoptera ####
plec <- subset(dat.tot, Genus=="Capnia" | Genus=="Leuctra" | Genus== "Nemoura" | Genus=="Protonemura")
plec$Species_All_1 <- plec$Genus
plec$Landscape<- factor(plec$Landscape, levels=(c("F", "A")))

plec <- plec %>%
  dplyr::group_by(Sample_ID)%>%
  dplyr::mutate(TotalDryMass_Plec = sum(TotalDryMass),
                MeanDryMass_Plec = mean(DryMass))
plec <- plec %>%
  distinct(Sample_ID, .keep_all=T) #155 bags with Plecoptera

dat.all.plec <- left_join(dat.all, plec, by=c("Sample_ID" ))
dat.all.plec <- subset(dat.all.plec, Order.x=="Plecoptera" & T=="sample")

dat.all.plec <- subset(dat.all.plec, !is.na(TotalDryMass_Plec))

dat.all.plec <- dat.all.plec %>%
  distinct(Sample_ID, .keep_all=T) ## 67 bags with Plecoptera

dat.all.plec$TotalDryMass_Plec
dat.all.plec$MeanDryMass_Plec


### prepare data
plecSTAN <-  dat.all.plec %>%
  dplyr::select(Sample_Name,Sample_ID,Stream.x, Region.x, Vegetation,Delta13_C_VPDB_corr,log_molar_CN, TotalDryMass_Plec, MeanDryMass_Plec) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Delta13_C_VPDB_corr,log_molar_CN,TotalDryMass_Plec, MeanDryMass_Plec), scale2) # scale numberic variables to 0

plecSTAN <-  plecSTAN[complete.cases(plecSTAN), ]


### model for Plecoptera
plemod1 <- brm(MeanDryMass_Plec ~ Vegetation+(1|Region.x/Stream.x) , data = plecSTAN, prior=prior1)
summary(plemod1)
plemod1B <- brm(TotalDryMass_Plec ~ Vegetation+(1|Region.x/Stream.x) , data = plecSTAN, prior=prior1)
summary(plemod1B)
plemod2 <- brm(log_molar_CN ~ MeanDryMass_Plec+Vegetation+(1|Region.x/Stream.x), data = plecSTAN, prior=prior1)
summary(plemod2)
plemod3 <- brm(Delta13_C_VPDB_corr ~ MeanDryMass_Plec+Vegetation+(1|Region.x/Stream.x), data = plecSTAN,  prior=prior1)
summary(plemod3)

plemod1.1 <- brmsformula(MeanDryMass_Plec ~ Vegetation+(1|Region.x/Stream.x))
plemod2.1 <- brmsformula(log_molar_CN ~ MeanDryMass_Plec+Vegetation+(1|Region.x/Stream.x))
plemod3.1 <- brmsformula(Delta13_C_VPDB_corr ~ MeanDryMass_Plec+Vegetation+(1|Region.x/Stream.x))

Stan_Plec <- brm(
  plemod1.1+
    plemod2.1+
    plemod3.1,
  chains=4, cores=4,iter=10000,
  family=c(gaussian,
           gaussian,
           gaussian),
  prior= prior1,
  inits=0,data = plecSTAN,
  save_all_pars=T,control = list(adapt_delta = 0.99,max_treedepth=20))


summary(Stan_Plec)
pp_check(Stan_Plec, resp = 'DryMassy') ### looks ok
pp_check(Stan_Plec, resp = 'logmolarCN') ### looks ok
pp_check(Stan_Plec, resp = 'Delta13CVPDB') ### looks ok

##### 5. Trichoptera ####
tri <- subset(dat.tot, Species=="Gruppe AURICOLLIS" | Genus=="Limnephilini" )
tri$Species_All_1 <- tri$Taxon
tri$Species_All_1<- gsub("Chaetopt.u.Stenoph. Gruppe AURICOLLIS", "Auricollis", tri$Species_All_1)
tri$Landscape<- factor(tri$Landscape, levels=(c("F", "A")))

tri <- tri %>%
  dplyr::group_by(Sample_ID)%>%
  dplyr::mutate(TotalDryMass_Tric = sum(TotalDryMass),
                MeanDryMass_Tric = mean(DryMass))

tri <- tri %>%
  distinct(Sample_ID, .keep_all=T)


dat.all.tri <- left_join(dat.all, tri, by=c("Species_All_1", "Sample_ID"))
dat.all.tri <- subset(dat.all.tri, !is.na(Sample_Name))
dat.all.tri <- subset(dat.all.tri, Order.x=="Trichoptera" & T=="sample")

dat.all.tri <- subset(dat.all.tri, !is.na(TotalDryMass_Tric))

dat.all.tri <- dat.all.tri %>%
  distinct(Sample_ID, .keep_all=T)


triSTAN <-  dat.all.tri %>%
  dplyr::select(Sample_Name,Sample_ID,Stream.x, Region.x, LeafLeafLeaf, Vegetation,Delta13_C_VPDB_corr,log_molar_CN, TotalDryMass_Tric, MeanDryMass_Tric) %>% ## subset only the variables of interest
  dplyr::mutate_at(vars(Delta13_C_VPDB_corr,log_molar_CN,TotalDryMass_Tric, MeanDryMass_Tric), scale2) # scale numberic variables to 0

triSTAN <- triSTAN[complete.cases(triSTAN), ]

### model for Trichoptera
trimod1 <- brm(MeanDryMass_Tric ~ Vegetation+(1|Region.x/Stream.x) , data = triSTAN, prior=prior1)
summary(trimod1)
trimod1B <- brm(TotalDryMass_Tric ~ Vegetation+(1|Region.x/Stream.x) , data = triSTAN, prior=prior1)
summary(trimod1B)

trimod2 <- brm(log_molar_CN ~ MeanDryMass_Tric+Vegetation+(1|Region.x/Stream.x), data = triSTAN, prior=prior1)
summary(trimod2)
trimod3 <- brm(Delta13_C_VPDB_corr ~ MeanDryMass_Tric+Vegetation+(1|Region.x/Stream.x), data = triSTAN, prior=prior1)
summary(trimod3)


trimod1.1 <- brmsformula(MeanDryMass_Tric ~ Vegetation+(1|Region.x/Stream.x))
trimod2.1 <- brmsformula(log_molar_CN ~ MeanDryMass_Tric+Vegetation+(1|Region.x/Stream.x))
trimod3.1 <- brmsformula(Delta13_C_VPDB_corr ~ MeanDryMass_Tric+Vegetation+(1|Region.x/Stream.x))

Stan_Tri <- brm(
  trimod1.1+
    trimod2.1+
    trimod3.1,
  chains=4, cores=4,iter=10000,
  family=c(gaussian,
           gaussian,
           gaussian),
  prior= prior1,
  inits=0,data = triSTAN,
  save_all_pars=T,control = list(adapt_delta = 0.99,max_treedepth=20))


summary(Stan_Tri)
pp_check(Stan_Tri, resp = 'DryMassy') ### looks ok
pp_check(Stan_Tri, resp = 'logmolarCN') ### looks ok
pp_check(Stan_Tri, resp = 'Delta13CVPDB') ### looks ok


save.image("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/all_data.RData")



## Figure  ####
rm(list=ls())

##### Figure data ####
dat.all<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.all.txt")
dat.erg<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.erg.txt")
dat.tot<- read.delim("/media/rebecca/Windows/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes/dat.tot.txt")

dat.erg <- dat.erg %>%
  mutate(Vegetation = ifelse(Country=="T0", "T0", Vegetation))

Shreds <- c(  "red", "darkorange","#762a83", "#5ab4ac",  "#addd8e")

bw=c("white", "black")

## quick plot on fungal biomass
ggplot(data=subset(dat.erg, T=="sample"))+theme_bw()+
  geom_boxplot(aes(x= Vegetation, y= Fungal_Biom_Leaf_mgg_noRR, alpha=Vegetation), fill="darkorange4", color="darkorange4")+
  facet_wrap(~Stream)+
  scale_alpha_manual(values=c(0.4,0.7))+
  stat_summary(aes(x= Vegetation, y= Fungal_Biom_Leaf_mgg_noRR), fun.y=mean, geom="point", shape=20, size=3, color="black", fill="black")

aggregate(Fungal_Biom_Leaf_mgg_noRR~Stream+Vegetation, data=dat.erg, mean)

ggplot(data=subset(dat.erg, T=="sample"))+theme_bw()+
  geom_boxplot(aes(x= LeafLeafLeaf, y= Fungal_Biom_Leaf_mgg_noRR, alpha=Vegetation, fill=LeafLeafLeaf, color=LeafLeafLeaf))+
  stat_summary(aes(x= LeafLeafLeaf, y= Fungal_Biom_Leaf_mgg_noRR, alpha=Vegetation, fill=LeafLeafLeaf), fun.y=mean, geom="point", shape=20, size=3, color="black",
               position = position_dodge(width = 0.75))+
  facet_wrap(~Stream)+
  scale_alpha_manual(values=c(0.4,0.7))+
  scale_fill_manual(values=Shreds)+
  scale_color_manual(values=Shreds)+
  scale_alpha_manual(values=c(0.4,0.7))+
  guides(alpha= guide_legend(override.aes = list(fill= bw )))+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.8, 0.1), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15))+
  ylab(expression(paste("Fungal biomass [mg/g]")))

##### Shredders ####
###### 1. Crustacea ####
gam <- subset(dat.tot, Taxon=="Gammarus fossarum")
gam <- gam[!duplicated(gam$Sample_ID),] ## delete duplicate Sample_ID ==> 47
gam <- select(gam, Sample_ID, DryMass, TotalDryMass, Order, Landscape)
gam$Landscape<- factor(gam$Landscape, levels=(c("F", "A")))

dat.all.gamm <- right_join(dat.all, gam, by="Sample_ID")
dat.all.gamm <- subset(dat.all.gamm, Order.x=="Crustacea")

###### 2. Plecoptera ####
plec <- subset(dat.tot, Genus=="Capnia" | Genus=="Leuctra" | Genus== "Nemoura" | Genus=="Protonemura")
## 292 times we found these Plecoptera in the leaf litter bags
plec$Species_All_1 <- plec$Genus
plec <- select(plec, Sample_ID, Total, Order,TotalDryMass, DryMass, Species_All_1 , Landscape)

# sum up the biomass of all P shredders and keep 1 row for each bag
plec <- plec %>%
  dplyr::group_by(Sample_ID)%>%
  dplyr::mutate(TotalDryMass_Plec = sum(TotalDryMass),
                MeanDryMass_Plec = mean(DryMass))
plec <- plec %>%
  distinct(Sample_ID, .keep_all=T) #155 bags with Plecoptera

dat.all.plec <- left_join(dat.all, plec, by=c("Sample_ID" ))
dat.all.plec <- subset(dat.all.plec, Order.x=="Plecoptera" & T=="sample")

dat.all.plec <- subset(dat.all.plec, !is.na(TotalDryMass_Plec))

dat.all.plec <- dat.all.plec %>%
  distinct(Sample_ID, .keep_all=T) ## 67 bags with Plecoptera

###### 3. Trichoptera ####
tri <- subset(dat.tot, Species=="Gruppe AURICOLLIS" | Genus=="Limnephilini" )
tri$Species_All_1 <- tri$Taxon
tri$Species_All_1<- gsub("Chaetopt.u.Stenoph. Gruppe AURICOLLIS", "Auricollis", tri$Species_All_1)

tri <- select(tri, Sample_ID, Total, Order, TotalDryMass, DryMass, Species_All_1, Landscape )

tri <- tri %>%
  dplyr::group_by(Sample_ID)%>%
  dplyr::mutate(TotalDryMass_Tric = sum(TotalDryMass),
                MeanDryMass_Tric = mean(DryMass))

tri <- tri %>%
  distinct(Sample_ID, .keep_all=T)


dat.all.tri <- left_join(dat.all, tri, by=c("Species_All_1", "Sample_ID"))
dat.all.tri <- subset(dat.all.tri, !is.na(Sample_Name))
dat.all.tri <- subset(dat.all.tri, Order.x=="Trichoptera" & T=="sample")

dat.all.tri <- subset(dat.all.tri, !is.na(TotalDryMass_Tric))

dat.all.tri <- dat.all.tri %>%
  distinct(Sample_ID, .keep_all=T)


###### C:N and d13C panels ####
dat.erg.shortC <- select(dat.erg, Sample_Name, Vegetation, LeafLeafLeaf,Delta13_C_VPDB, log_molar_CN, Stream)
colnames(dat.erg.shortC)[colnames(dat.erg.shortC) == "LeafLeafLeaf"] ="Order"
colnames(dat.erg.shortC)[colnames(dat.erg.shortC) == "Delta13_C_VPDB"] ="Delta13_C_VPDB_corr"


dat.gam.shortC <- select(dat.all.gamm, Sample_Name,Vegetation,Order.y, Delta13_C_VPDB_corr, log_molar_CN, Stream)
colnames(dat.gam.shortC)[colnames(dat.gam.shortC) == "Order.y"] ="Order"

dat.plec.shortC <- select(dat.all.plec, Sample_Name, Vegetation,Order.y, Delta13_C_VPDB_corr, log_molar_CN, Stream)
colnames(dat.plec.shortC)[colnames(dat.plec.shortC) == "Order.y"] ="Order"

dat.tric.shortC <- select(dat.all.tri, Sample_Name, Vegetation,Order.y, Delta13_C_VPDB_corr, log_molar_CN, Stream)
colnames(dat.tric.shortC)[colnames(dat.tric.shortC) == "Order.y"] ="Order"


datC <- rbind(dat.erg.shortC, dat.gam.shortC,dat.plec.shortC , dat.tric.shortC)
datC$Order<- factor(datC$Order, levels=(c("Alder", "Ash", "Crustacea", "Plecoptera", "Trichoptera")))
datC$Vegetation<- factor(datC$Vegetation, levels=(c("T0", "forested", "non-forested")))


datC<- datC%>%
  mutate(Color = ifelse(Vegetation=="T0", "grey",
                        ifelse(Order=="Alder", "red",
                               ifelse(Order=="Ash", "darkorange",
                                      ifelse(Order=="Crustacea","#762a83",
                                             ifelse(Order=="Plecoptera", "#5ab4ac", "#addd8e" ))))))
datC<- datC%>%
  mutate(Alpha = ifelse(Vegetation=="T0", 0.1,
                        ifelse(Vegetation=="forested", 0.4, 0.7)))

datC$Combination<- as.factor(paste(datC$Order,datC$Vegetation, sep=""))
datC$Combination<- factor(datC$Combination, levels=c("AlderT0",
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
Shred3 <- c("#762a83",  "#5ab4ac", "#addd8e")

bw=c("grey18", "white", "black")

c<-ggplot(datC)+theme_bw()+
  geom_boxplot(aes(x=Order, y=Delta13_C_VPDB_corr, fill=Combination, color=Combination, alpha=Vegetation))+
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
  ylab(expression(paste(delta^{13}, "C [\u2030]")))+
  annotate("text", y = -25, x= 3, label = "*", size=8)+
  annotate("text", y = -25, x= 4, label = "*", size=8)+
  annotate("text", y = -28, x= 1.875, label = "*", size=8)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))+
  scale_y_continuous(limits=c(-35, -21), breaks=c(-35, -28,-21))
c

bw=c("grey18", "white", "black")
cn<-ggplot(datC)+theme_bw()+
  geom_boxplot(aes(x=Order, y=log_molar_CN, fill=Combination, color=Combination, alpha=Vegetation))+
  scale_fill_manual(values=Shred2)+
  scale_color_manual(values=Shred2)+
  scale_alpha_manual(values=c(0.2,0.4,0.7))+
  guides(alpha= guide_legend(override.aes = list(fill= bw )))+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(legend.position = c(0.7, 0.7), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.direction="vertical")+
  labs(alpha="Riparian vegetation:")+
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15))+
  ylab(expression(paste(C:N)))+
  annotate("text", y = 3.35, x= 0.875, label = "*", size=8)+
  annotate("text", y = 3.35, x= 1.125, label = "*", size=8)+
  annotate("text", y = 3.35, x= 1.875, label = "*", size=8)+
  annotate("text", y = 3.35, x= 2.125, label = "*", size=8)+
  annotate("text", y = 2.25, x= 5, label = "*", size=8)+
  theme(axis.text.x = element_blank())+
  scale_y_continuous(limits=c(1.25, 3.5), breaks=c(1.5, 2 ,2.5,3, 3.5 ))
cn


###### ellipse plot ####
bw=c( "white", "black")
all<- ggplot(subset(datC, Vegetation!="T0"))+theme_bw()+
  geom_point(aes(x=Delta13_C_VPDB_corr, y=log_molar_CN, color=Order, alpha=Vegetation), size=3)+
  stat_ellipse(aes(x=Delta13_C_VPDB_corr, y=log_molar_CN, color=Order, alpha=Vegetation))+
  scale_fill_manual(values=Shreds)+
  scale_color_manual(values=Shreds)+
  scale_alpha_manual(values=c(0.4,0.7))+
  guides(alpha= guide_legend(override.aes = list(fill= bw )))+
  guides(fill = FALSE, size = FALSE)+
  guides(color = FALSE, size = FALSE)+
  theme(axis.text=element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        plot.title = element_text(size=15))+
  ylab(expression(paste(C:N)))+
  theme(legend.position = "none")+
  xlab(expression(paste(delta^{13}, "C [\u2030]")))+facet_wrap(~Vegetation)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text.x = element_text(size = 15))+
  scale_y_continuous(limits=c(1.25, 3.5), breaks=c(1.5, 2 ,2.5,3, 3.5 ))+
  scale_x_continuous(limits=c(-35, -21), breaks=c(-35, -28,-21))
all


###### 3 panel plot ####
ggarrange(ggarrange(all, ncol = 1, labels = c("A")),
          ggarrange(cn,c, labels =c("B","C"), nrow = 2, heights = c(1, 1.25), align="v"))

ggsave("food-web-2.png", width = 14, height = 8)

