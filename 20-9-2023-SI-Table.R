########################## Table 1 / Figure 1
############################### April 2023 
############################ Rebecca Oester

rm(list=ls())

setwd("C:/Users/Rebecca/switchdrive/Private/PhD/Data/Stable Isotopes")
###1. load packages and data ####
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

### Data
dat.all <- read.delim("C:/Users/rebec/switchdrive/Private/PhD/Data/dat.all.txt") 
dat.tot <- read.delim("C:/Users/rebec/switchdrive/Private/PhD/Data/dat.tot.txt")


### 2. Crustacea ####
gam <- subset(dat.tot, Taxon=="Gammarus fossarum")
gam <- gam[!duplicated(gam$Sample_ID),] ## delete duplicate Sample_ID ==> 47

dat.all.gamm <- right_join(dat.all, gam, by="Sample_ID")
dat.all.gamm <- subset(dat.all.gamm, Order.x=="Crustacea")

#### Table S1 ####
aggregate(Total ~Landscape, dat.all.gamm, mean)
aggregate(Total ~Landscape, dat.all.gamm, sd)
aggregate(DryMass ~Landscape, dat.all.gamm, mean)
aggregate(DryMass ~Landscape, dat.all.gamm, sd)


### 3. Plecoptera ####
#### 3.1 Capnia ####
cap <- subset(dat.tot, Genus=="Capnia")
cap$Species_All_1 <- cap$Genus
cap$Landscape<- factor(cap$Landscape, levels=(c("F", "A")))
dat.all.cap <- right_join(dat.all, cap, by=c("Species_All_1", "Sample_ID"))
dat.all.cap <- subset(dat.all.cap, !is.na(Sample_Name))

#### Table S1 ####
aggregate(Total ~Landscape, dat.all.cap, mean)
aggregate(Total ~Landscape, dat.all.cap, sd)
aggregate(DryMass ~Landscape, dat.all.cap, mean)
aggregate(DryMass ~Landscape, dat.all.cap, sd)


#### 3.1 Leuctra ####
leu <- subset(dat.tot, Genus=="Leuctra")
leu$Species_All_1 <- leu$Genus
leu$Landscape<- factor(leu$Landscape, levels=(c("F", "A")))
dat.all.leu <- right_join(dat.all, leu, by=c("Species_All_1", "Sample_ID"))
dat.all.leu <- subset(dat.all.leu, !is.na(Sample_Name))

#### Table S1 ####
aggregate(Total ~Landscape, dat.all.leu, mean)
aggregate(Total ~Landscape, dat.all.leu, sd)
aggregate(DryMass ~Landscape, dat.all.leu, mean)
aggregate(DryMass ~Landscape, dat.all.leu, sd)



#### 3.3 Nemoura ####
nem <- subset(dat.tot, Genus=="Nemoura")
nem$Species_All_1 <- nem$Genus
nem$Landscape<- factor(nem$Landscape, levels=(c("F", "A")))
dat.all.nem <- right_join(dat.all, nem, by=c("Species_All_1", "Sample_ID"))
dat.all.nem <- subset(dat.all.nem, !is.na(Sample_Name))

#### Table S1 ####
aggregate(Total ~Landscape, dat.all.nem, mean)
aggregate(Total ~Landscape, dat.all.nem, sd)
aggregate(DryMass ~Landscape, dat.all.nem, mean)
aggregate(DryMass ~Landscape, dat.all.nem, sd)


#### 3.4 Protoproura ####
pro <- subset(dat.tot, Genus=="Protonemura")
pro$Species_All_1 <- pro$Genus
pro$Landscape<- factor(pro$Landscape, levels=(c("F", "A")))
dat.all.pro <- right_join(dat.all, pro, by=c("Species_All_1", "Sample_ID"))
dat.all.pro <- subset(dat.all.pro, !is.na(Sample_Name))

#### Table S1 ####
aggregate(Total ~Landscape, dat.all.pro, mean)
aggregate(Total ~Landscape, dat.all.pro, sd)
aggregate(DryMass ~Landscape, dat.all.pro, mean)
aggregate(DryMass ~Landscape, dat.all.pro, sd)


### 4. Trichoptera ####
#### 4.1 Limnephilini ####
lim <- subset(dat.tot, Genus=="Limnephilini" )
lim$Species_All_1 <- lim$Genus
lim$Landscape<- factor(lim$Landscape, levels=(c("F", "A")))
dat.all.lim <- right_join(dat.all, lim, by=c("Species_All_1", "Sample_ID"))
dat.all.lim <- subset(dat.all.lim, !is.na(Sample_Name))

#### Table S1 ####
aggregate(Total ~Landscape, dat.all.lim, mean)
aggregate(Total ~Landscape, dat.all.lim, sd)
aggregate(DryMass ~Landscape, dat.all.lim, mean)
aggregate(DryMass ~Landscape, dat.all.lim, sd)

#### 4.2 Auricollis ####
aur <- subset(dat.tot, Species=="Gruppe AURICOLLIS" )
aur$Species_All_1 <- aur$Taxon
aur$Species_All_1<- gsub("Chaetopt.u.Stenoph. Gruppe AURICOLLIS", "Auricollis", aur$Species_All_1)
aur$Landscape<- factor(aur$Landscape, levels=(c("F", "A")))
dat.all.aur <- right_join(dat.all, aur, by=c("Species_All_1", "Sample_ID"))
dat.all.aur <- subset(dat.all.aur, !is.na(Sample_Name))

#### Table S1 ####
aggregate(Total ~Landscape, dat.all.aur, mean)
aggregate(Total ~Landscape, dat.all.aur, sd)
aggregate(DryMass ~Landscape, dat.all.aur, mean)
aggregate(DryMass ~Landscape, dat.all.aur, sd)
