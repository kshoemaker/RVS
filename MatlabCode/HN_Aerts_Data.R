### Case Study Data Processing ###
  
# setup 
rm(list = ls())
library(MASS)
library(BAMMtools)
library(tidyverse)
library(gtools)


#set.seed(12345678) ## seed good for stages?s
set.seed(1234)
#setwd("~/Google Drive/Classes/CP_Project/LungProjectData")
#setwd("C:/Users/KAShoemaker/Box Sync/Stats + radiomics project")
setwd("~/Documents/GitHub/ReliableRadiomicsVS/CaseStudyData/")
# 78 subjects
hn_data <- read.csv("CTPatientDataForKate.csv", stringsAsFactors = T)
# hn_data <- hn_data[c(-81,-82),]

## load in reliability information 
reliability_data <- read.csv("LungProjectData/TotalSDForKate.csv", header = T, stringsAsFactors = F)
#reliability_data$Control.Total.SD <- log(reliability_data$Control.Total.SD)






hn_aerts <- filter(hn_data,DataType =="A") ## ,  HPV.status == "Negative")
# hn_aerts <- hn_aerts[c(-81,-82),]
# hn_aerts <- hn_data
hn_clinical <- hn_aerts[,1:14]
 hn_radiomics <- hn_aerts[,15:186]
# hn_radiomics <- hn_aerts[,58:100]
hn_radiomics <- hn_radiomics[,which(colnames(hn_radiomics) %in% reliability_data$PreFeat)]


# par(mfrow = c(1,2))
# for(i in 1:64){
#   hist(hn_radiomics[,i], main = paste("Histogram of Variable ",i))
#   if(min(hn_radiomics[,i] < 0)){
#     hist(log(hn_radiomics[,i] - min(hn_radiomics[,i]) + 0.1),main = paste("Histogram of Variable ",i))}
#   else {hist(log(hn_radiomics[,i]),main = paste("Histogram of Variable ",i)) }
# }

#hn_radiomics[,160] <- hn_radiomics[,160]*100000
library(car)
trans_hn_radiomics <- hn_radiomics
for (i in 1:(ncol(hn_radiomics))){
  vec <- hn_radiomics[,i]
  if (min(vec) <= 0){
    fam = "bcnPower"
  } else {
      gam <- 0
      fam <- "bcPower"
    }
  if (mean(vec)<0){ print(paste("vector", i,"has negative mean"))}
  else {
  pT <- powerTransform(vec, family = fam)
    if (pT$roundlam != 1){
    lambda <- pT$roundlam
    gamma <- max(pT$gamma,0)
    trans_hn_radiomics[,i] <- bcnPower(vec, lambda = lambda, gamma = gamma)
    print(paste("transformed vector" , i,"with lambda = ", lambda))
    } else { print(paste("vector",i,"doesn't need transform"))}
  }
  hist(trans_hn_radiomics[,i])
}

var_aerts <- apply(hn_radiomics,2,var)

hn_radiomics <- trans_hn_radiomics
hn_radiomics$ShapeVolume <- hn_aerts$ShapeVolume  
hn_radiomics$Age <- hn_aerts$Age
hn_radiomics <- scale(hn_radiomics,T,T)

# N <- N[which(reliability_data$PreFeat %in% colnames(hn_radiomics))]
N <- reliability_data$Control.Total.SD[match(colnames(hn_radiomics), reliability_data$PreFeat)]
# N[41:42] <- mean(N, na.rm = T)
 N[161:162] <- mean(N, na.rm = T)
N <- log(N)
N <- abs(N- max(N))
N <- N/max(N)


## NO GENES
n <- nrow(hn_radiomics)
gene_data <- matrix(rep(0,n*4),ncol = 4)


train_size <- floor(n*3/4)
test_size <- n - train_size
train <- sample(1:n,train_size,replace = F) #seq(1,101,by = 2)  #

            
            
#make groups
groups <- kmeans(hn_clinical$OStime,2)
type <- groups$cluster[train] - 1 
test_type <- groups$cluster[-train] - 1
# trainingOStime <- hn_clinical$OStime[train]
# a <-  max(trainingOStime[groups$cluster == 1])
# b <- min(trainingOStime[groups$cluster == 2])
# border <- mean(a,b)
# type <- groups$cluster - 1 
# test_type <- as.numeric(hn_clinical$OStime[-train] > border)

# what if we do stage?
type2 <- as.integer(hn_clinical$AJCC.stage[train]) -1
test_type2 <- as.integer(hn_clinical$AJCC.stage[-train]) - 1

# what if we do HPV status? 
type3 <- as.integer(hn_clinical$HPV.status[train]) -1
test_type3 <- as.integer(hn_clinical$HPV.status[-train]) - 1


X <- (hn_radiomics[train,])
Xf <- (hn_radiomics[-train,])
Y <- type3
Yf <- test_type3
Z <- as.matrix(gene_data[train,])
Zf <- as.matrix(gene_data[-train,])
N <- t(N)


