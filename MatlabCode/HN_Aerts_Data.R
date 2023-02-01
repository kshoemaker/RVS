### Case Study Data Processing ###
  
# setup 
rm(list = ls())
library(MASS)
library(BAMMtools)
library(tidyverse)
library(gtools)
library(R.matlab)



set.seed(1234)

setwd("../CaseStudyData/")

# 78 subjects
hn_data <- read.csv("CTPatientDataForKate.csv", stringsAsFactors = T)


## load in reliability information 
reliability_data <- read.csv("LungProjectData/TotalSDForKate.csv", header = T, stringsAsFactors = F)




hn_aerts <- filter(hn_data,DataType =="A") ## ,  HPV.status == "Negative")
hn_clinical <- hn_aerts[,1:14] # select the clinical variables
 hn_radiomics <- hn_aerts[,15:186] # select the radiomic variables
  ## chooses only the variables we have reliability data for
hn_radiomics <- hn_radiomics[,which(colnames(hn_radiomics) %in% reliability_data$PreFeat)] 


# uses the suggested box cox power transform on variables
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


N <- reliability_data$Control.Total.SD[match(colnames(hn_radiomics), reliability_data$PreFeat)]

 N[161:162] <- mean(N, na.rm = T) # set N to the mean value for the two clinical features
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


#  HPV status
type <- as.integer(hn_clinical$HPV.status[train]) -1
test_type <- as.integer(hn_clinical$HPV.status[-train]) - 1


X <- (hn_radiomics[train,])
Xf <- (hn_radiomics[-train,])
Y <- type
Yf <- test_type
N <- t(N)


writeMat(con = "ProcessedData.mat", X = X, Xf = Xf,Y = Y, Yf = Yf, N = N)

