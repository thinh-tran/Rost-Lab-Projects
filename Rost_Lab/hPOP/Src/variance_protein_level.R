##Calculate variance and select for most varied peptides

##Set Directory 

## Clear Data
rm(list = ls())
library(dplyr)
library(magrittr)

#import data
hpopData <- read.csv("~/Desktop/Rost_Lab/hPOP/Data/hpop2_aligned.mat.csv", header=TRUE)

l <- length(hpopData[1,])-3

# #sort data by protein name then by peptide
# hpopData <- hpopData[order(hpopData$Protein, hpopData$RT_mean, na.last = TRUE, decreasing = TRUE),]
# 
# #pick the top 3 peptides by means 
# peptide <- hpopData[ave(hpopData$RT_mean, hpopData$Protein, FUN = seq_along) <= 3L, ]

#create new dataframe for normalized values containing same # of rows and columns
hpopNorm = as.data.frame(lapply(hpopData, function(x) rep.int(NA, length(x))))
hpopNorm[,1:2] <- hpopData[,1:2]

#median normalize data for each sample:
i <- 3
while (i >= 3 && i <= l-1) {
  med <- median(hpopData[,i], na.rm = TRUE)
  hpopMed <- (hpopData[2:length(hpopData[,i]),i])/med
  hpopNorm[2:length(hpopNorm$Peptide),i] <- hpopMed
  i = i + 1
}

boxplot(log(hpopNorm[,3:l]), outline = FALSE)


##Pick the top 3 peptides of each protein
expmean <- rowMeans(hpopData[,3:l], na.rm = TRUE)
hpopNorm$MeanExpression <- expmean

peptide <- hpopNorm %>% 
  arrange(desc(MeanExpression)) %>% 
  group_by(Protein) %>% 
  slice(1:3)

# boxplot(peptide[, 3:l], outline = FALSE)

#Sum peptide values
peptidevalue_sum <- aggregate(peptide[, 3:l], peptide["Protein"], sum)

# #count and plot the numbers of non-NA peptides in each run
# na_count <- as.numeric(apply(hpopData, 2, function(x) sum(!is.na(x))))
# dotchart(na_count)

#create new column to store variance
hpopVar <- integer(length(peptidevalue_sum[,1]))

#calculate variance for each peptide/row
for (i in c(1:length(peptidevalue_sum$Protein))){
  hpopVar[i] <- var(as.numeric(peptidevalue_sum[i,3:l-2]), na.rm=TRUE)
}

#create new dataframe with peptide name and variance
peptideVariance = data.frame(peptidevalue_sum$Protein, hpopVar)
sorted <- peptideVariance[order(peptideVariance$hpopVar, na.last = TRUE, decreasing = TRUE),]
write.csv(sorted, "peptideVariance_sorted_proteinlevel.csv")


source("http://bioconductor.org/biocLite.R")
biocLite("topGO")


