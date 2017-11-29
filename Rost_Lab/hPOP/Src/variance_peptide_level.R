##Calculate variance and select for most varied peptides

##Set Directory 

## Clear Data
rm(list = ls())

#import data
hpopData <- read.csv("~/Desktop/Rost_Lab/hPOP/Data/hpop2_aligned.mat.csv", header=TRUE)

l <- length(hpopData[1,])-3

#plot peptide abundance per sample
boxplot(log(hpopData[,3:l]), outline = FALSE)

#create new dataframe for normalized values containing same # of rows and columns

hpopNorm = as.data.frame(lapply(hpopData, function(x) rep.int(NA, length(x))))
hpopNorm[,1:2] <- hpopData[,1:2]

#median normalize data for each sample:
i <- 3
while (i >= 3 && i <= l) {
  med <- median(hpopData[,i], na.rm = TRUE)
  hpopMed <- (hpopData[2:length(hpopData[,i]),i])/med
  hpopNorm[2:length(hpopNorm$Peptide),i] <- hpopMed
  i = i + 1
}

boxplot(log(hpopNorm[,3:l]), outline = FALSE)

# #count and plot the numbers of NA peptides in each run
# na_count <- as.numeric(apply(hpopData, 2, function(x) sum(!is.na(x))))
# na_count
# dotchart(na_count)

#
#create new column to store variance
hpopVar <- integer(length(hpopData[,1]))

#calculate variance for each peptide/row
for (i in c(1:length(hpopData$Peptide))){
  hpopVar[i] <- var(as.numeric(hpopNorm[i,3:l]), na.rm=TRUE)
}

#calculate covariance for each peptide/row
hpopCovar <- integer(length(hpopData[,1]))
for (i in c(1:length(hpopData$Peptide))){
  hpopCovar[i] <- cov(as.numeric(hpopNorm[i,3:l]))
}
  

#create new dataframe with peptide name and variance
peptideVariance = data.frame(hpopData$Peptide, hpopData$Protein, hpopVar)
sorted <- peptideVariance[order(peptideVariance$hpopVar, na.last = TRUE, decreasing = TRUE),]
write.csv(sorted, "peptideVariance_sorted_peptidelevel.csv")
