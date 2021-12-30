# Hemimethylation in Breast Cancer Cell Lines
# Noah Ledbetter
# Date of Last Update: 12/1/2021

# Notes about running
# Data should be in your working directory
# This will probably take about 2 hours to run start to finish


# Libraries

library(exactRankTests)
library(coin)
library(glue)
library(ggplot2)
library(RcppRoll)
library(Rcpp)



# Defining Functions that will be used in the program
delete.na <- function(DF, n=4) {
  # deletes rows that have more than n missing values for either forward or reverse data
  DF[rowSums(is.na(DF[,3:9])) <= n & rowSums(is.na(DF[,11:17])) <= n,]
}


pval_meandif_calculation <- function(df) {
  #Caclulates the pvalue with wilcoxsign_test and the mean difference
  vec = vector(length = (dim(df)[1]))
  for (row in 1:dim(df)[1]){
    tryCatch(
      {
        out = pvalue(wilcoxsign_test(as.numeric(df[row,3:9]) ~ as.numeric(df[row,11:17]), alternative = "two.sided", zero.method = c("Wilcoxon"), distribution = exact()))
      }, error=function(e) {return(conditionMessage(e))}, warning = function(w) {return(conditionMessage(w))})
    vec[row] <- out
  }
  df$pval = vec
  df$meandiff =rowMeans(df[,3:9],na.rm = TRUE) - rowMeans(df[row,11:17],na.rm = TRUE)
  return (df)
}


cluster_hist <- function(df){
  #creates a histogram of the length of a cpg cluster
  cpgSiteLength = vector(length = length(split(df,df$groupNumber)))
  index = 1
  for (group in split(df,df$groupNumber)){
    if (group$groupNumber[1] != 0){
      cpgSiteLength[index]<-group[1]$POSITION[length(group[1]$POSITION)] - group[1]$POSITION[1]
      index = index + 1
    }
  }
  dfName = df$CHR.x[1]
  main = glue("{dfName} Histogram of Cluster Length")
  cpgClusters = unlist(cpgSiteLength[which(cpgSiteLength > 0)])
  jpeg(file = glue("{main}.jpg"))
  hist(cpgClusters, breaks = 20, main = main, xlab = "Length")
  dev.off()
}


cluster_hist_zoomed <- function(df){
  #creates a histogram of the length of a cpg cluster, removing clusters > 100 in length
  cpgSiteLength = vector(length = length(split(df,df$groupNumber)))
  index = 1
  for (group in split(df,df$groupNumber)){
    if (group$groupNumber[1] != 0){
      cpgSiteLength[index]<-group[1]$POSITION[length(group[1]$POSITION)] - group[1]$POSITION[1]
      index = index + 1
    }
  }
  dfName = df$CHR.x[1]
  main = glue("{dfName} Histogram of Cluster Length < 100")
  cpgClusters = unlist(cpgSiteLength[which(cpgSiteLength > 0 & cpgSiteLength <= 100 )])
  jpeg(file = glue("{dfName}zoomed.jpg"))
  hist(cpgClusters, breaks = 20, main = main, xlab = "Length")
  dev.off()
}


grouping <- function(df){
  # Creates groups with ungrouped items assigned to group 0
  df$adj <- rownames(df)
  groupNumber = 1
  df$groupNumber = df$adj
  for (group in (split(rownames(df), cumsum(c(1,diff(as.integer(rownames(df))) != 1))))){
    if (length(group) > 1) {
      df[df$adj %in% group,'groupNumber'] = groupNumber
      groupNumber = groupNumber + 1
    } else{
      df[df$adj %in% group,'groupNumber'] = 0
    }
  }
  return(df)
}


size_of_df <- function(df,nmd,np,rmd,rp){
  # Returns the size of a data frame when filtered under certain conditions
  size = dim(subset(df, (meandiff > nmd & pval <np) | (meandiff > rmd & pval <rp & groupNumber != 0)))[1]
  return(size)
}

load_in_data <- function(chrm){
  # Loads in chr data
  fileFor <- read.table(paste0("cancer.hg19.for.chr", chrm,".txt"), header=TRUE, sep = "")
  fileRev <- read.table(paste0("cancer.hg19.rev.chr", chrm,".txt"), header=TRUE, sep = "")
  merged <- merge(fileFor, fileRev, by = "POSITION")
  return(merged)
}

filter_data <- function(df){
  # Filters the data
  newDf = subset(df, (pval< .05 & meandiff >.4) | (pval< .13 & meandiff >.4 & groupNumber != 0))
  return(newDf)
}



composite <- function(df){
  # Runs through all of the previous funtions
  print("Removing missing values")
  df <- delete.na(df)
  print("calculating mean dif")
  df <- pval_meandif_calculation(df)
  print("Creating Groups")
  df <- grouping(df)
  df2 <- filter_data(df)
  if (dim(df2)[1] > 1 & dim(df2)[2] > 1){
    print("Making Hist")
    cluster_hist(df)
    cluster_hist_zoomed(df)
  }
  return(df)
}



# Main
sizeListMaster = data.frame(index = 1:651)
for (i in 1:22){
  print(i)
  merged <- load_in_data(i)
  merged <- composite(merged)
  sizeList = list()
  pvalList = list()
  meandifList = list()
  for (pval in seq(.05,.2,.005)){
    for (meandif in seq(.4-.2,.4,.01)){
      pvalList=append(pvalList,pval)
      meandifList=append(meandifList,meandif)
      sizeList=append(sizeList,size_of_df(merged,.4,.05,meandif,pval))
    }
  }
  sizeListMaster[, ncol(sizeListMaster) + 1] <- unlist(sizeList)
}

## Seperate CHR X
merged <- load_in_data("x")
merged <- composite(merged)
sizeList = list()
pvalList = list()
meandifList = list()
for (pval in seq(.05,.2,.005)){
  for (meandif in seq(.4-.2,.4,.01)){
    pvalList=append(pvalList,pval)
    meandifList=append(meandifList,meandif)
    sizeList=append(sizeList,size_of_df(merged,.4,.05,meandif,pval))
  }
}
sizeListMaster[, ncol(sizeListMaster) + 1] <- unlist(sizeList)
colnames(sizeListMaster) <-c("index",paste0("chr",1:22), "chrx")
sizeListMaster$pval <- as.numeric(pvalList)
sizeListMaster$MeanDifference <- as.numeric(meandifList)

## Making Filter Chart
jpeg(file = "Cluster Criteria.jpeg")
ggplot(data = sizeListMaster, aes(x = pval, color = MeanDifference)) + 
  geom_point(shape = 1, aes(y = chr1))+
  geom_point(shape = 2, aes(y = chr2))+
  geom_point(shape = 3, aes(y = chr3)) +
  geom_point(shape = 4, aes(y = chr4))+
  geom_point(shape = 5, aes(y = chr5)) +
  geom_point(shape = 6, aes(y = chr6))+
  geom_point(shape = 7, aes(y = chr7)) +
  geom_point(shape = 8, aes(y = chr8))+
  geom_point(shape = 9, aes(y = chr9)) +
  geom_point(shape = 10, aes(y = chr10))+
  geom_point(shape = 11, aes(y = chr11)) +
  geom_point(shape = 12, aes(y = chr12))+
  geom_point(shape = 13, aes(y = chr13)) +
  geom_point(shape = 14, aes(y = chr14))+
  geom_point(shape = 15, aes(y = chr15)) +
  geom_point(shape = 16, aes(y = chr16))+
  geom_point(shape = 17, aes(y = chr17)) +
  geom_point(shape = 18, aes(y = chr18))+
  geom_point(shape = 19, aes(y = chr19)) +
  geom_point(shape = 20, aes(y = chr20))+
  geom_point(shape = 21, aes(y = chr21)) +
  geom_point(shape = 22, aes(y = chr22))+
  geom_point(shape = 23, aes(y = chrx))+
  geom_vline(xintercept = .13, color = "red")+
  ggtitle("Cluster Criteria") +
  ylab("Size")
dev.off()
