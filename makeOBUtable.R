###KS###
#make a empty dataframe called KS_OBUtable#
ks_otutable_IDisASV <- read.table(file="ks_otutable_IDisASV.csv", sep=';', row.names=1, header=TRUE)

clustertable <- read.table(file="clustertable.csv", sep=';', header=TRUE)

KS_OBUtable = data.frame()

clusters <- unique(clustertable$cluster_number)


for (i in 1:length(clusters))
{
  ASVs <- (clustertable[clustertable$cluster_number == clusters[i],"ASV_ID"])
  
  OTUtable.cluster <- ks_otutable_IDisASV[row.names(ks_otutable_IDisASV) %in% ASVs,]
  
  newrow <- colSums(OTUtable.cluster)
  
  newrow <- data.frame(t(newrow))
  
  row.names(newrow) <- clusters[i]
  
  names(newrow) <- names(ks_otutable_IDisASV)
  
  KS_OBUtable <- rbind(KS_OBUtable,newrow)
}

write_csv(KS_OBUtable, "KS_OBUtable.csv")


###AD###
AD_otutable_IDisASV <- read.table(file="AD_otutable_IDisASV.csv", sep=';', row.names=1, header=TRUE)

clustertable <- read.table(file="clustertable.csv", sep=';', header=TRUE)

AD_OBUtable = data.frame()

clusters <- unique(clustertable$cluster_number)


for (i in 1:length(clusters))
{
  ASVs <- (clustertable[clustertable$cluster_number == clusters[i],"ASV_ID"])
  
  OTUtable.cluster <- AD_otutable_IDisASV[row.names(AD_otutable_IDisASV) %in% ASVs,]
  
  newrow <- colSums(OTUtable.cluster)
  
  newrow <- data.frame(t(newrow))
  
  row.names(newrow) <- clusters[i]
  
  names(newrow) <- names(AD_otutable_IDisASV)
  
  AD_OBUtable <- rbind(AD_OBUtable,newrow)
}

write_csv(AD_OBUtable, "AD_OBUtable.csv")
