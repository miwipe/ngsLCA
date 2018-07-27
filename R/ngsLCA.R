######################################
######################################
#This R script is for downstream processing of ngsLAC output datasets
#Developed and tested in R version 3.5.1
#Require all lca files copied into a working directory as imputs
#Then move blank control lca files into a directory named "blanks" under the working directory
#lca file names are recomanded to be as "sample_info.lca", e.g. "age.lca" (8500BP.lca) or "location.lca" (NZ_103.lca)
################


################
#install and request R packages
################
if (!require("devtools")) install.packages("devtools")
if (!require("vegan") | !require("gplots") | !require("ComplexHeatmap") | !require("reshape") | !require("analogue")) {
  install.packages("vegan")
  install.packages("gplots")
  devtools::install_github("jokergoo/ComplexHeatmap")
  install.packages("reshape")
  install.packages("analogue")
}
library(vegan)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(RColorBrewer)
library(analogue)


################
#1
#Define variables
################
rm(list=(ls())) #clean up Global Environment, optional
Path = "/Path/to/your/working/directory/" #The path to working directory
Path = "/Users/yuchengwang/Desktop/ngsLCA_test/"
Thr1 = 2 #for each sample, lowest reads number representing a taxa that considered to be authentic
Thr2 = 5 #for combined sample taxonomic profile, lowest summed reads number representing a taxa that considered to be authentic
Taxa.re = c("1:root","33090:Viridiplantae", "50451:Arabis", "71240:eudicotyledons", "3398:Magnoliophyta")
#A manually edited taxa list that will be removed from sample taxa profile because it is, 
#1) a empirical contamination taxa
#or 2) a taxa that has no damage signal
#or 3) a high taxanomic unit having no sence ecologically
#taxa in the list MUST be the same format as in sample taxa profile ("taxID:name")
Sample.re = c() #fill in lca file names if any that the sample will be removed from final taxa profile
group.name = c("2:Bacteria", "33630:Alveolata", "33682:Euglenozoa", "4751:Fungi", "33208:Metazoa", 
               "33090:Viridiplantae", "10239:Viruses") #taxonomic unit names that will be used for grouping taxa
Taxa.level = "genus" #taxonomic level that will be used for clustering taxa profile for heatmap


################
#2
#Data read in and pre-process
################
#function for data read in and remove taxa represented by 
#reads NO. less than pre-setting threshold in each sample, then merged into one table
ReadIn = function(FileList, Blank){
  DF1 = data.frame(taxa=character(), stringsAsFactors=F)
  if(Blank){
    PathF = paste(Path, "blanks/", sep="")
  } else{
    PathF = Path
  }
  for (i in 1:length(FileList)) { 
    DF2 = read.csv(paste(PathF, FileList[i], sep=""), header=F, sep="\n", stringsAsFactors=F)
    for (j in 1:dim(DF2)[1]) {
      DF2[j,1] = paste(strsplit(DF2[j,1], "\t")[[1]][-1], collapse=",")
    }
    DF2$count = rep(1, dim(DF2)[1])
    DF3 = aggregate(DF2[,2]~DF2[,1], data=DF2, FUN=sum) #count reads NO. for each taxa
    DF4 = subset(DF3, DF3[ ,2] >= Thr1) 
    colnames(DF4)=c("taxa", FileList[i])
    DF1 = merge(DF1, DF4, by="taxa", all=T)
  }
  DF1[is.na(DF1)] = 0
  for(i in 2:length(colnames(DF1))) #remove the ".lca" in column names
  {
    colnames(DF1)[i] = sub(".lca", "", colnames(DF1)[i])
  }
  return(DF1)
}
#samples and coontrols read in seperately, might take some time depending on file number and size
file.list1 = dir(Path, pattern = ".lca") #list sample lca files
X1.1 = ReadIn(FileList = file.list1, Blank = F)
file.list2 = dir(paste(Path, "blanks/",sep = ""), pattern = ".lca") #list control lca files
X1.2 = ReadIn(FileList = file.list2, Blank = T)
#merge samples and controls, produce full taxonomic branch for each taxa
X1 = merge(X1.1, X1.2, by="taxa", all=T)
X1[is.na(X1)] = 0
X2 = data.frame(matrix(ncol = 35, nrow = 0), stringsAsFactors=F)
for (i in 1:dim(X1)[1]) {
  X3 = strsplit(X1[i,1], ",")[[1]]
  X2[i,] = c(X3, rep(NA, 35-length(X3)))
}
dir.create(paste(Path, "results", sep=""))
write.table(X2, file = paste(Path, "results/", "taxa_branch.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=",")
#extract the lowest taxonomic unit for each taxa
ExtacTaxa = function(DF){
  for (i in 1:dim(DF)[1]) {
    DF[i,1] = paste(strsplit(DF[i,1], ":")[[1]][1:2], collapse=":")
  }
  return(DF)
}
#write combined taxonomic profile for samples and controls
write.table(ExtacTaxa(X1), file = paste(Path, "results/", "taxa_profile(with_controls).txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=",")
#write combined taxonomic profile for samples
write.table(ExtacTaxa(X1.1), file = paste(Path, "results/", "taxa_profile.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=",")
#write potential contamination list
write.table(ExtacTaxa(X1.2), file = paste(Path, "results/", "contambation_list.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=",")


################
#3
#filter and clean up sample taxonomic profile 
################
X1 = read.csv(paste(Path, "results/", "taxa_profile.txt", sep=""), quote="", stringsAsFactors=F)
cont.list = read.csv(paste(Path, "results/", "contambation_list.txt", sep=""), quote="", stringsAsFactors=F)
#remove taxa in contamanation list
X2 = X1[-which(X1$taxa %in% cont.list$taxa),]
#remove lower representing taxa
X3 = X2[,-1]
X4 = 0
for(i in 1:dim(X3)[1])
{
  if(sum(X3[i,]) >= Thr2)
  {
    X4 = c(X4, i)
  }
}
X4 = X4[-1]
X5 = X2[X4,]
#remove taxa in the Taxa.re list
X6 = X5[!(X5$taxa %in% Taxa.re),]
#remove samples in the sample removing list
X7 = X6[,!(colnames(X6) %in% Sample.re)]
write.table(X7, file = paste(Path, "results/", "clean_taxa_profile.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=",")


################
#4
#NMDS on complete taxa profile
################
dir.create(paste(Path, "results/NMDS", sep=""))
#function for NMDS, pure numeric dataframe imputing required
NMDSc = function(DF, OutName){
  NMDS.dimention = as.integer((pmin(dim(DF)[1],dim(DF)[2])-1)/2)+1 #caculate the NMDS dimerntion
  if (NMDS.dimention >= 2) {
    NMDS = metaMDS(DF, k=NMDS.dimention, trymax=100) #NMDS cluster
    #output NMDS_stressplot figure
    pdf(paste(Path, "results/NMDS/", OutName, "_NMDS_stressplot.pdf", sep=""), width=12, height=12) 
    stressplot(NMDS)
    dev.off()
    #output NMDS figure
    pdf(paste(Path, "results/NMDS/", OutName, "_NMDS.pdf", sep=""), width=12, height=12) 
    ordiplot(NMDS$points, type="n", xlab = "NMDS1", ylab = "NMDS2")
    points(as.data.frame(NMDS$species), pch=1, cex=3, col="black")
    orditorp(NMDS,display="species",air=0.01)
    dev.off()
  } else {
    print(paste(OutName, "Number of samples or taxa is less than 2, insufficient for NMDS",sep = " ::: "))
  }
}
#based on non-filtered profile containing samples and controls, might take some time depending on size of dataset, optional
X1 = read.csv(paste(Path, "results/", "taxa_profile(with_controls).txt", sep=""), quote="", stringsAsFactors=F)
row.names(X1) = X1[,1]
X1 = X1[,-1]
NMDSc(DF = X1, OutName = "complete_NMDS(withControls)")
#based on clean sample taxa profile 
X1 = read.csv(paste(Path, "results/", "clean_taxa_profile.txt", sep=""), quote="", stringsAsFactors=F)
row.names(X1) = X1[,1]
X1 = X1[,-1]
NMDSc(DF = X1, OutName = "complete_NMDS")


################
#5
#Devide the cleaned taxa profile into groups then NMDS on each group
################
dir.create(paste(Path, "results/groups", sep=""))
X1 = read.csv(paste(Path, "results/", "clean_taxa_profile.txt", sep=""), quote="", stringsAsFactors=F)
tax.branch = read.csv(paste(Path, "results/", "taxa_branch.txt", sep=""), quote="", stringsAsFactors=F)
#function for NMDS clustering
#function for taxa grouping
GroupTaxa = function(DF, TaxaUnit){
  L1 = NULL
  L1[group.name] = list(NULL)
  for (i in 1:dim(DF)[1]) {
    for (j in 1:length(TaxaUnit)) {
      if (length(grep(TaxaUnit[j], tax.branch[grep(DF[i,1], tax.branch[,1]),])) != 0) {
        L1[[j]] = c(L1[[j]],i)
      }
    }
  }
  for (i in 1:length(L1)) {
    if (length(L1[[i]]) != 0) {
      DF1 = DF[L1[[i]],]
      write.table(DF1, file = paste(Path, "results/groups/", names(L1)[i], ".taxa.txt", sep=""), quote=F, 
                  row.names=F, col.names=T, sep=",")
      row.names(DF1) = DF1[,1]
      DF1 = DF1[,-1]
      NMDSc(DF = DF1, OutName = names(L1)[i])
    }
  }
}
##
GroupTaxa(DF = X1, TaxaUnit = group.name)


################
#6 
#Rarefaction
################
dir.create(paste(Path, "results/rarefy", sep=""))
#function for rarefaction
Raref = function(DF, OutName){
  if (pmin(dim(DF)[1],dim(DF)[2]) >= 2) {
    DF1 = t(DF[,-1])
    S = specnumber(DF1) # observed number of taxa
    raremax = min(rowSums(DF1))
    Srare = rarefy(DF1, raremax)
    #
    pdf(paste(Path, "results/rarefy/", OutName, "_rarefy_scaling.pdf", sep=""), width=12, height=8) 
    plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
    abline(0, 1)
    dev.off()
    #
    pdf(paste(Path, "results/rarefy/", OutName, "_rarefy_rarecurve.pdf", sep=""), width=12, height=8) 
    rarecurve(DF1, step = 20, sample = raremax, col = "blue", cex = 1.0)
    dev.off()
  } else {
    print(paste(OutName, "Number of samples or taxa is less than 2, insufficient for rarefaction",sep = " ::: "))
  }
}
#for the clean taxa profile
X1 = read.csv(paste(Path, "results/", "clean_taxa_profile.txt", sep=""), quote="",stringsAsFactors=F)
rownames(X1) = X1[,1]
Raref(DF=X1, OutName="all_taxa")
#for the groups in step 5
file.list = dir(paste(Path, "results/groups/", sep=""), pattern = ".txt") #list sample lca files
for (i in 1:length(file.list)) {
  X1 = read.csv(paste(Path, "results/groups/", file.list[i], sep=""), quote="",stringsAsFactors=F)
  rownames(X1) = X1[,1]
  Raref(DF = X1, OutName = sub(".txt", "", file.list[i]))
}


################
#7 
#cluster the taxa into a specific taxonomic level then produce heatmaps and stratplot
################
dir.create(paste(Path, "results/heatmap&stratplot", sep=""))
tax.branch = read.csv(paste(Path, "results/", "taxa_branch.txt", sep=""), quote="", stringsAsFactors=F)
Taxa.L = paste(":",Taxa.level,sep="")
#function for clustering taxa
Taxa.cluster = function(DF, OutName){
  DF1 = data.frame(matrix(vector(), 0, dim(DF)[2], dimnames=list(c(), colnames(DF))),stringsAsFactors=F)
  j=1
  for (i in 1:dim(DF)[1]) {
    V1 = grep(Taxa.L, tax.branch[grep(DF[i,1], tax.branch[,1]),])
    if (length(V1) !=0) {
      V2 = strsplit(tax.branch[grep(DF[i,1], tax.branch[,1]),V1],":")[[1]][2]
      DF1[j,] = c(V2,DF[i,2:dim(DF)[2]])
      j = j+1
    }
  }
  DF2 = aggregate(. ~ taxa, data = DF1, sum)
  write.table(DF2, file = paste(Path, "results/heatmap&stratplot/", OutName, Taxa.level, ".txt", sep=""), quote=F, 
              row.names=F, col.names=T, sep=",")
  return(DF2)
}
#function for heatmap
HeatMap = function(DF){
  #transfer data into percentage, subset the top 100 abundant taxa
  rownames(DF) = DF[,1]
  DF = DF[,-1]
  for (i in 1:dim(DF)[2]) {
    DF[,i] = DF[,i]/sum(DF[,i])
  }
  DF$sum = rowSums(DF)
  DF = DF[order(-DF$sum),]
  DF = DF[1:100,-dim(DF)[2]]
  #heatmap
  DF1 = melt(DF)
  F1 = Heatmap(DF, column_dend_height = unit(1.5, "cm"), row_dend_width = unit(3, "cm"), show_row_names = T, show_column_names = T,
          row_names_gp = gpar(cex=0.5), column_names_gp = gpar(cex=0.5), cluster_rows = T, cluster_columns= F, 
          clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",
          col = colorRamp2(c(max(DF1$value), quantile(DF1$value, c(3)/4), quantile(DF1$value, c(2)/4), 
                             quantile(DF1$value, c(1)/4), min(DF1$value)), brewer.pal(n=5, name="Spectral")), 
          heatmap_legend_param = list(title = "abundance", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8), 
                                      legend_height = unit(10, "cm"), color_bar = "continous"))
  return(F1)
}
#Function for stratplot, only valid when the original file names are in ages
StratPlot = function(DF,OutName){
  #transfer data into percentage, subset the top 100 abundant taxa
  rownames(DF) = DF[,1]
  DF = DF[,-1]
  for (i in 1:dim(DF)[2]) {
    DF[,i] = DF[,i]/sum(DF[,i])
  }
  DF$sum = rowSums(DF)
  DF = DF[order(-DF$sum),]
  DF = DF[1:100,-dim(DF)[2]]
  #stratplot
  DF1 = t(DF)
  Ages = as.numeric(rownames(DF2))
  pdf(paste(Path, "results/heatmap&stratplot/", OutName, "_stratplot.pdf", sep=""), width=15, height=7)
  Stratiplot(x = chooseTaxa(DF1, n.occ = 2, max.abun = 0.05), y = Ages, type = "poly", sort = "wa", xlab="Percentage", ylab = "Ages")
  dev.off()
}
#process the complete taxa profile
X1 = read.csv(paste(Path, "results/", "clean_taxa_profile.txt", sep=""), quote="",stringsAsFactors=F)
X2 = Taxa.cluster(DF = X1, OutName = "com_taxa")
StratPlot(X2,OutName = "com_taxa")
X3 = HeatMap(X2)
pdf(paste(Path, "results/heatmap&stratplot/com_taxa_heatmap.pdf", sep=""), width=8, height=13)
X3
dev.off()
#for the groups in step 5
file.list = dir(paste(Path, "results/groups/", sep=""), pattern = ".txt") #list sample lca files
for (i in 1:length(file.list)) {
  X1 = read.csv(paste(Path, "results/groups/", file.list[i], sep=""), quote="",stringsAsFactors=F)
  X2 = Taxa.cluster(DF = X1, OutName = sub(".txt", "", file.list[i]))
  StratPlot(X2,OutName = sub(".txt", "", file.list[i]))
  X3 = HeatMap(X2)
  pdf(paste(Path, "results/heatmap&stratplot/", sub(".txt", "", file.list[i]), "_heatmap.pdf", sep=""), width=8, height=13)
  X3
  dev.off()
}


################
#8
#produce files for megan and krona
################
#for megan
dir.create(paste(Path, "results/megan&krona", sep=""))
X1 = read.csv(paste(Path, "results/", "clean_taxa_profile.txt", sep=""), quote="",stringsAsFactors=F)
for (i in 1:dim(X1)[1]) {
  X1[i,1]  = strsplit(X1[i,1],":")[[1]][2]
}
colnames(X1)[1] = "#dataset"
write.table(X1, file = paste(Path, "results/megan&krona/megan_com.txt", sep=""), quote=F, row.names=F, col.names=T, sep=",")
#for krona
#produce a krona style taxomomic path dataframe
tax.branch = read.csv(paste(Path, "results/", "taxa_branch.txt", sep=""), quote="", stringsAsFactors=F)
krona.branch = data.frame(matrix(nrow = dim(tax.branch)[1], ncol = 2), stringsAsFactors=F)
for (i in 1:dim(tax.branch)[1]) {
  krona.branch[i,1] = paste(strsplit(tax.branch[i,1],":")[[1]][1:2], collapse = ":")
  Taxa.B = strsplit(krona.branch[j,1],":")[[1]][2]
  for (j in 2:35) {
    if (!is.na(tax.branch[i,j])){
      Taxa.B = paste(strsplit(tax.branch[i,j],":")[[1]][2], Taxa.B, collapse = "/t")
    }
  }
  krona.branch[i,2] = Taxa.B
}
#output the top 100 abandanced taxa as krona txt file
X1 = read.csv(paste(Path, "results/", "clean_taxa_profile.txt", sep=""), quote="",stringsAsFactors=F)
for (i in 2:dim(X1)[2]) {
  X2 = data.frame(X1[,i], X1[,1], stringsAsFactors=F)
  X2 = X2[X2[,1] != 0,]
  X2 = X2[order(-X2[,1]),]
  X2 = X2[1:100,]
  for (j in 1:dim(X2)[1]){
    X2[j,2] = krona.branch[which(krona.branch[,1] == X2[j,2], arr.ind=T), 2]
  }
  write.table(X2, file = paste(Path, "results/megan&krona/krona_", colnames(X1)[i] ,".txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
}


################
#End
######################################
######################################
