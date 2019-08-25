###############################
#This R script processes and visualizes outputs from ngsLCA 
#Developed and tested in R version 3.5.1 under unix
#Description can be found at https://github.com/miwipe/ngsLCA
#Please report bugs and/or suggestions to ycwang@snm.ku.dk or wyc661217@gmail.com


########## Important ##########
#As input, this script requires
# 1) all lca files (.lca) copied into a working directory
# 2) and negetive control lca files moved into a sub-directory named "blanks" under the working directory.
# A metadata text (sample names, ages, locations, etc.) is optional as input. If provided, the reletive metadata will be illustrated in the results instead of file names
# see https://github.com/miwipe/ngsLCA for an example of metadata format


###############################
#install and request R packages
#this step might require manually handling for the first running

cat("\n ->  loading R libraries \n\n")

if (!require("devtools", quietly=T)) {
  install.packages("devtools")
}



if (!require("vegan", quietly=T) | 
    !require("gplots", quietly=T) | 
    !require("ComplexHeatmap", quietly=T) | 
    !require("circlize", quietly=T) | 
    !require("reshape", quietly=T) | 
    !require("RColorBrewer", quietly=T) | 
    !require("analogue", quietly=T) | 
    !require("readr",quietly=T)) {
  
  install.packages("vegan")
  install.packages("gplots")
  devtools::install_github("jokergoo/ComplexHeatmap")
  install.packages("circlize")
  install.packages("reshape")
  install.packages("RColorBrewer")
  install.packages("analogue")
  install.packages("readr")
  
}


###############################
#0
#arguments read-in
rm(list=(ls())) #clean up global environment
options(warn=-1) #turn off Warning messages

l<-commandArgs(TRUE)

getArgs<-function(x,l){
  unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
}


Args<-function(l,args){
  if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
    cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument\n")
    q("no")
  }
  arguments<-list()
  for(a in names(args))
    arguments[[a]]<-getArgs(a,l)
  
  if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
    cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
    q("no")
  }
  for(a in names(args))
    if(is.null(arguments[[a]]))
      arguments[[a]]<-args[[match(a,names(args))]]
  
  arguments
}

print.args<-function(args,des){
  if(missing(des)){
    des<-as.list(rep("",length(args)))
    names(des)<-names(args)
  }
  cat("->  Needed arguments:\n")
  mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
  cat("->  Optional arguments (defaults):\n")
  mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
  q("no")
}

#NULL is an non-optional argument – please fill, NA is an optional argument, others are the default arguments
args<-list(
  path = NULL,
  func = c("NMDS", "group", "rarefy", "heatmap"),
  metadata = NA,
  thr1 = 2,
  thr2 = 5,
  taxa.re = NA,
  sample.re = NA,
  group.name = c("2:Bacteria", "33630:Alveolata", "33682:Euglenozoa", "4751:Fungi", "33208:Metazoa", "33090:Viridiplantae", "10239:Viruses") ,
  top.abundance = 100
)

##if no argument are given this will print the needed arguments and the default settings for the optional arguments
des<-list(
  path = "working directory contains all lca fiels",
  func = "\n\t\tfunctions that will be applied",
  metadata = "a comma separated file with first column for file names and second column for corresponding metadata, an example at https://github.com/miwipe/ngsLCA",
  thr1 = "minimum reads number representing a taxa that considered to be authentic for each sample",
  thr2 = "minimum summed reads number representing a taxa that considered to be authentic across all samples",
  taxa.re = "\n\t\t high taxonomic ranks and empirical contamination taxa that will be removed",
  sample.re = "file names that will be removed from final taxa profile",
  group.name = "\n\t\t taxonomic unit names that will be used for grouping taxa, format is taxID:scientific name",
  top.abundance = "number of abundant taxa illustrated in figs"
)


#Load arguments and add it to workspace
if(length(l)==0) print.args(args,des)
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
  cat(" Arguments: output prefix\n")
  q("no")
}


cat("path = ", path,"\n")
cat("func = ", func,"\n")
cat("metadata = ", metadata,"\n")
cat("thr1 = ", thr1,"\n")  
cat("thr2 = ", thr2,"\n")
cat("taxa.re = ", taxa.re,"\n")
cat("sample.re = ", sample.re,"\n")
cat("group.name = ", group.name,"\n")
cat("top.abundance = ", top.abundance,"\n")


######################################
######################################
#define all functions

#lca file read-in, merging, counting and filtering
ReadIn = function(FileList, Blank){
  
  DF1 = data.frame(taxa=character(), stringsAsFactors=F)
  
  if(Blank){
    pathF = paste(path, "blanks/", sep="")
  } else{
    pathF = path
  }
  
  for (i in 1:length(FileList)) {
    DF2.1 =  read.csv(paste(pathF, FileList[i], sep=""), header=F, sep="\t", stringsAsFactors=F, fill=T, col.names = paste0("V",seq_len(80)), comment.char = "#")
     if(dim(DF2.1)[1]>0){
    DF2 = matrix(nrow=dim(DF2.1)[1],ncol=2)
    DF2[,2] = rep(1, dim(DF2.1)[1])
    DF2=as.data.frame(DF2)
    DF2[,1] = DF2.1[,2]
    DF3 = aggregate(DF2[,2]~DF2[,1], data=DF2, FUN=sum) #count reads NO. for each taxa
    DF4 = subset(DF3, DF3[ ,2] >= thr1) 
    colnames(DF4)=c("taxa", FileList[i])
    DF1 = merge(DF1, DF4, by="taxa", all=T)
      }
  }
  
  DF1[is.na(DF1)] = 0
  
  for(i in 2:length(colnames(DF1))) #remove the ".lca" in column names
  {
    colnames(DF1)[i] = sub(".lca", "", colnames(DF1)[i])
  }

  #Removes rows with mistakes
  DF2 = DF1[,-1]
  X1 = 0
  for(i in 1:dim(DF2)[1])
  {
    if(any(is.na(as.numeric(DF2[i,]))))
    {
      X1 = c(X1, i)
    }
  }
  
  if (length(X1)>1) {
    X1 = X1[-1]
    DF1 = DF1[-X1,]
  }
  
  return(DF1)
}

#concatenate 2 dataframes
concatenate <- function(DF1, DF2) {
  d1.names <- names(DF1)
  d2.names <- names(DF2)
  d2.add <- setdiff(d1.names, d2.names)
  d1.add <- setdiff(d2.names, d1.names)
  
  if(length(d2.add) > 0) {
    for(i in 1:length(d2.add)) {
      DF2[d2.add[i]] <- NA
    }
  }
  
  if(length(d1.add) > 0) {
    for(i in 1:length(d1.add)) {
      DF1[d1.add[i]] <- NA
    }
  }
  
  return(rbind(DF1, DF2))
}

#taxa profile read-in and pre-process, replacing file names with corresponding metadata in the merged table
datain = function(metadata, full){
  
  if(full){
    X1 = read.csv(paste(path, "results/intermediate/", "taxa_profile(with_controls).txt", sep=""), sep=";", quote="", stringsAsFactors=F, check.names=F)
  } else {
    X1 = read.csv(paste(path, "results/intermediate/", "clean_taxa_profile.txt", sep=""), sep=";", quote="", stringsAsFactors=F, check.names=F)
  }
  
  row.names(X1) = X1[,1]
  X1 = X1[,-1]
  
  if(!is.na(metadata)){
    
    Mdata = read.csv(metadata, quote="", stringsAsFactors=F, header=F)
    Mdata[,1] = sub(".lca", "", Mdata[,1])
    
    if(any(!colnames(X1)%in%Mdata[,1])){
      cat("\n\n\n Error ->   the supplied metadata file does not cover all inputing files, please notice that the controls in \"blanks\" also require supplied metadata")
      q("no")
    }
    
    #Reordering the dataset
    Mdata.1 = Mdata[which(Mdata[,1]%in%colnames(X1)),]
    X1 = X1[,c(Mdata.1[,1])]
    colnames(X1) = Mdata.1[,2]
  }
  
  return(X1)
  
}

#group taxa into the user defiened groups
GroupTaxa = function(DF, TaxaUnit, metadata){
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
      write.table(DF1, file = paste(path, "results/taxa_groups/", names(L1)[i], ".taxa.txt", sep=""), quote=F, 
                  row.names=F, col.names=T, sep="\t")
    }
  }
}


#Function for clustering taxa into user-defined taxa ranks
Taxa.cluster = function(DF, OutName){
  
  DF1 = data.frame(matrix(vector(), 0, dim(DF)[2]+1, dimnames=list(c(), c("taxa",colnames(DF)))),stringsAsFactors=F)
  colnames(DF1)[2:dim(DF1)[2]] = colnames(DF)
  
  j=1
  for (i in 1:dim(DF)[1]) {
    
    V1 = grep(Taxa.L, tax.branch[grep(rownames(DF)[i], tax.branch[,1]),])
    if (length(V1)>1) {
      V1=V1[1]
    }
    
    if (length(V1) !=0) {
      V2 = strsplit(tax.branch[grep(rownames(DF)[i], tax.branch[,1]),V1],":")[[1]][2]
      DF1[j,] = c(V2,DF[i,1:dim(DF)[2]])
      j = j+1
    }
  }
  
  if (dim(DF1)[1] != 0) {
    
    DF2 = aggregate(. ~ taxa, data = DF1, sum)
    write.table(DF2, file = paste(path, "results/taxa_ranks/", OutName, "_",taxa.level, ".txt", sep=""), quote=F, 
                row.names=F, col.names=T, sep="\t")
    
  } else {
    cat("\n ->  ", OutName, ":::no taxa clustered into specific taxa ranks, omitted  \n")
  }
}

#function for NMDS
NMDSc = function(DF, OutName){
  
  NMDS.dimention = as.integer((pmin(dim(DF)[1],dim(DF)[2])-1)/2)+1 #caculate the NMDS dimerntion
  
  if (NMDS.dimention >= 2) {
    NMDS = metaMDS(DF, k=NMDS.dimention, trymax=100) #NMDS cluster
    #output NMDS_stressplot figure
    pdf(paste(path, "results/NMDS/", OutName, "_NMDS_stressplot.pdf", sep=""), width=12, height=12) 
    stressplot(NMDS)
    dev.off()
    #output NMDS figure
    pdf(paste(path, "results/NMDS/", OutName, "_NMDS.pdf", sep=""), width=12, height=12) 
    ordiplot(NMDS$points, type="n", xlab = "NMDS1", ylab = "NMDS2")
    points(as.data.frame(NMDS$species), pch=1, cex=3, col="black")
    orditorp(NMDS,display="species",air=0.01)
    dev.off()
  } else {
    cat("\n\n ->  ", OutName, ":::Number of samples or taxa is less than 2, insufficient for NMDS \n\n")
  }
}


#Function for rarefaction
Raref = function(DF, OutName){
  
  if (pmin(dim(DF)[1],dim(DF)[2]) >= 2) {
    
    DF1 = t(DF[,-1])
    S = specnumber(DF1) # observed number of taxa
    raremax = min(rowSums(DF1))
    Srare = rarefy(DF1, raremax)
    #
    pdf(paste(path, "results/rarefy/", OutName, "_rarefy_scaling.pdf", sep=""), width=12, height=8) 
    plot(S, Srare, xlab = "Observed No. of taxa", ylab = "Rarefied No. of taxa")
    abline(0, 1)
    dev.off()
    #
    pdf(paste(path, "results/rarefy/", OutName, "_rarefy_rarecurve.pdf", sep=""), width=12, height=8) 
    rarecurve(DF1, step = 20, sample = raremax, col = "blue", cex = 1.0)
    dev.off()
  } else {
    cat("\n\n ->  ", OutName, ":::Number of samples or taxa is less than 2, insufficient for rarefaction \n\n")
  }
}

#function for heatmap
HeatMap = function(DF){
  
  #transfer data into percentage and subset
  for (i in 1:dim(DF)[2]) {
    if(sum(DF[,i]) != 0){
      DF[,i] = DF[,i]/sum(DF[,i])
    }
  }
  DF$sum = rowSums(DF)
  DF = DF[order(-DF$sum),]
  
  if (dim(DF)[1] > top.abundance) {
    DF = DF[1:top.abundance,-dim(DF)[2]]
  }
  
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


#Function for stratplot
StratPlot = function(DF,OutName){
  
  #transfer data into percentage and subset
  for (i in 1:dim(DF)[2]) {
    if(sum(DF[,i]) != 0){
      DF[,i] = DF[,i]/sum(DF[,i])
    }
  }
  DF$sum = rowSums(DF)
  DF = DF[order(-DF$sum),]
  
  if (dim(DF)[1] > top.abundance) {
    DF = DF[1:top.abundance,-dim(DF)[2]]
  }
  
  #Stratplot variables
  DF1 = t(DF)
  Ages = as.numeric(rownames(DF1))
  
  if(!any(is.na(Ages))){
    
    pdf(paste(path, "results/heatmap&stratplot/", OutName, "_stratplot.pdf", sep=""), width=15, height=7)
    Stratiplot(x = chooseTaxa(DF1, n.occ = 2, max.abun = 0.05), y = Ages, type = "poly", sort = "wa", xlab="Percentage", ylab = "Ages")
    dev.off()
    
  } else {
    cat("\n\n ->  Metadata file does not contain numeric ages, Stratiplot abandoned \n\n")
  }
  
}


#analysis on the taxa groups 
groupin = function(metadata, Funct){
  
  file.list = dir(paste(path, "results/taxa_groups/", sep=""), pattern = ".txt") #list files
  
  for (i in 1:length(file.list)) {
    
    X1 = read.csv(paste(path, "results/taxa_groups/", file.list[i], sep=""), sep="\t", quote="",stringsAsFactors=F, check.names=F)
    row.names(X1) = X1[,1]
    X1 = X1[,-1]
    
    if (Funct == "NMDS"){
      NMDSc(DF = X1, OutName = sub(".txt", "", file.list[i]))
    }
    
    if (Funct == "rarefy"){
      Raref(DF = X1, OutName = sub(".txt", "", file.list[i]))
    }
    
    if (Funct == "heatmap"){
      X2 = Taxa.cluster(DF = X1, OutName = sub(".txt", "", file.list[i]))
      pdf(paste(path, "results/heatmap&stratplot/", sub(".txt", "", file.list[i]), "_heatmap.pdf", sep=""), width=8, height=13)
      print({
        HeatMap(X2)
      })
      dev.off()
    }
    
    if (Funct == "stratplot"){
      X2 = Taxa.cluster(DF = X1, OutName = sub(".txt", "", file.list[i]))
      StratPlot(DF = X2, OutName = sub(".txt", "", file.list[i]))
    }
    
  }
}


######################################
######################################
#1
#data read-in and pre-processes

#genarate results directory
cat("\n\n ->  Data read-in and pre-process, please be patient, it might take some time depending on file number and size\n\n")
dir.create(paste(path, "results", sep=""))
dir.create(paste(path, "results/intermediate", sep=""))

#samples and coontrols read in seperately
file.list1 = dir(path, pattern = ".lca$") #list sample lca files
file.list2 = dir(paste(path, "blanks/",sep = ""), pattern = ".lca$") #list control lca files

if(length(file.list1)==0){
  cat("\n\n\n Error ->  no lca files found in the appointed working directory\n\n\n")
  q("no")
}

if(length(file.list2)==0){
  cat("\n\n ->  no \"blanks\" under working directory or no lca files found in the \"blanks\" directory, contamanation detecting and removing function denied\n\n")
}

X1.1 = ReadIn(FileList = file.list1, Blank = F)

if(length(file.list2)!=0){
X1.2 = ReadIn(FileList = file.list2, Blank = T)

#merges samples and controls
X1 = merge(X1.1, X1.2, by="taxa", all=T)
X1[is.na(X1)] = 0
}else{
  X1 = X1.1
  X1[is.na(X1)] = 0
}

#writes combined taxonomic profile for each sample and control
write.table(X1, file = paste(path, "results/intermediate/", "taxa_profile(with_controls).txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=";")

#writes a combined taxonomic profile for all samples
write.table(X1.1, file = paste(path, "results/intermediate/", "taxa_profile.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=";")

if(length(file.list2)!=0){
#writes potential contamination list
write.table(X1.2, file = paste(path, "results/intermediate/", "contambation_list.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=";")
}

#produces a full taxonomic branch for each taxa that will be used for taxa clustering afterwards
DF2 = read.csv(paste(path, file.list1[1], sep=""), header=F, sep="\t", stringsAsFactors=F, fill=T, col.names = paste0("V",seq_len(60)))
DF2 = DF2[,-1]

for (i in 2:length(file.list1)) {
  DF2.1 =  read.csv(paste(path, file.list1[i], sep=""), header=F, sep="\t", stringsAsFactors=F, fill=T, col.names = paste0("V",seq_len(60)))
  DF2.1 = DF2.1[,-1]
  
  DF2 = concatenate(DF2, DF2.1)
  DF2 = DF2[!duplicated(DF2[,1]),]
}

write.table(DF2, file = paste(path, "results/intermediate/", "taxa_branch.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=";")

######################################
######################################
#2
#Filters and cleans up the sample taxonomic profile 
cat("\n\n ->  Filter and clean up taxonomic profiles \n\n")

X1 = read.csv(paste(path, "results/intermediate/", "taxa_profile.txt", sep=""), sep=";", quote="", stringsAsFactors=F, check.names=F)

if(length(file.list2)!=0){
  cont.list = read.csv(paste(path, "results/intermediate/", "contambation_list.txt", sep=""), sep=";", quote="", stringsAsFactors=F, check.names=F)
  if (dim(cont.list)[1]>0){
    #Removes taxa from contamanation list
    X2 = X1[-which(X1$taxa %in% cont.list$taxa),]
  }else{
    X2 =X1
  }
}else{
  X2 =X1
}

#Removes lower representing taxa
X3 = X2[,-1]
X3[is.na(X3)] = 0
X4 = 0
for(i in 1:dim(X3)[1])
{
  if(sum(as.numeric(X3[i,])) >= thr2)
  {
    X4 = c(X4, i)
  }
}
X4 = X4[-1]
X5 = X2[X4,]

#Removes taxa in the taxa.re list
X6 = X5[!(X5$taxa %in% taxa.re),]

#Removes samples in the sample removing list
X7 = X6[,!(colnames(X6) %in% sample.re)]

#Writes results to new table
write.table(X7, file = paste(path, "results/intermediate/", "clean_taxa_profile.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep=";")

#read-in metadata if available and replace the column names
if(is.na(metadata)){
  cat("\n\n ->  Metadata file not supplied and file names will be illustrated in figs \n\n")
}

if(!is.na(metadata)){
  cat("\n\n ->  Metadata supplied and will be used in illustrations \n\n")
}

X8 = datain(metadata = metadata, full =F)
X8$taxa = row.names(X8)
X8 = X8[,c(dim(X8)[2],1:(dim(X8)[2]-1))]

write.table(X8, file = paste(path, "results/", "clean_taxa_profile.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep="\t")


######################################
######################################
#3
#Cluster the cleaned taxa profile into groups and taxa ranks
#groups
if ("group"%in%func) {
  
  cat("\n\n ->  grouping taxa \n\n")
  dir.create(paste(path, "results/taxa_groups", sep=""))
  
  tax.branch = read.csv(paste(path, "results/intermediate/", "taxa_branch.txt", sep=""), sep=";", quote="", stringsAsFactors=F)
  X1 = read.csv(paste(path, "results/", "clean_taxa_profile.txt", sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
  
  GroupTaxa(DF = X1, TaxaUnit = group.name, metadata = metadata)
}

#taxa ranks
dir.create(paste(path, "results/taxa_ranks", sep=""))

X1 = datain(metadata = metadata, full=F)
tax.branch = read.csv(paste(path, "results/intermediate/", "taxa_branch.txt", sep=""), sep=";",quote="", stringsAsFactors=F)

for (taxa.level in c("species","genus","family")) {
  
  Taxa.L = paste(":",taxa.level,sep="")
  X2 = Taxa.cluster(DF = X1, OutName = "com_taxa")
  
}

if ("group"%in%func) {
  
  file.list = dir(paste(path, "results/taxa_groups/", sep=""), pattern = ".txt") #list files
  
  for (i in 1:length(file.list)) {
    
    X1 = read.csv(paste(path, "results/taxa_groups/", file.list[i], sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
    rownames(X1)=X1[,1]
    X1=X1[,-1]
    OutName = sub(".taxa.txt", "", file.list[i])
    
    for (taxa.level in c("species","genus","family")) {
      
      Taxa.L = paste(":",taxa.level,sep="")
      X2 = Taxa.cluster(DF = X1, OutName = OutName)
      
    }
  }
}

######################################
######################################
#4
#NMDS
if ("NMDS"%in%func) {
  
  cat("\n\n ->  NMDS on combined taxa profile, will take some time depending on input file number and size \n\n")
  dir.create(paste(path, "results/NMDS", sep=""))
  
  # NMDS on taxa profiles for all samples
  X1 = datain(metadata = metadata, full=F)
  NMDSc(DF = X1, OutName = "complete_NMDS")
  
  if ("group"%in%func){
    groupin(metadata = metadata, Funct = "NMDS")
  }
}


######################################
######################################
#5
#rarefaction
if ("rarefy"%in%func){
  
  cat("\n\n ->  random rarefaction \n\n")
  dir.create(paste(path, "results/rarefy", sep=""))
  X1 = datain(metadata = metadata, full=F)
  
  Raref(DF=X1, OutName="all_taxa")
  
  #for the groups produced in step 5
  if ("group"%in%func) {
    groupin(metadata = metadata, Funct = "rarefy")
  }
  
}


######################################
######################################
#6 heatmap and stratplot
dir.create(paste(path, "results/heatmap&stratplot", sep=""))

#heatmap
if ("heatmap"%in%func){
  
  cat("\n ->  drawing heatmap \n")
  file.list = dir(paste(path, "results/taxa_ranks/", sep=""), pattern = ".txt")
  
  for (i in 1:length(file.list)) {
    
    X1 = read.csv(paste(path, "results/taxa_ranks/", file.list[i], sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
    row.names(X1) = X1[,1]
    X2 = X1[,-1]
    Name = sub(".txt", "", file.list[i])
    
    pdf(paste(path, "results/heatmap&stratplot/",Name,"_heatmap.pdf", sep=""), width=8, height=13)
    print({
      HeatMap(X2)
    })
    dev.off() 
  }
}


#stratplot
if ("stratplot"%in%func){
  
  cat("\n ->  drawing stratplot \n")
  file.list = dir(paste(path, "results/taxa_ranks/", sep=""), pattern = ".txt")
  
  for (i in 1:length(file.list)) {
    
    X1 = read.csv(paste(path, "results/taxa_ranks/", file.list[i], sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
    row.names(X1) = X1[,1]
    X2 = X1[,-1]
    Name = sub(".txt", "", file.list[i])
    
    StratPlot(X2,OutName = Name)
    
  }
}


######################################
######################################
#7
#Genarate files for MEGAN and Krona

#For megan
dir.create(paste(path, "results/megan&krona", sep=""))
X1 = read.csv(paste(path, "results/", "clean_taxa_profile.txt", sep=""), sep="\t", quote="",check.names=F,stringsAsFactors=F)
cat("\n ->  generate files for MEGAN and Krona \n")

for (i in 1:dim(X1)[1]) {
  X1[i,1]  = strsplit(X1[i,1],":")[[1]][2]
}

colnames(X1)[1] = "#dataset"
write.table(X1, file = paste(path, "results/megan&krona/megan_com.txt", sep=""), quote=F, row.names=F, col.names=T, sep=",")



#Generate krona compatible output

#generate a krona style taxomomic path dataframe
tax.branch = read.csv(paste(path, "results/intermediate/", "taxa_branch.txt", sep=""), sep=";", quote="", stringsAsFactors=F)
krona.branch = data.frame(matrix(nrow = dim(tax.branch)[1], ncol = 2), stringsAsFactors=F)

for (i in 1:dim(tax.branch)[1]) {
  
  krona.branch[i,1] = paste(strsplit(tax.branch[i,1],":")[[1]][1:2], collapse = ":")
  
  Taxa.B = strsplit(krona.branch[i,1],":")[[1]][2]
  
  for (j in 2:dim(tax.branch)[2]) {
    if (!is.na(tax.branch[i,j])){
      Taxa.B = paste(strsplit(tax.branch[i,j],":")[[1]][2], Taxa.B, collapse = "/t")
    }
  }
  
  krona.branch[i,2] = Taxa.B
}

#Outputs the most abundant taxa as krona txt file
X1 = read.csv(paste(path, "results/", "clean_taxa_profile.txt", sep=""), quote="",sep="\t",check.names=F,stringsAsFactors=F)

for (i in 2:dim(X1)[2]) {
  X2 = data.frame(X1[,i], X1[,1], stringsAsFactors=F)
  X2 = X2[X2[,1] != 0,]
  X2 = X2[order(-X2[,1]),]
  
  if (dim(X2)[1]>0) {
    
    if (dim(X2)[1] >top.abundance ) {
      X2 = X2[1:top.abundance,]
    }
    
    for (j in 1:dim(X2)[1]) {
      X2[j,2] = paste(strsplit(X2[j,2],":")[[1]][c(1,2)], collapse = ":")
    }
    
    
    for (j in 1:dim(X2)[1]){
      X2[j,2] = krona.branch[which(krona.branch[,1] == X2[j,2], arr.ind=T), 2]
    }
    
    write.table(X2, file = paste(path, "results/megan&krona/krona_", colnames(X1)[i] ,".txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
    
  }
}
#End
######################################
######################################