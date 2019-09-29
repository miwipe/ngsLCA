###############################
###############################
# This R script processes and visualizes outputs of ngsLCA 
# Developed and tested in R version 3.6.1
# Wiki can be found at https://github.com/miwipe/ngsLCA
# Please report bugs and/or suggestions to wyc661217@gmail.com


###############################
###############################
#install and request R packages
packs <- c("vegan","gplots","circlize","reshape","RColorBrewer","analogue","readr","BiocManager","devtools","tidyr","ggplot2")
biocondpacks <- c("ComplexHeatmap")
libs <- .libPaths()
installedPacks <- rownames(installed.packages(lib.loc=libs))
needPacks <- packs[!packs%in%installedPacks]
needPacks2 <- biocondpacks[!biocondpacks%in%installedPacks]

if(length(c(needPacks,needPacks2))>0){
  afile=tempfile()
  cat('\n\n\t->  You need to install packages run the below command\n')
  if(length(needPacks)>0){
    cat("install.packages(c(",file=afile)
    for(i in 1:length(needPacks)){
      cat("\"",needPacks[i],"\"",file=afile,append=T)
      if(i<length(needPacks))
        cat(",",file=afile,append=T)
    }
    cat("))\n",file=afile,append=T)
  }
  if(length(needPacks2)>0){
    cat("BiocManager::install(c(",file=afile,append=T)
    for(i in 1:length(needPacks2)){
      cat("\"",needPacks2[i],"\"",file=afile,append=T)
      if(i<length(needPacks2))
        cat(",",file=afile,append=T)
    }
    cat("))\n",file=afile,append=T)
  }
  tmp<-gsub(" ","",readLines(afile))
  if(length(tmp)>0)
    cat(" \t",tmp[1],"\n")
  if(length(tmp)>1)
    cat("\n\t",tmp[2],"\n")
  unlink(afile)
  quit()
}


cat("\n\n\t-> loading R libraries \n\n")
tmp<-lapply(c(packs,biocondpacks), require, character.only = TRUE, quietly = T)


###############################
###############################
#0
#Parameters read-in
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
  func = "NMDS,group,rarefy,heatmap",
  metadata = NA,
  thr1 = 2,
  thr2 = 5,
  taxa.re = NA,
  sample.re = NA,
  group.name = "10239:Viruses,2157:Archaea,2:Bacteria,4751:Fungi,33090:Viridiplantae,33208:Metazoa",
  thr3 = 0,
  top.abundance = 100
)

##if no argument are given this will print the needed arguments and the default settings for the optional arguments
des<-list(
  path = "working directory containing all lca files",
  func = "\n\t\tfunctions that will be performed; default: NMDS, group, rarefy, heatmap; other option: stratplot (recommend only when metadata are ages)",
  metadata = "full path to your metadata, optional, an example at https://github.com/miwipe/ngsLCA",
  thr1 = "minimum reads number representing a taxa in each sample that considered to be authentic",
  thr2 = "minimum summed reads number across all samples representing a taxa that considered to be authentic",
  taxa.re = "\n\t\t a list of NCBI taxaID representing the taxa that will be removed from final results, taxaID can be found at https://www.ncbi.nlm.nih.gov/Taxonomy",
  sample.re = "a list of file names representing the samples that will be removed from the final results",
  group.name = "\n\t\t higher taxonomic ranks that will be used for grouping the taxa, format: \"NCBI taxaID:Scientific name\"",
  thr3 = "minimum percentage of the reads number of a taxon to the total reads numbers of the group, range from 0 to 1",
  top.abundance = "how many most abundant taxa will be illustrated in figs"
)



#Load arguments
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
cat("thr3 = ", thr3,"\n")
cat("top.abundance = ", top.abundance,"\n")

if (!is.na(func)) {
  func = strsplit(func,",")[[1]]
}

if (!is.na(taxa.re)) {
  taxa.re = as.numeric(strsplit(taxa.re,",")[[1]])
}

if (!is.na(sample.re)) {
  sample.re = strsplit(sample.re,",")[[1]]
}

if (!is.na(group.name)) {
  group.name = strsplit(group.name,",")[[1]]
}

thr1=as.numeric(thr1)
thr2=as.numeric(thr2)
thr3=as.numeric(thr3)


######################################
######################################
#define all functions
#Lca file read in, and replace the file name with metadata if supplied
ReadIn = function(FileList, Blank){
  #creat a black dataframe
  X1 = data.frame(taxa=character(), stringsAsFactors=F)
  #change the read in directory if stating blank is T
  if(Blank){
    path = paste(path, "blanks/", sep="")
  } else{
    path = path
  }
  #read in
  for (i in 1:length(FileList)) {
    
    X2.1 =  read.csv(paste(path, FileList[i], sep=""), header=F, sep="\t", stringsAsFactors=F, fill=T, col.names = paste0("V",seq_len(60)),comment.char = "#")
    
    if(dim(X2.1)[1]>0){
      
      X2.2 = data.frame(taxa=X2.1[,2],
                        count=rep(1, dim(X2.1)[1]),
                        stringsAsFactors = F)
      
      X2.3 = aggregate(X2.2[,2]~X2.2[,1], data=X2.2, FUN=sum)
      colnames(X2.3) = c("taxa",sub(".lca","",FileList[i]))
      X1 = merge(X1, X2.3, by="taxa", all=T)
      
    }
  }
  
  X1[is.na(X1)] = 0
  #Removes samples in the sample removing list
  if (!is.na(sample.re)){
    sample.re = sub(".lca", "", sample.re)
    X1 = X1[,!(colnames(X1) %in% sample.re)]
  }
  #replcing the file names with metadata
  if(!is.na(metadata)){
    
    Mdata = read.csv(metadata, quote="", stringsAsFactors=F, header=F,sep = "\t")
    Mdata$V1 = sub(".lca","",Mdata$V1)
    Mdata = Mdata[order(Mdata$V3),]
    
    N = 0
    for (i in 2:dim(X1)[2]) {
      
      if (length(which(colnames(X1)[i] %in% Mdata$V1))>0) {
        colnames(X1)[i] = Mdata$V2[which(Mdata$V1 %in% colnames(X1)[i])]
        N = c(N,i)
      }
    }
    
    N = N[-1]
    M = which(!2:dim(X1)[2] %in% N) + 1
    X1 = X1[,c(1,N,M)]
    
  }
  #Removes taxa in the taxa.re list
  if (!is.na(taxa.re)) {
    
    TAXAID = as.numeric(sapply(strsplit(X1$taxa, split = ":"),function(x) x[1]))
    if (length(which(TAXAID %in% taxa.re)) >0 ) {
      X1 = X1[-which(TAXAID %in% taxa.re),]
    }
  }
  return(X1)
}


#group taxa into the user defiened groups
GroupTaxa = function(DF, TaxaUnit){
  
  L1 = NULL
  L1[group.name] = list(NULL)
  N=0
  
  for (i in 1:dim(DF)[1]) {
    for (j in 1:length(TaxaUnit)) {
      if (length(grep(TaxaUnit[j], tax.branch[which(tax.branch[,1] %in%  DF[i,1]),])) != 0) {
        L1[[j]] = c(L1[[j]],i)
        N=c(N,i)
      }
    }
  }
  N=N[-1]
  DF_ALL = DF[-N,]
  
  name = data.frame(taxa=names(L1))
  name = name  %>% separate(taxa, sep=":", c("taxa_ID","taxa_name"))
  
  reads_3 = as.data.frame(matrix(nrow =dim(DF)[2]-1, ncol = (dim(name)[1]+2)))
  taxa_3 = as.data.frame(matrix(nrow =dim(DF)[2]-1, ncol = (dim(name)[1]+2)))
  colnames(reads_3) = c("samples","passed_thr3",name$taxa_name)
  colnames(taxa_3) = c("samples","passed_thr3",name$taxa_name)
  reads_3$samples = colnames(DF)[-1]
  taxa_3$samples = colnames(DF)[-1]
  
  for (i in 1:length(L1)) {
    if (length(L1[[i]]) != 0) {
      
      DF1 = DF[L1[[i]],]
      
      if (dim(DF1)[2]>2) {
        for (j in 1:dim(DF1)[1]) {
          DF1[j, which(DF1[j,-1] <= colSums(DF1[,-1]) * thr3)+1] = 0
        }
      }else{
        for (j in 1:dim(DF1)[1]) {
          DF1[j, which(DF1[j,-1] <= sum(DF1[,-1]) * thr3)+1] = 0
        }
      }
      
      DF2 = DF1[rowSums(as.data.frame(DF1[,-1])) != 0,]
      write.table(DF2, file = paste(path, "R_results/intermediate/taxa_groups/", name$taxa_name[i], ".txt", sep=""), quote=F, 
                  row.names=F, col.names=T, sep="\t")
      
      DF2.1 = DF2 %>% separate(taxa, sep=":", c("taxa_ID","taxa_name","taxonomic_rank"))
      write.table(DF2.1, file = paste(path, "R_results/taxa_groups/", name$taxa_name[i], ".txt", sep=""), quote=F, 
                  row.names=F, col.names=T, sep="\t")
      
      DF_ALL = rbind(DF2,DF_ALL)
      
      if (dim(DF1)[2]>2) {
        reads_3[,i+2] = colSums(DF2[,-1])
        DF2[DF2 >0] = 1
        taxa_3[,i+2] = colSums(DF2[,-1])
      }else{
        reads_3[,i+2] = sum(DF2[,-1])
        DF2[DF2 >0] = 1
        taxa_3[,i+2] = sum(DF2[,-1])
      }
    }
  }
  
  write.table(DF_ALL, file = paste(path, "R_results/intermediate/", "taxa_profile_v3.txt", sep=""), quote=F, 
              row.names=F, col.names=T, sep="\t")
  DF_ALL = DF_ALL %>% separate(taxa, sep=":", c("taxa_ID","taxa_name","taxonomic_rank"))
  write.table(DF_ALL, file = paste(path, "R_results/", "clean_taxa_profile.txt", sep=""), quote=F,
              row.names=F, col.names=T, sep="\t")
  
  reads_3[is.na(reads_3)] = 0
  taxa_3[is.na(taxa_3)] = 0
  
  if (dim(DF_ALL)[2]>4) {
    reads_3$passed_thr3 = colSums(DF_ALL[,-c(1:3)])
    DF_ALL[DF_ALL>0] = 1
    taxa_3$passed_thr3 = colSums(DF_ALL[,-c(1:3)])
  }else{
    reads_3$passed_thr3 = sum(DF_ALL[,-c(1:3)])
    DF_ALL[DF_ALL>0] = 1
    taxa_3$passed_thr3 = sum(DF_ALL[,-c(1:3)])
  }
  
  ReadNO[1:dim(reads_3)[1],(dim(ReadNO)[2]+1):(dim(ReadNO)[2]+dim(reads_3)[2]-1)] = reads_3[,-1]
  TaxaNO[1:dim(taxa_3)[1],(dim(TaxaNO)[2]+1):(dim(TaxaNO)[2]+dim(taxa_3)[2]-1)] = taxa_3[,-1]
  
  write.table(ReadNO, file = paste(path, "R_results/readsNO/", "reads_counts.txt", sep=""), quote=F, 
              row.names=F, col.names=T, sep="\t")
  write.table(TaxaNO, file = paste(path, "R_results/readsNO/", "taxa_counts.txt", sep=""), quote=F, 
              row.names=F, col.names=T, sep="\t")
  
}


#Function for clustering taxa into user-defined taxa ranks
Taxa.cluster = function(DF, OutName){
  
  DF1 = data.frame(matrix(vector(), 0, dim(DF)[2]+1, dimnames=list(c(), c("taxa",colnames(DF)))),stringsAsFactors=F)
  colnames(DF1)[2:dim(DF1)[2]] = colnames(DF)
  
  j=1
  for (i in 1:dim(DF)[1]) {
    V1 = grep(Taxa.L, tax.branch[which(tax.branch[,1] %in% rownames(DF)[i]),])
    if (length(V1)>1) {
      V1=V1[1]
    }
    if (length(V1) !=0) {
      DF1[j,] = c(tax.branch[which(tax.branch[,1] %in% rownames(DF)[i]),V1], as.numeric(DF[i,1:dim(DF)[2]]))
      j = j+1
    }
  }
  
  for (i in 2:dim(DF1)[2]) {
    DF1[,i] = as.numeric(DF1[,i])
  }
  
  if (dim(DF1)[1] != 0) {
    DF2 = aggregate(. ~ taxa, data = DF1, sum)
    DF2 = DF2 %>% separate(taxa, sep=":", c("taxa_ID","taxa_name","taxonomic_rank"))
    write.table(DF2, file = paste(path, "R_results/taxa_ranks/", OutName, "_",taxa.level, ".txt", sep=""), quote=F, 
                row.names=F, col.names=T, sep="\t")
  }else{
    cat("\n\n\t-> for", OutName, ":::no taxa clustered into", taxa.level, "\n\n")
  }
}

#function for NMDS
NMDSc = function(DF, OutName){
  
  NMDS.dimention = as.integer((pmin(dim(DF)[1],dim(DF)[2])-1)/2)+1 #caculate the NMDS dimerntion
  
  if (NMDS.dimention >= 2) {
    NMDS = metaMDS(DF, k=NMDS.dimention, trymax=100) #NMDS cluster
    saveRDS(NMDS, file = paste(path, "R_results/NMDS/", OutName, "_NMDS_data.rda", sep=""))
    #output NMDS_stressplot figure
    pdf(paste(path, "R_results/NMDS/", OutName, "_NMDS_stressplot.pdf", sep=""), width=12, height=12) 
    stressplot(NMDS)
    dev.off()
    #output NMDS figure
    X3 = as.data.frame(NMDS$points)[,c(1,2)]
    X4 = as.data.frame(NMDS$species)[,c(1,2)]
    
    p1 = ggplot(X4, aes(MDS1, MDS2))+
      geom_point(size=4,colour="#7F7F7F")+
      geom_text(aes(label = rownames(X4)),colour="black",size=5)+
      labs(title=paste(OutName, "_NMDS", sep=""),
           y="NMDS2", x="NMDS1") +
      theme(axis.title.x = element_text(size = 17),
            title = element_text(size = 19, face="bold"),
            axis.title.y = element_text(size = 17),
            axis.text.x = element_text(size = 13, colour = "black"),
            axis.text.y = element_text(size = 13, colour = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 17),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black",size=0.5),
            panel.border = element_rect(colour = "black",fill = "transparent",size=1.8))+
      geom_density_2d()+
      geom_point(data = X3, 
                 mapping = aes(x = MDS1, y = MDS2),size=1.5,pch=17,colour="grey")
    
    ggsave(plot=p1,height=8,width=10,dpi=500, filename=paste(path, "R_results/NMDS/", OutName, "_NMDS.pdf", sep=""), useDingbats=FALSE)
  } else {
    cat("\n\n\t-> for", OutName, ":::Number of samples or taxa is less than 2, insufficient for NMDS \n\n")
  }
}


#Function for rarefaction
Raref = function(DF, OutName){
  
  if (pmin(dim(DF)[1],dim(DF)[2]) >= 2) {
    
    DF1 = t(DF)
    S = specnumber(DF1) # observed number of taxa
    raremax = min(rowSums(DF1))
    Srare = rarefy(DF1, raremax)
    #
    pdf(paste(path, "R_results/rarefy/", OutName, "_rarefy_scaling.pdf", sep=""), width=12, height=8) 
    plot(S, Srare, xlab = "Observed No. of taxa", ylab = "Rarefied No. of taxa")
    abline(0, 1)
    dev.off()
    #
    pdf(paste(path, "R_results/rarefy/", OutName, "_rarefy_rarecurve.pdf", sep=""), width=12, height=8) 
    rarecurve(DF1, step = 20, sample = raremax, col = "blue", cex = 1.0)
    dev.off()
  } else {
    cat("\n\n\t-> for", OutName, ":::Number of samples or taxa is less than 2, insufficient for rarefaction \n\n")
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
  } else{
    DF = DF[,-dim(DF)[2]]
  }
  
  #heatmap
  F1 = Heatmap(DF, column_dend_height = unit(1.5, "cm"), row_dend_width = unit(3, "cm"), show_row_names = T, show_column_names = T,
               row_names_gp = gpar(cex=0.8), column_names_gp = gpar(cex=0.5), cluster_rows = T, cluster_columns= F, 
               clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",
               col = colorRamp2(c(0, 0.001, 0.01, 0.03, 0.06, 0.2), 
                                c("white", "cornflowerblue", "yellow", "#FD8D3C","#E31A1C","#B10026")),
               rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.5),
               heatmap_legend_param = list(title = "abundance", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8), 
                                           legend_height = unit(10, "cm"), color_bar = "continous"))
  
  return(F1)
}


#analysis on the taxa groups 
groupin = function(Funct){
  
  file.list = dir(paste(path, "R_results/intermediate/taxa_groups/", sep=""), pattern = ".txt") #list files
  
  for (i in 1:length(file.list)) {
    
    DF = read.csv(paste(path, "R_results/intermediate/taxa_groups/", file.list[i], sep=""), sep="\t", quote="",stringsAsFactors=F, check.names=F)
    rownames(DF) = DF$taxa
    DF1 = as.data.frame(DF[,-1])
    colnames(DF1) = colnames(DF)[-1]
    rownames(DF1) = rownames(DF)
    
    if (Funct == "NMDS"){
      NMDSc(DF = DF1, OutName = sub(".txt", "", file.list[i]))
    }
    
    if (Funct == "rarefy"){
      Raref(DF = DF1, OutName = sub(".txt", "", file.list[i]))
    }
    
  }
}



######################################
#1
#data read-in and pre-processes
cat("\n\n\t-> Data read-in and pre-process. Please be patient, it might take some time depending on file number and size\n\n")
dir.create(paste(path, "R_results", sep=""))
dir.create(paste(path, "R_results/intermediate", sep=""))

file.list1 = dir(path, pattern = ".lca$")
file.list2 = dir(paste(path, "blanks/",sep = ""), pattern = ".lca$")

if(length(file.list1)==0){
  cat("\n\n\t-> ERROR: no lca files found in the appointed working directory\n\n")
  q("no")
}

if(length(file.list2)==0){
  cat("\n\n\t-> No \"blanks\" directory under the working directory or no lca files in it, contamination detecting function denied \n\n")
}

if(is.na(metadata)){
  cat("\n\n\t-> Metadata file not supplied and file names will be illustrated \n\n")
} else {
  cat("\n\n\t-> Metadata supplied and will be used \n\n")
}

DF1 = ReadIn(FileList=file.list1, Blank=F)
write.table(DF1, file = paste(path, "R_results/intermediate/", "taxa_profile_v1.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep="\t")


######################################
#2
#Filters and cleans up the sample taxonomic profile 
cat("\n\n\t-> Filter and clean up taxonomic profiles \n\n")
DF2.1 = as.data.frame(DF1[,-1])
colnames(DF2.1) = colnames(DF1)[-1]
DF2.1[DF2.1 < thr1] = 0
DF2.1$taxa  = DF1$taxa
DF2 = DF2.1[,c(dim(DF2.1)[2],2:dim(DF2.1)[2]-1)]
if (dim(DF2)[2]>2) {
  DF2 = DF2[which(rowSums(DF2[,-1]) !=0),]
  DF3 = DF2[which(rowSums(DF2[,-1]) >= thr2),]
}else{
  DF2 = DF2[which(DF2[,2] !=0),]
  DF3 = DF2[which(DF2[,2] >= thr2),]
}


if (length(file.list2)!=0) {
  
  CON = ReadIn(FileList=file.list2, Blank=T)
  CON.1 = as.data.frame(CON[,-1])
  colnames(CON.1) = colnames(CON)[-1]
  CON.1[CON.1 < thr1] = 0
  CON.1$taxa  = CON$taxa
  CON.2 = CON.1[,c(dim(CON.1)[2],2:dim(CON.1)[2]-1)]
  
  if (dim(CON.2)[2]>2) {
    CON.2 = CON.2[which(rowSums(CON.2[,-1]) !=0),]
  }else{
    CON.2 = CON.2[which(CON.2[,2] !=0),]
  }
  
  write.table(CON.2, file = paste(path, "R_results/intermediate/", "contambation_list.txt", sep=""), quote=F, 
              row.names=F, col.names=T, sep="\t")
  
  DF4 = DF3[-which(DF3$taxa %in% CON.2$taxa),]
  
}else{
  DF4 = DF3
}

write.table(DF4, file = paste(path, "R_results/intermediate/", "taxa_profile_v2.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep="\t")


if (length(file.list2)!=0){
  ReadNO = data.frame(sample=c(colnames(DF1)[-1],colnames(CON.2)[-1]),
                      original_Reads=NA,
                      passed_thr1=NA,
                      passed_thr2=NA,
                      contaminate_removed=NA,
                      stringsAsFactors = F)
  
  TaxaNO = data.frame(sample=c(colnames(DF1)[-1],colnames(CON.2)[-1]),
                      taxa_in_original_lca=NA,
                      passed_thr1=NA,
                      passed_thr2=NA,
                      contaminate_removed=NA,
                      stringsAsFactors = F)}else{
                      
                      ReadNO = data.frame(sample=c(colnames(DF1)[-1]),
                                          original_Reads=NA,
                                          passed_thr1=NA,
                                          passed_thr2=NA,
                                          contaminate_removed=NA,
                                          stringsAsFactors = F)
                      
                      TaxaNO = data.frame(sample=c(colnames(DF1)[-1]),
                                          taxa_in_original_lca=NA,
                                          passed_thr1=NA,
                                          passed_thr2=NA,
                                          contaminate_removed=NA,
                                          stringsAsFactors = F)}

if (dim(DF1)[2]>2) {
  
  ReadNO[1:(dim(DF1)[2]-1),2] = colSums(DF1[,-1])
  ReadNO[1:(dim(DF1)[2]-1),3] = colSums(DF2[,-1])
  ReadNO[1:(dim(DF1)[2]-1),4] = colSums(DF3[,-1])
  ReadNO[1:(dim(DF1)[2]-1),5] = colSums(DF4[,-1])
  
  for (i in 2:dim(DF1)[2]) {
    TaxaNO[i-1,2] = length(which(DF1[,i] != 0))
    TaxaNO[i-1,3] = length(which(DF2[,i] != 0))
    TaxaNO[i-1,4] = length(which(DF3[,i] != 0))
    TaxaNO[i-1,5] = length(which(DF4[,i] != 0))
  }
} else{
  
  ReadNO[1:(dim(DF1)[2]-1),2] = sum(DF1[,-1])
  ReadNO[1:(dim(DF1)[2]-1),3] = sum(DF2[,-1])
  ReadNO[1:(dim(DF1)[2]-1),4] = sum(DF3[,-1])
  ReadNO[1:(dim(DF1)[2]-1),5] = sum(DF4[,-1])
  
  for (i in 2:dim(DF1)[2]) {
    TaxaNO[i-1,2] = length(which(DF1[,i] != 0))
    TaxaNO[i-1,3] = length(which(DF2[,i] != 0))
    TaxaNO[i-1,4] = length(which(DF3[,i] != 0))
    TaxaNO[i-1,5] = length(which(DF4[,i] != 0))
  }
}

if (length(file.list2)!=0) {
  
  if (dim(CON)[2]>2) {
    ReadNO[dim(DF1)[2]:dim(ReadNO)[1],2] = colSums(CON[,-1])
    ReadNO[dim(DF1)[2]:dim(ReadNO)[1],3] = colSums(CON.2[,-1])
    
    for (i in 2:dim(CON)[2]) {
      TaxaNO[dim(DF1)[2]-2+i,2] = length(which(CON[,i] != 0))
      TaxaNO[dim(DF1)[2]-2+i,3] = length(which(CON.2[,i] != 0))
    }
  }else{
    ReadNO[dim(DF1)[2]:dim(ReadNO)[1],2] = sum(CON[,-1])
    ReadNO[dim(DF1)[2]:dim(ReadNO)[1],3] = sum(CON.2[,-1])
    
    for (i in 2:dim(CON)[2]) {
      TaxaNO[dim(DF1)[2]-2+i,2] = length(which(CON[,i] != 0))
      TaxaNO[dim(DF1)[2]-2+i,3] = length(which(CON.2[,i] != 0))
    }
  }
}


PATH = read.csv(paste(path, file.list1[1], sep=""), header=F, sep="\t", stringsAsFactors=F, fill=T, col.names = paste0("V",seq_len(60)),comment.char = "#")
PATH = PATH[,-1]
if (length(file.list1)>1) {
  for (i in 2:length(file.list1)) {
    PATH.1 = read.csv(paste(path, file.list1[i], sep=""), header=F, sep="\t", stringsAsFactors=F, fill=T, col.names = paste0("V",seq_len(60)),comment.char = "#")
    PATH.1 = PATH.1[,-1]
    PATH = rbind(PATH,PATH.1)
    PATH = PATH[!duplicated(PATH[,1]),]
  }
}else{
  PATH = PATH[!duplicated(PATH[,1]),]
}
write.table(PATH, file = paste(path, "R_results/intermediate/", "taxa_branch.txt", sep=""), quote=F, 
            row.names=F, col.names=T, sep="\t")


######################################
#3
#Cluster the cleaned taxa profile into groups and ranks
#groups
dir.create(paste(path, "R_results/readsNO", sep=""))

if ("group"%in%func) {
  
  cat("\n\n\t-> grouping taxa\n\n")
  dir.create(paste(path, "R_results/taxa_groups", sep=""))
  dir.create(paste(path, "R_results/intermediate/taxa_groups", sep=""))
  
  tax.branch = read.csv(paste(path, "R_results/intermediate/", "taxa_branch.txt", sep=""), sep="\t", quote="", stringsAsFactors=F)
  DF1 = read.csv(paste(path, "R_results/intermediate/", "taxa_profile_v2.txt", sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
  
  GroupTaxa(DF = DF1, TaxaUnit = group.name)
}else{
  
  X1 = read.csv(paste(path, "results/intermediate/", "clean_taxa_profile_v2.txt", sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
  write.table(X1, file = paste(path, "results/intermediate/", "clean_taxa_profile_v3.txt", sep=""), quote=F, 
              row.names=F, col.names=T, sep="\t")
  
  X2 = X1 %>% separate(taxa, sep=":", c("taxa_ID","taxa_name","taxonomic_rank"))
  write.table(X2, file = paste(path, "results/", "clean_taxa_profile.txt", sep=""), quote=F,
              row.names=F, col.names=T, sep="\t")
  
}

#taxa ranks
dir.create(paste(path, "R_results/taxa_ranks", sep=""))
cat("\n\n\t-> clustering taxa into different taxonomic ranks\n\n")

DF = read.csv(paste(path, "R_results/intermediate/", "taxa_profile_v3.txt", sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
rownames(DF) = DF$taxa
DF1 = as.data.frame(DF[,-1])
colnames(DF1) = colnames(DF)[-1]
rownames(DF1) = rownames(DF)

tax.branch = read.csv(paste(path, "R_results/intermediate/", "taxa_branch.txt", sep=""), sep="\t",quote="", stringsAsFactors=F)

for (taxa.level in c("species","genus","family")) {
  
  Taxa.L = paste(":",taxa.level,sep="")
  X2 = Taxa.cluster(DF = DF1, OutName = "all_taxa")
  
}

if ("group"%in%func) {
  
  file.list = dir(paste(path, "R_results/intermediate/taxa_groups/", sep=""), pattern = ".txt") #list files
  
  for (i in 1:length(file.list)) {
    
    DF = read.csv(paste(path, "R_results/intermediate/taxa_groups/", file.list[i], sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
    rownames(DF) = DF$taxa
    DF1 = as.data.frame(DF[,-1])
    colnames(DF1) = colnames(DF)[-1]
    rownames(DF1) = rownames(DF)
    
    OutName = sub(".txt", "", file.list[i])
    
    for (taxa.level in c("species","genus","family")) {
      
      Taxa.L = paste(":",taxa.level,sep="")
      X2 = Taxa.cluster(DF = DF1, OutName = OutName)
      
    }
  }
}

######################################
######################################
#4
#NMDS
if ("NMDS"%in%func) {
  
  cat("\n\n\t-> NMDS, will take some time depending on input file number and size\n\n")
  dir.create(paste(path, "R_results/NMDS", sep=""))
  
  # NMDS on taxa profiles for all samples
  DF = read.csv(paste(path, "R_results/intermediate/", "taxa_profile_v3.txt", sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
  rownames(DF) = DF$taxa
  DF1 = as.data.frame(DF[,-1])
  colnames(DF1) = colnames(DF)[-1]
  rownames(DF1) = rownames(DF)
  
  NMDSc(DF = DF1, OutName = "all_taxa")
  
  if ("group"%in%func){
    groupin(Funct = "NMDS")
  }
}


######################################
######################################
#5
#rarefaction
if ("rarefy"%in%func){
  
  cat("\n\n\t-> random rarefaction \n\n")
  dir.create(paste(path, "R_results/rarefy", sep=""))
  DF = read.csv(paste(path, "R_results/intermediate/", "taxa_profile_v3.txt", sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
  rownames(DF) = DF$taxa
  DF1 = as.data.frame(DF[,-1])
  colnames(DF1) = colnames(DF)[-1]
  rownames(DF1) = rownames(DF)
  
  Raref(DF=DF1, OutName="all_taxa")
  
  #for the groups produced in step 5
  if ("group"%in%func) {
    groupin(Funct = "rarefy")
  }
  
}


######################################
#6 heatmap 
if ("heatmap"%in%func){
  
  dir.create(paste(path, "R_results/heatmap", sep=""))
  cat("\n\n\t-> draw heatmap \n\n")
  file.list = dir(paste(path, "R_results/taxa_ranks/", sep=""), pattern = ".txt")
  
  for (i in 1:length(file.list)) {
    
    X1 = read.csv(paste(path, "R_results/taxa_ranks/", file.list[i], sep=""),sep="\t", quote="", check.names=F,stringsAsFactors=F)
    row.names(X1) = X1[,2]
    X2 = as.data.frame(X1[,-c(1:3)])
    colnames(X2) = colnames(X1)[-c(1:3)]
    rownames(X2) = rownames(X1)

    Name = sub(".txt", "", file.list[i])
    
    if (dim(X2)[2]>1) {
      pdf(paste(path, "R_results/heatmap/",Name,"_heatmap.pdf", sep=""), width=8, height=13)
      print({
        HeatMap(X2)
      })
      dev.off() 
    } else {
      cat("\n\n\t-> only one sample detected, not drawing heatmap \n\n")
    }
  }
}


######################################
#7
#Genarate files for MEGAN and Krona
dir.create(paste(path, "R_results/megan", sep=""))
dir.create(paste(path, "R_results/krona", sep=""))
X1 = read.csv(paste(path, "R_results/", "clean_taxa_profile.txt", sep=""), sep="\t", quote="",check.names=F,stringsAsFactors=F)
cat("\n\n\t-> generate files for MEGAN and Krona \n\n")
X1 = X1[,-c(1,3)]
colnames(X1)[1] = "#dataset"
write.table(X1, file = paste(path, "R_results/megan/all_taxa.txt", sep=""), quote=F, row.names=F, col.names=T, sep=",")


file.list = dir(paste(path, "R_results/intermediate/taxa_groups/", sep=""), pattern = ".txt")
for (i in 1:length(file.list)) {
  X1 = read.csv(paste(path, "R_results/taxa_groups/", file.list[i], sep=""), sep="\t", quote="",check.names=F,stringsAsFactors=F)
  X1 = X1[,-c(1,3)]
  colnames(X1)[1] = "#dataset"
  write.table(X1, file = paste(path, "R_results/megan/", file.list[i],sep=""), quote=F, row.names=F, col.names=T, sep=",")
  
}

#Generate krona compatible output
#generate a krona style taxomomic path dataframe
tax.branch = read.csv(paste(path, "R_results/intermediate/", "taxa_branch.txt", sep=""), sep="\t", quote="", stringsAsFactors=F)
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
X1 = read.csv(paste(path, "R_results/intermediate/", "taxa_profile_v3.txt", sep=""), quote="",sep="\t",check.names=F,stringsAsFactors=F)

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
    
    write.table(X2, file = paste(path, "R_results/krona/krona_", colnames(X1)[i] ,".txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
    
  }
}
#End
######################################
######################################

