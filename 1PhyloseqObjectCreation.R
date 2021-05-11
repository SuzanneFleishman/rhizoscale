#Spatial window data analysis initiated 11/10/19
#following prior data exploration, hopefully final tests and visualization

#prior analyses in scripts: DADA2_Pipeline_spatialwindow.R; SpatialAnalysis.R;RDAPhylum.R
#All saved in "spatial window" folder

#Beginning here, located in "spatialDownstream"

#this script preps spatial data and phyloseq objects for downstream analyses

#### Setup and Data Import ####

setwd("/Users/suzannefleishman/Google\ Drive/PSU/Rcodes/spatialDownstream/outputs")

### Clear workspace ###
rm(list=ls())

## load packages we may need
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("doParallel","phyloseq", "ggplot2","vegan","microbiome","microbiomeSeq", "sfsmisc","fastDummies",
              "dplyr","sf","spdep","ape","ade4","adegraphics","adespatial","tmap","dada2","viridis", "corncob","metagMisc")
ipak(packages)



###load data###
#These initial ITS and 16s phyloseq (PS) objects include only metadata as Sample IDs

#load ITS object called "ps.its.ALL"
load("/Users/suzannefleishman/Google\ Drive/PSU/Rcodes/spatialDownstream/DADA2/ITS/ITSspatialPSobj3.2.21.Rdata")

#load 16s object called "ps.16s.bac"
load("/Users/suzannefleishman/Google Drive/PSU/Rcodes/spatialDownstream/DADA2/16s/PSobj16sBacOnly.Rdata")

#load metadata file with coordinates, root trait information, and cluster assignments that map to the Sample IDs
meta<-as.data.frame(read.csv("/Users/suzannefleishman/Google\ Drive/PSU/Rcodes/spatialDownstream/DADA2/metaUpdate3.2.21.csv")) 

####ensure proper metadata for both 16s and its####

sample_names(ps.16s.bac)<-meta$id
sample_names(ps.its.ALL)

sample_sums(ps.its.ALL)


#### Add MEMs to metadata ####
#create dataframe of just x and y coordinates
xadj.tmp<-meta$x
yadj.tmp<-meta$y

cord<-cbind(xadj.tmp,yadj.tmp)

#create MEMs
mem.temp<-dbmem(cord,silent=FALSE,thresh=NULL)
mem<-as.data.frame(mem.temp)


#newmeta
meta.2<-cbind(meta,mem)

meta.3<-meta.2[,-3]

####Make columns with clusters as dummy variables ####

#create dummy variables dataframe
cluster<-dummy_cols(meta$cluster)[,2:9] #all unassigned (NA) are treated as "zero" by selecting just 2:9
colnames(cluster)<-c("a","b","c","d","e","f","g","h")

#
meta.4<-cbind(meta.3,cluster)
meta.5 <- data.frame(meta.4[,-1], row.names=meta.4[,1])

rownames(meta.5)


####add new meta to ps objects####

#convert 16s dataframe to sample data
#no samples need to be removed, because all were sequenced
meta.16s<-as.data.frame(meta.5)
meta.16s<-sample_data(meta.16s)

#check that sample names match between metadata and PS object
sample_names(meta.16s)==sample_names(ps.16s.bac)

#add in metadata and ensure sample names match
ps.16s<-ps.16s.bac
sample_names(meta.16s)<- sample_names(ps.16s)
ps.16s@sam_data<-meta.16s
sample_names(ps.16s)<-sample_names(ps.16s@otu_table)
sample_names(ps.16s@sam_data)<-sample_names(ps.16s@otu_table)

#ITS
#convert ITS dataframe to sample data
#first remove samples that were not sequenced (74, 97-103)

sample_names(ps.its.ALL)

meta.its<-as.data.frame(meta.5[-75,])
meta.its<-as.data.frame(meta.its[1:95,])
meta.its<-sample_data(meta.its)


ps.its<-ps.its.ALL
sample_names(ps.its)<- sample_names(meta.its)
ps.its@sam_data<-sample_data(meta.its)

####prune low samples####
#by looking at sample numbers, it was clear that anything <50 sequences was insufficiently sequenced compared to the others

sample_sums(ps.16s)
ps.16s = prune_samples(sample_sums(ps.16s)>50, ps.16s)

sample_sums(ps.its)
ps.its = prune_samples(sample_sums(ps.its)>50, ps.its)


#### Create taxonomic levels #####
ps.16s.gen<-taxa_level(ps.16s, "Genus")
ps.its.gen<-taxa_level(ps.its, "Genus")

ps.16s.phy<-taxa_level(ps.16s, "Phylum")
ps.its.phy<-taxa_level(ps.its, "Phylum")


ps.16s.cla<-taxa_level(ps.16s, "Class")
ps.its.cla<-taxa_level(ps.its, "Class")

#### Filter for rare ####

#plot then use curve to visualized the distribution of abundance vs prevalence of sequence counts for each ASV
quartz()
phyloseq_prevalence_plot(ps.16s, prev.trh = NULL, taxcolor = NULL,
                         facet = FALSE, point_alpha = 0.7, showplot = T)
quartz()
phyloseq_prevalence_plot(ps.its, prev.trh = NULL, taxcolor = NULL,
                         facet = FALSE, point_alpha = 0.7, showplot = T)



#After testing a few filter options in downstream analyses, a prevalence threshold of 10% was used to remove rare taxa
ps.16s.10 <- phyloseq_filter_prevalence(ps.16s, prev.trh = 0.1, abund.trh = 0.1,
                                        threshold_condition = "AND", abund.type = "total")
ps.its.10 <- phyloseq_filter_prevalence(ps.its, prev.trh = 0.1, abund.trh = 0.1,
                                        threshold_condition = "AND", abund.type = "total")



#### create PA object of 16s and its combined ####

#convert to presence-absence for 10% prevalence and abundance
ps.its.pa<-phyloseq_standardize_otu_abundance(ps.its.10, method = "pa")

ps.16s.pa<-phyloseq_standardize_otu_abundance(ps.16s.10, method = "pa")

#merge taxa tables
ps.merged<-merge_phyloseq(ps.its.pa, ps.16s.pa)

#remove samples that are not present for both ITS and 16S 
remove<-c("sa17","sa41","sa48","sa70", "sa72","sa75","sa76","sa97","sa98","sa99","sa100","sa101","sa102","sa103")

ps.merged.filt<-ps.merged

for(i in 1:length(remove))

{
  ps.merged.filt <- prune_samples(sample_names(ps.merged.filt) != remove[i], ps.merged.filt)
  
}


ps.PA<-ps.merged.filt

####save####

save(ps.16s, ps.16s.10,
     ps.its,ps.its.10,
     ps.16s.gen,ps.16s.phy,ps.16s.cla,
     ps.its.gen,ps.its.phy,ps.its.cla,
     ps.PA,
     file = "PS_rhizoscale.Rdata")
