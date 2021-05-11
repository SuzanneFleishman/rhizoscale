#This code computes variation partitioning analysis for microbiome analyses
#project repository: https://github.com/SuzanneFleishman/rhizoscale
#publication code created for: (under review, will be updated)

#Code compiled by Suzanne Fleishman; this version saved May 11, 2021
#Based on Borcard, D., Gillet, F., and Legendre, P. (2018). Numerical Ecology with R (Springer). (https://doi.org/10.1007/978-1-4419-7976-6)


#Inputs are Phyloseq objects created in "1PhyloseqObjectCreation"

#To complete this code, run it separately for each phyloseq object (select at line 42)
#most of the code can run straight through, but needs attention at adjust at line #108 and lines 206-307. Adjustments depend on which factors are significant in forward selection


#### Setup and Data Import ####

setwd("xxx")

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
packages <- c("doParallel","phyloseq", "ggplot2","vegan","microbiome","microbiomeSeq", "metagMisc",
              "dplyr","sf","spdep","ape","ade4","adegraphics","adespatial","tmap","dada2","viridis", "corncob","fastDummies")
ipak(packages)


##load augmented PS and spatial data ##
load("xxx/PS_rhizoscale.Rdata")

####.....data prep ####

#select ps object for this run
ps=ps.16s.phy

#remove unclassified genera
gen<-as.data.frame(ps@otu_table)
gen.class<-select(gen,-"unclassified")
ps@otu_table<-otu_table(gen.class, taxa_are_rows = FALSE)

#removeunclassifield class
cla<-as.data.frame(ps@otu_table)
cla.class<-select(cla,-"unclassified")
ps@otu_table<-otu_table(cla.class, taxa_are_rows = FALSE)


#remove unclassifield phyla
phy<-as.data.frame(ps@otu_table)
phy.class<-select(phy,-"unclassified")
ps@otu_table<-otu_table(phy.class, taxa_are_rows = FALSE)





#ensure removal of samples without taxa
ps = prune_samples(sample_sums(ps)>0, ps)


#break down phyloseq objects

counts<-as.data.frame(ps@otu_table)
meta<-as.data.frame(ps@sam_data)


Order2<-dummy_cols(meta$type)[,3] #order 2 is represented by a 1
Dark<-dummy_cols(meta$colorcat)[,2] #dark is represented by a 1
Old<-dummy_cols(meta$agecat)[,2] #old is represented by a 1
data<-cbind(Order2,Dark,Old,meta[,8])



cluster<-as.data.frame(as.matrix(meta[,35:42]))

cord<-cbind(meta$x,meta$y)

colnames(cord)<-c("x","y")

###First test data and see if need to be transformed (reccomended in book)
counts.hel<-decostand(counts,"hellinger") #can be applied to presence absence
decorana(counts.hel)


####Forward selection via RDA on each factor category###

####.....XY RDA####
XY.rda<-rda(counts.hel,cord)
anova(XY.rda)
xy.R2a<-RsquareAdj(XY.rda)$adj.r.squared


#FORWARD SELECTION
XY.fwd<-
  forward.sel(counts.hel,cord,
              adjR2thresh=xy.R2a,
              #adjR2thresh=.99, alpha = 0.1,#use this if xy not significant to double check
              nperm=9999)
XY.sign<-sort(XY.fwd$order)
XY.red<-colnames(cord)[c(XY.sign)]

XY.sign
XY.red

####DECISION########################################################################################################

#if x or y is p<0.10, detrend with these lines
cord.df=as.data.frame(cord)
counts.hel.det<-resid(lm(as.matrix(counts.hel)~.,data=cord.df))


###if not significant, relabel un-detrended coordinates so code can continue
#counts.hel.det<-counts.hel


####.....root traits RDA #####
#Forward selection of env variables
env.rda<-rda(counts.hel~.,data)
env.R2a<-RsquareAdj(env.rda)$adj.r.squared
env.P<-anova(env.rda)


set.seed(1)
env.fwd<-
  forward.sel(counts.hel,data,
              adjR2thresh=env.R2a,
              nperm=9999)

env.sign<-sort(env.fwd$order)
red<-data[c(env.sign)]
env<-colnames(data[env.sign])



####.....Root cluster RDA####

data.cluster<-as.data.frame(as.matrix(cluster))


#Forward selection of cluster variables
cluster.rda<-rda(counts.hel~.,data.cluster)
cluster.R2a<-RsquareAdj(cluster.rda)$adj.r.squared
anova(cluster.rda)


cluster.fwd<-
  forward.sel(counts.hel,data.cluster,
              adjR2thresh=cluster.R2a,
              nperm=9999)

cluster.sign<-sort(cluster.fwd$order)
red.cluster<-data.cluster[c(cluster.sign)]
colnames(data.cluster[cluster.sign])


####.....DBMEM#####

dbmem<-as.data.frame(c(meta[,9:34]))


det.dbmem.rda<-rda(counts.hel.det~.,dbmem)
anova(det.dbmem.rda) 


#if analysis is significant, compute the adjusted R2
det.dbmemR2a<-RsquareAdj(det.dbmem.rda)$adj.r.squared


det.dbmem.fwd<-
  forward.sel(counts.hel.det,
              as.matrix(dbmem),
              adjR2thresh=det.dbmemR2a,
              nperm=9999)

nb.sig.dbmem<-nrow(det.dbmem.fwd)
#id significant dbMEM sorted in increasing order-2
dbmem.sign<-sort(det.dbmem.fwd$order)
#write the significant dbMEM to a new object (reduced set)
dbmem.red<-dbmem[,c(dbmem.sign)]



####variation partitioning ####
#first ensure all matricies are present and record pvals in table
#test each for significance and create pvalue table
set.seed(1)
a.cord<-anova(XY.rda)[1,4]

set.seed(1)
a.rc<-anova(env.rda)[1,4]

set.seed(1)
a.clus<-anova(cluster.rda)[1,4]

set.seed(1)
a.mem<-anova(det.dbmem.rda) [1,4]

initialp<-c(a.cord,a.rc,a.clus,a.mem)

initialp



####VARITION PARTITIONING####

#Done differently depending on how many factor categories are different

                  ###if all factor categories are significant###
                  gen.varpart<-
                    varpart(counts.hel,red.cluster, red, cord, dbmem.red)
                  #venndiagram
                  quartz()
                  showvarparts(4,bg=c("red","blue","yellow","green"))
                  quartz()
                  plot(gen.varpart,
                       digits=2,
                       bg=c("red","blue","yellow","green"))
                  


                  ###IF 2 factor categories are significant###
                  gen.varpart<-
                    varpart(counts.hel,red.cluster, cord, dbmem.red)
                  #venndiagram
                  quartz()
                  showvarparts(3,bg=c("red","blue","yellow"))
                  quartz()
                  plot(gen.varpart,
                       digits=2,
                       bg=c("red","blue","yellow"))
 
                 ###IF 2 factor categories are significant###
                  gen.varpart<-
                    varpart(counts.hel,red.cluster, cord)
                  
                  #venndiagram
                  quartz()
                  showvarparts(2,bg=c("red","blue"))
                  quartz()
                  plot(gen.varpart,
                       digits=2,
                       bg=c("red","blue","yellow"))
                

          


#### Tests of the unique fractions using partial RDAs ####
                  #only do these using factor categories that were used in variation partitioning
a<-rda(counts.hel,red,cbind(cord, red.cluster, dbmem.red))
set.seed(1)
Pchar<-anova(a)$`Pr(>F)`[1]
#FRACTION A - root characteristics

b<- rda(counts.hel,cord,cbind(red, red.cluster, dbmem.red))
set.seed(1)
Pxy<-anova(b)$`Pr(>F)`[1]
#FRaction B - coordinates

c<-rda(counts.hel, red.cluster,cbind(dbmem.red,cord, red))
set.seed(1)
Pcluster<-anova(c)$`Pr(>F)`[1]
#FRaction c - clusters

d<-rda(counts.hel,dbmem.red,cbind( red, cord, red.cluster))
#d<-rda(counts.hel,dbmem.red,cbind(cord,red.cluster))
set.seed(1)
Pmem<-anova(d)$`Pr(>F)`[1]
#FRaction d - MEM




####plot entire RDA####
#way to visualize significant factors along with taxa they may be associated with

cord<-as.data.frame(cord)
colnames(cord)<-c("x","y")

all<-rda(counts,cbind(dbmem.red,cord, red, red.cluster))#HERE remove any factor categories that are not significant
set.seed(1)
anova(all)$`Pr(>F)`[1]
RsquareAdj(all)$adj.r.squared



scl=2 #scaling #2 so that all vectors are coorelated
quartz()
plot(all, type = "n", scaling = scl)
with(cord, points(all, display = "sites", col = "grey",
                       scaling = scl, pch = 21, bg = "grey"))
text(all, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
text(all, display = "bp", scaling = scl, cex = 1, col = "red2")
with(all, legend("topright", legend = "Samples", bty = "n",
                 col = "grey", pch = 21, pt.bg = "grey"))

#use forward selection to see if all factors are significant

fwd<-
  forward.sel(counts.hel,cbind(dbmem.red,cord,red, red.cluster),
              nperm=9999)
fwd.sign<-sort(fwd$order)
fwd.red<-cord[c(fwd.sign)]

fwd
