#This code computes shannon diversity and tests it using a mixed model
#project repository: https://github.com/SuzanneFleishman/rhizoscale
#publication code created for: (under review, will be updated)

#Code compiled by Suzanne Fleishman; this version saved May 11, 2021
#Based on Borcard, D., Gillet, F., and Legendre, P. (2018). Numerical Ecology with R (Springer). (https://doi.org/10.1007/978-1-4419-7976-6)


#Inputs are Phyloseq objects created in "1PhyloseqObjectCreation"

#To complete this code, run it separately for all phyloseq objectss; select at L41
#to compare bacterial and fungal diversity, see L137
#code needs attention at each step

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
packages <- c("nlme","lme4", "effects","corncob","magrittr","car", "interactions","RVAideMemoire","tibble","doParallel","phyloseq", "ggplot2","vegan","microbiome","microbiomeSeq",
              "dplyr","Biostrings","ShortRead","ape","ade4","adegraphics","adespatial","tmap","dada2","jtools","rcompanion","ggnewscale","RColorBrewer","eemeans")
ipak(packages)



##load augmented PS and spatial data ##
load("xxx/PS_rhizoscale.Rdata")


####Select PS object for this iteration ####
ps<-ps.PA


# PS objects were not ultimately rarefied for this study, but it was done to double check results were similar (they were)
#set.seed(500)
#ps.rare<-rarefy_even_depth(ps.pruned)
#sample_sums(ps.rare)

#reassign PS object, since it was not rarefied
ps.rare<-ps.16s

####breakdown ps rare into components####

counts<-as.data.frame(ps.rare@otu_table)

meta<-as.data.frame(ps.rare@sam_data)

#### Calculate alphadiversity ####
rich <- estimate_richness(ps.rare,measures="Shannon")
richtable<-cbind(meta,rich)



#### Examine for normality ####

quartz()
qqnorm(richtable$Shannon, pch = 1, frame = FALSE)
qqline(richtable$Shannon, col = "steelblue", lwd = 2)


####exploratory plots ####

#manipulate factors to first visualize data
ggplot(richtable,aes(agecat, Shannon, col=agecat))+
  geom_boxplot()+
  geom_point()


#### test ####
#for this study, only MEMs significant in "2VariationPartitionAnalysis" were used to prevent overfitting

#reverse selection was used with anything under 0.25 remaining in the final model

reg.sh=lme(Shannon~agecat+type+MEM9, random=~1|cluster, data = richtable) #cluster used as random factor in order to prevent colinearity
sh.pval<-as.data.frame(anova.lme(reg.sh))#type 3 sums of squares, because of unequal sample sizes
sh.pval


####plotsignificant factors####

#boxplot for categorical
quartz()
ggplot(richtable,aes(type, Shannon, fill=type))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75),dotsize=.3)+
  scale_fill_manual(values=c("darkgrey", "lightgrey","snow4","white"))


#line for continuous

qu.pval<-as.data.frame(anova.lme(reg.sh))#type 3 sums of squares
qu.e<-allEffects(reg.sh,xlevels=8)

qu.1<-as.data.frame(qu.e[[2]]) #0.95 default confidence
qu.1$var<-"x(cm)"

colnames(qu.1)<-c("val","fit","se","lower","upper","var")

pont<-as.data.frame(cbind(richtable$MEM8,richtable$Shannon))
colnames(pont)<-c("x","Shannon")

g<-ggplot(qu.1,aes(val, fit,col=var)) + 
  geom_line(data=qu.1,aes(x=val,y=fit),linetype = 1) +
  geom_ribbon(data=qu.1,aes(val, fit, ymin = lower, ymax = upper),#,fill=var), 
              linetype="blank",
              alpha = .1,show.legend=FALSE)+
  scale_color_manual(values=c("#D55E00"))+
  scale_fill_manual(values=c("#D55E00"))+
  geom_point(data=pont,aes(x=x,y=Shannon),color="#D55E00")

quartz()
g


####fungi vs bacteria####
##run code 2x separately to get a table for its and 16s and assign "richtable" to new dataframe name

richtable.its.gen<-richtable
richtable.16s.gen<-richtable

##merge tables
x<-as.data.frame(richtable.its.gen)#fungi
y<-as.data.frame(richtable.16s.gen) #bacteria

richtable.combined<-merge(x, y,
                          by.x = 'traceid', by.y = 'traceid', all = FALSE)

ggplot(richtable.combined,aes(Shannon.x, Shannon.y, col=cluster.x))+
  geom_point()

ggplot(richtable.combined,aes(agecat.x, Shannon.x, col=colorcat.x))+
  geom_boxplot()

richtable.combined2<-richtable.combined[richtable.combined$Shannon.y!=0,]

reg.sh=lme(Shannon.x~ Shannon.y, random=~1|cluster.x, data = richtable.combined2)
sh.pval<-as.data.frame(anova.lme(reg.sh))#type 3 sums of squares

e<-coef(summary(reg.sh))
e

pi.e<-allEffects(reg.sh,xlevels=20)
pi.e.df<-as.data.frame(pi.e[[1]]) #0.95 default confidence

quartz()
ggplot(pi.e.df,aes(Shannon.y, fit)) + 
  geom_line()+
  geom_point(data=richtable.combined2, aes(Shannon.y, Shannon.x),alpha = 0.6,size=2) +
  geom_ribbon(data=pi.e.df,aes(Shannon.y, ymin = lower, ymax = upper), alpha = .2,show.legend=FALSE,colour = NA)+
  xlab("Bacterial Diversity")+
  ylab("Fungal Diversity")

####test response to ratio ####
richtable.combined$ratio<-richtable.combined$Shannon.y/richtable.combined$Shannon.x
  #higher value = higher bacteria and lower fungi

