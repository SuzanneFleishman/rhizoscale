#This code computes shannon diversity and tests it using a mixed model
#project repository: https://github.com/SuzanneFleishman/rhizoscale
#publication code created for: (under review, will be updated)

#Code compiled by Suzanne Fleishman; this version saved May 11, 2021
#Based on Borcard, D., Gillet, F., and Legendre, P. (2018). Numerical Ecology with R (Springer). (https://doi.org/10.1007/978-1-4419-7976-6)


#Inputs are Phyloseq objects created in "1PhyloseqObjectCreation"

#To complete this code, run it separately for all phyloseq objects; select at L44
#code needs attention at each step



#### Setup ####

setwd("xxx")


### Clear workspace ###
rm(list=ls())

### Packages ###
#Corncob Tutorial and more info: https://rdrr.io/github/bryandmartin/corncob/f/vignettes/corncob-intro.Rmd
#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")
# devtools::install_github("bryandmartin/corncob")
### Load Packages ###
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("corncob","phyloseq", "ggplot2","DESeq2","microbiomeSeq","metagMisc","dplyr","fastDummies","stringr","vegan","tidyr")
ipak(packages)

### Load Phyloseq Objects ###
load("xxx/PS_rhizoscale.Rdata")

####Select PS object ####

ps = ps.its


#### Calculate bray curtis dissimilarity ####
BC <- phyloseq::distance(otu_table(ps), "bray")

####calculate beta dispersion ####

meta<-as.data.frame(ps@sam_data) #metadata dataframe

beta.gen<-betadisper(BC, meta$age) #betadispersion based on age

set.seed(1)
permutest(beta.gen,strata=cluster) #permutest, stratified by cluster

####visualize####

quartz()
g<-boxplot(beta.gen)

quartz()
plot(beta.gen)

####create plot####

#create table of means and standard errors
  mod<-beta.gen
  
  ## Compute mean distance to centroid per group
  me<-tapply(mod$distances, meta$age, mean)
  
  ## variance
  sd<-tapply(mod$distances, meta$age, var)
  
  #run first for bacteria
  bac<-as.data.frame(rbind(me,sd))
    bac<-as.data.frame(t(bac))
    bac$sp<-"bac"
  
  #then run for fungi
  fung<-as.data.frame(rbind(me,sd))
    fung<-as.data.frame(t(fung))
    fung$sp<-"fung"

#create dataframe for bac and fungi, then  plot

tab<-rbind(bac,fung)
tab$age<-c("old","young","old","young")


p<- ggplot(tab, aes(x=sp, y=me, fill=age)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=me-sd, ymax=me+sd), width=.2,
                position=position_dodge(.9)) +
  ylim(0, 1.25)

quartz()
p

