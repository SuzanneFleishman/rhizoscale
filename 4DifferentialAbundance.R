#This code computes differential abundance and tests it using the corncob package
#project repository: https://github.com/SuzanneFleishman/rhizoscale
#publication code created for: (under review, will be updated)

#Code compiled by Suzanne Fleishman; this version saved May 11, 2021
#Corncob package citation: Bryan D. Martin, Daniela Witten, and Amy D. Willis. (2020). Modeling microbial abundances and dysbiosis with beta-binomial regression. Annals of Applied Statistics, Volume 14, Number 1, pages 94-115.DOI: 10.1214/19-AOAS1283
#Corncob Tutorial and more info: https://rdrr.io/github/bryandmartin/corncob/f/vignettes/corncob-intro.Rmd

#Inputs are Phyloseq objects created in "1PhyloseqObjectCreation"

#For this study, this code was only run on the filtered datasets (ps.16s.10 and ps.its.10) in order to reduce false significance due to high numbers of rare taxa

#code needs attention at each step

#### Setup ####

setwd("xxx")


### Clear workspace ###
rm(list=ls())

### Packages ###

### Load Packages ###
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("corncob","phyloseq", "ggplot2","DESeq2","microbiomeSeq","metagMisc","dplyr","fastDummies","stringr")
ipak(packages)

### Prep Phyloseq Objects ###
load("xxx/PS_rhizoscale.Rdata")

####CORNCOB ####

ps.sub = ps.its.10

####
#for this study, only MEMs significant in "2VariationPartitionAnalysis" were used to prevent overfitting

ps.sub@sam_data

set.seed(1)
da_rootchar<- differentialTest(formula = ~ x+y+colorcat+agecat+type+len+a+b+c+d+e+f+g+h+MEM2+MEM6, #the DA formula - if you have a controlling factor it must be included. 
                                phi.formula = ~ x+y+colorcat+agecat+type+len+a+b+c+d+e+f+g+h+MEM2+MEM6, #the DV formula - if you have a controlling factor it must be included
                                formula_null = ~ 1, # DA controlling factors
                                phi.formula_null = ~x+y+colorcat+agecat+type+len+a+b+c+d+e+f+g+h+MEM2+MEM6, #DV controlling factors
                                test = "Wald", boot = FALSE, #wald test is "standard"
                                data = ps.sub, #PS object
                                fdr_cutoff = 0.05) #pval false discovery value
#explore significant taxa
l<-as.data.frame(da_rootchar$significant_taxa)
models<-da_rootchar$significant_models

p<-plot(da_rootchar)

#visualize
quartz()
p


d<-as.data.frame(p$data)

####create dataframe for improved plot####

#first, create a taxonomic table of information wanted for the plot
all.tax<-as.data.frame(ps.sub@tax_table@.Data)
phy.gen.all<-cbind(rownames(all.tax), all.tax$sequence, all.tax$Phylum,all.tax$Family,all.tax$Genus)
phy.gen<-unique(phy.gen.all)
colnames(phy.gen)<-c("asv","seq", "Phylum", "Family","Genus")

#then, extract the taxa names in the same format to plot
d$seq2<-str_extract(d$taxa, "[^_]+")

#merge the plots of significant factors with taxonomic information
joined_df <- merge(d, phy.gen, by.x = "seq2", 
                   by.y = "seq", all.x = TRUE, all.y = FALSE)

#create new phylum column based on significance
joined_df %>% arrange(Phylum)
joined_df$sigNA<-joined_df$Phylum
joined_df$sigNA[(joined_df$xmax*joined_df$xmin)<0]<-"NS" #anything that overlaps "0" is not significant


#reformat for plotting
joined_df_up_order<-joined_df
joined_df_up_order$Genus <- factor(joined_df_up_order$Genus, levels=unique(joined_df_up_order$Genus ))
joined_df_up_order<-joined_df[order(joined_df$Genus),]
joined_df_up_order$Genus <- factor(joined_df_up_order$Genus, levels=unique(joined_df_up_order$Genus ))
joined_df_up_order$taxphy<-paste("(",joined_df_up_order$Genus,")",joined_df_up_order$asv)



####Remove Taxa without significant factors####
#corncob will still select models that are significant, even if just the intercept is significant
#this is a loop code that filters out models for taxa that only have a significant intercept

#create vector of taxa that only have significant intercepts
nonsig<-c("F_asv63","F_asv1")


joined_df_up_order2<-joined_df_up_order


for (i in 1:length(nonsig)) { joined_df_up_order2<-subset(joined_df_up_order2, !(joined_df_up_order2$asv==nonsig[i]))

}
joined_df_up_order2



####plot####

#ensure the color "lightgrey" matches up with "NS" in the legend
#other colors should map to the Phylum

g =ggplot(data=joined_df_up_order2,
    aes(x = taxphy, y = x, ymin = xmin, ymax = xmax ))+
  geom_bar(stat="identity",aes(fill=sigNA),width=0.8)+
  geom_errorbar(aes(x=taxphy, ymin=xmin, ymax=xmax),width=0,size=.2)+
  xlab('Genus')+ ylab("Log2fold Change")+
  facet_wrap(~variable,strip.position="top",ncol=8,scales = "free_x") +
  scale_fill_manual(values = c(  "lightgrey","darkgreen","#F0E442",
                                 "#0072B2",  "#D55E00", "#CC79A7",
                                 "#E69F00", "#009E73", "#F0E442","#009E73")) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"))+
  coord_flip()

quartz()
g

