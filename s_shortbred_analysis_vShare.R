
#######Script for interpreting shortbred ARG burden predictions
##Usage: Figures 2E, 2F


library(phyloseq)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(dplyr)
library(vegan)


'%not in%' <- function(x,y)!('%in%'(x,y))

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


arg <- read.csv('shortbred_combined.txt', sep='\t')


length(rownames(arg))
length(unique(arg$Sample)) 

#"count" is the normalized RPKM (reads per kilobase of reference sequence per million sample reads)

#want to first generate a data frame like species.df 
args.df <-dcast(arg[c("Sample","Family","Count")],Sample~Family,fill=0.0, aggregate = sum)

head(args.df)
markers <- colnames(args.df)
keepers <- c()

n <- length(markers)
for (i in 1:n){
  temp <- substr(markers[i],1,1)
  if (temp!='D'){
    keepers <- c(keepers, i)
  }
  
  
}

length(keepers) #2416

args.df.filt <- args.df[,keepers]

rownames(args.df.filt) <- args.df.filt$Sample  

meta <- read.csv('stl_broad_HH_metadata.csv')
meta_filt <- meta[meta$Sample_ID %in% rownames(arg.alpha),]
length(rownames(meta_filt)) #146
args.df.filt <- args.df.filt[args.df.filt$Sample %in% meta_filt$Sample_ID,]
args.df.filt$Sample <- NULL

#ALPHA DIVERSITY 
arg.alpha<-data.frame(
  shannon=diversity(args.df.filt,index="shannon"),
  simpson=diversity(args.df.filt,index="simpson"),
  richness=rowSums(args.df.filt>0.0),
  RPKM=rowSums(args.df.filt)
)

##add metadata to alpha.diversity df 
n <- length(rownames(arg.alpha))
arg.alpha$cohort <- c('.')
arg.alpha$status_fixed <- c('.')


for (i in 1:n){
  temp <- as.character(rownames(arg.alpha)[i])
  temp2 <- as.character(meta_filt[meta_filt$Sample_ID==temp,]$Study)
  temp3 <- as.character(meta_filt[meta_filt$Sample_ID==temp,]$status_fixed)
  arg.alpha$cohort[i] <- temp2
  arg.alpha$status_fixed[i] <- temp3
  
}

arg.alpha$status_2 <- arg.alpha$status_fixed
arg.alpha$status_2 <- gsub('R', 'N', arg.alpha$status_2)

arg.alpha.filt <- arg.alpha[arg.alpha$cohort!='UMB',]

###########

library(ggpubr)


##Figure 2E 
ggplot(aes(x=status_2,y=log(RPKM)),  
       data=arg.alpha)+ stat_compare_means(method="kruskal")+
  geom_boxplot(outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1)) +theme_bw()+
  geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
             pch=21, aes(fill=factor(cohort))) +
  scale_fill_manual(values = c("UMB" = "blue",
                               "STL" = "red",
                               "HH" = "yellow"))

##Figure 2F
ggplot(aes(x=status_2,y=richness), 
       data=arg.alpha)+ stat_compare_means(method="kruskal")+
  geom_boxplot(outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1)) +theme_bw()+
  geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
             pch=21, aes(fill=factor(cohort))) +
  scale_fill_manual(values = c("UMB" = "blue",
                               "STL" = "red",
                               "HH" = "yellow"))




