#####Script for microbiome analyses across cohorts
##Usage: Figure S1, Figure 2A-D

##load packages

library(phyloseq)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(dplyr)
library(vegan)
library(ggpubr)


'%not in%' <- function(x,y)!('%in%'(x,y))

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}



##########################################################################
##load samples from STL cohort (n=644): this part is identical to previous script 

sub_metaphlan <- read.csv('UTI_metaphlan3_merged.txt', sep='\t') 

rownames(sub_metaphlan) <- sub_metaphlan$clade_name
sub_metaphlan$clade_name <- NULL
colnames(sub_metaphlan) <- gsub('\\.', '-', colnames(sub_metaphlan))

#take only rows that include s__ (species level annotation)
sub_metaphlan_species <- sub_metaphlan[grepl('s__', rownames(sub_metaphlan)), ]
#now get rid of rows that also include the 'taxa' level
sub_metaphlan_species <- sub_metaphlan_species[!grepl('t__', rownames(sub_metaphlan_species)), ]
sub_metaphlan_species$NCBI_tax_id <- NULL
samples <- colnames(sub_metaphlan_species)


#MAKE TAXA TABLE FOR PHYLOSEQ FROM METAPHLAN TAXA NAMES#
sub_species.names <- data.frame(Names = rownames(sub_metaphlan_species))
sub_species.names <- data.frame(do.call('rbind', strsplit(as.character(sub_species.names$Names),'|',fixed=TRUE)))
rownames(sub_species.names) <- rownames(sub_metaphlan_species)
colnames(sub_species.names) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
sub_species.names <- as.matrix(sub_species.names)

#load metadata
meta <- read.csv('master_metadata.csv')

rownames(meta) <- meta$Sample_ID
meta$status_fixed <- as.factor(as.character(meta$status_fixed))

#make phyloseq object
stl_ps <- phyloseq(otu_table(sub_metaphlan_species, taxa_are_rows = TRUE),
                   sample_data(meta),
                   tax_table(sub_species.names))

##########################################################
#filter stl_ps to only have E1-S1 sample 

wu_meta <- as.matrix(sample_data(stl_ps))
wu_meta <- as.data.frame(wu_meta)
#View(wu_meta)

patients <- unique(wu_meta$Patient_ID) #106 

wu_meta$Sample_ID <- as.character(wu_meta$Sample_ID)
wu_meta$end <- substrRight(wu_meta$Sample_ID,5)

keepers <- wu_meta[wu_meta$end=='E1-S1',] #96 

stl_ps_filt <- subset_samples(stl_ps, sample_names(stl_ps)%in%keepers$Sample_ID)

ntaxa(stl_ps_filt)  #527 
nsamples(stl_ps_filt)  #96 


######################################################################################
#Broad: load all data first, subset to just 1st sample per individual (total)

brd_metaphlan <- read.csv('Broad_31_metaphlan3.csv', sep=',') #new data from metaphlan3
rownames(brd_metaphlan) <- brd_metaphlan$X
brd_metaphlan$X <- NULL


#FIND ONLY ROWS AT THE SPECIES LEVEL
#take only rows that include s__ (species level annotation)
brd_metaphlan_species <- brd_metaphlan[grepl('s__', rownames(brd_metaphlan)), ]
#now get rid of rows that also include the 'taxa' level
brd_metaphlan_species <- brd_metaphlan_species[!grepl('t__', rownames(brd_metaphlan_species)), ]
brd_metaphlan_species$NCBI_tax_id <- NULL
samples <- colnames(brd_metaphlan_species)

#MAKE TAXA TABLE FOR PHYLOSEQ FROM METAPHLAN TAXA NAMES#
brd_species.names <- data.frame(Names = rownames(brd_metaphlan_species))
brd_species.names <- data.frame(do.call('rbind', strsplit(as.character(brd_species.names$Names),'|',fixed=TRUE)))
rownames(brd_species.names) <- rownames(brd_metaphlan_species)
colnames(brd_species.names) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
brd_species.names <- as.matrix(brd_species.names)

#load metadata
brd_meta <- read.csv('broad_metadata.csv')
rownames(brd_meta) <- brd_meta$Sample_ID
brd_meta$status_fixed <- as.factor(as.character(brd_meta$status_fixed))

#build phyloseq object
broad_ps <- phyloseq(otu_table(brd_metaphlan_species, taxa_are_rows = TRUE),
                     sample_data(brd_meta),
                     tax_table(brd_species.names))

ntaxa(broad_ps)  #342 
nsamples(broad_ps)  #31



###############################################################################
#Healthy Humans

HH_metaphlan <- read.csv('HH_metaphlan3_merged.txt', sep='\t') #new data from metaphlan3
rownames(HH_metaphlan) <- HH_metaphlan$clade_name
HH_metaphlan$clade_name <- NULL
colnames(HH_metaphlan) <- gsub('\\.', '-', colnames(HH_metaphlan))
colnames(HH_metaphlan) <- gsub('X', '\\', colnames(HH_metaphlan))

#FIND ONLY ROWS AT THE SPECIES LEVEL
#take only rows that include s__ (species level annotation)
HH_metaphlan_species <- HH_metaphlan[grepl('s__', rownames(HH_metaphlan)), ]
#now get rid of rows that also include the 'taxa' level
HH_metaphlan_species <- HH_metaphlan_species[!grepl('t__', rownames(HH_metaphlan_species)), ]
HH_metaphlan_species$NCBI_tax_id <- NULL
samples <- colnames(HH_metaphlan_species)

#MAKE TAXA TABLE FOR PHYLOSEQ FROM METAPHLAN TAXA NAMES#
HH_species.names <- data.frame(Names = rownames(HH_metaphlan_species))
HH_species.names <- data.frame(do.call('rbind', strsplit(as.character(HH_species.names$Names),'|',fixed=TRUE)))
rownames(HH_species.names) <- rownames(HH_metaphlan_species)
colnames(HH_species.names) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
HH_species.names <- as.matrix(HH_species.names)


HH_meta <- read.csv('HH_metadata.csv')
rownames(HH_meta) <- HH_meta$Sample_ID
HH_meta$status_fixed <- as.factor(as.character(HH_meta$status_fixed))


HH_ps <- phyloseq(otu_table(HH_metaphlan_species, taxa_are_rows = TRUE),
                  sample_data(HH_meta),
                  tax_table(HH_species.names))

ntaxa(HH_ps)  #312
nsamples(HH_ps)  #20

################################################################################################
###now we have three phyloseq objects: stl_ps_filt, broad_ps_filt, HH_ps
##let's merge them! 

GP3 = merge_phyloseq(stl_ps_filt, broad_ps, HH_ps) 
ntaxa(GP3) #707
nsamples(GP3) #147

GP3.df <- psmelt(GP3)

#filter for 10% presence
combined_ps = filter_taxa(GP3, function(x) sum(x > 0) > (0.1*length(x)), TRUE) ## this makes colSums not 100 because I didn't add in "Others" as taxon 

ntaxa(combined_ps) #133
nsamples(combined_ps) #147

###load metadata 
meta <- read.csv('stl_broad_HH_metadata.csv')
rownames(meta) <- meta$Sample_ID

length(rownames(meta)) #695

meta_filt <- meta[meta$Sample_ID %in% as.data.frame(sample_data(combined_ps))$Sample_ID,]
#View(meta_filt)

length(rownames(meta_filt)) #147
sample_data(combined_ps) <- meta_filt
melt <- psmelt(combined_ps) 

species.df<-dcast(melt[c("Sample","Species","Abundance")],Sample~Species,fill=0.0, aggregate = sum)
rownames(species.df) <- species.df$Sample
species.df$Sample <- NULL


#ALPHA DIVERSITY 
alpha.diversity<-data.frame(
  shannon=diversity(species.df,index="shannon"),
  simpson=diversity(species.df,index="simpson"),
  richness=rowSums(species.df>0.0)
)

##add metadata to alpha.diversity df 
n <- length(rownames(alpha.diversity))
alpha.diversity$cohort <- c('.')
alpha.diversity$status_fixed <- c('.')
alpha.diversity$sex <- c('.')
alpha.diversity$age <- c('.')

for (i in 1:n){
  temp <- as.character(rownames(alpha.diversity)[i])
  temp2 <- as.character(meta_filt[meta_filt$Sample_ID==temp,]$Study)
  temp3 <- as.character(meta_filt[meta_filt$Sample_ID==temp,]$status_fixed)
  temp4 <- as.character(meta_filt[meta_filt$Sample_ID==temp,]$sex)
  temp5 <- as.character(meta_filt[meta_filt$Sample_ID==temp,]$age)
  alpha.diversity$cohort[i] <- temp2
  alpha.diversity$status_fixed[i] <- temp3
  alpha.diversity$sex[i] <- temp4
  alpha.diversity$age[i] <- temp5
}

###########
#group by status... N: non-rUTI, R: rUTI, H: Healthy
##merging non-rUTI and rUTI as "UTI" for cross-cohort comparison
alpha.diversity$status_fixed <- gsub('N', 'R', alpha.diversity$status_fixed)


##Figure 2A: richness by cohort 


ggplot(aes(x=status_fixed,y=richness), 
       data=alpha.diversity)+ stat_compare_means(method="kruskal")+
  geom_boxplot(outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1)) +theme_bw()+
  geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
             pch=21, aes(fill=factor(cohort))) +
  scale_fill_manual(values = c("UMB" = "blue",
                               "STL" = "red",
                               "HH" = "yellow"))


#just 2 categories : UTI vs Healthy

combined_2groups <- combined_ps
sample_data(combined_2groups)$status_fixed <- gsub('N', 'R', sample_data(combined_2groups)$status_fixed)


sample_data(combined_2groups)$status_fixed <- as.factor(as.character(sample_data(combined_2groups)$status_fixed))
sample_data(combined_2groups)$Study <- as.factor(as.character(sample_data(combined_2groups)$Study))

combined_bray2 <- phyloseq::distance(combined_2groups, method = "bray")

adonis2(combined_bray2 ~ as.factor(Study)+ as.factor(status_fixed), data = data.frame(sample_data(combined_2groups)))
adonis2(combined_bray2 ~ as.factor(status_fixed)+ as.factor(Study), data = data.frame(sample_data(combined_2groups)))
adonis2(combined_bray2 ~ as.factor(Study)+ as.numeric(age)+as.factor(status_fixed), data = data.frame(sample_data(combined_2groups)))


#adonis(combined_bray2 ~ as.factor(Study)+ as.factor(status_fixed), data = data.frame(sample_data(combined_2groups)))


##visualize! 

##PCoA
GP.ord <- ordinate(combined_2groups, "PCoA", "bray")
p1 = plot_ordination(combined_2groups, GP.ord, color="status_fixed")
p1+stat_ellipse(type = "norm", linetype=2)+theme_bw()

ps.bray <- phyloseq::distance(combined_2groups, method='bray')
ps.cap.ord <- ordinate(combined_2groups, 
                       method = 'CAP',
                       distance = ps.bray,
                       formula = distance ~ Study+status_fixed,
                       na.action = na.omit)

ps.cap.ord



#cap: Figure 2B
cap_plotb <- plot_ordination(physeq = combined_2groups, ordination = ps.cap.ord, color = 'Study')+
  geom_point(size=2)+
  stat_ellipse(type="norm",linetype=2)+
  theme_classic()

cap_plotb


#cap: Figure 2C
cap_plot <- plot_ordination(physeq = combined_2groups, ordination = ps.cap.ord, color = 'status_fixed')+
  geom_point(size=2)+
  stat_ellipse(type="norm",linetype=2)+
  theme_classic()

cap_plot


#Access CAP axis. 
cap.df <- cap_plot$data

#vegan's anova for significance of the constraints, by term
#https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/anova.cca

cap.anova <- anova(ps.cap.ord, by='term') 
cap.anova

cap.anova.mar <- anova(ps.cap.ord, by='mar') 
cap.anova.axis <- anova(ps.cap.ord, by='axis') 



######################### additional testing ######################
## "match" cohorts for sample sizes, sex, and age

stl <- subset_samples(combined_2groups, Study=="STL")
nsamples(stl) #96
stl_young <- subset_samples(stl, age<42)

nsamples(stl_young) #20
unique(sample_data(stl_young)$sex) ##all are F 


combined_2groups_filt <- subset_samples(combined_2groups, Study!="STL")
combined_2groups_filt <- subset_samples(combined_2groups_filt, sex=="F")

nsamples(combined_2groups_filt) #41 (10 HH, 31 UMB)
temp <- as.data.frame(sample_data(combined_2groups_filt))
temp$age
mean(temp$age) #30.29
mean(sample_data(stl_young)$age) #32.35 



#combine phyloseq objects 

new_subset <- merge_phyloseq(stl_young, combined_2groups_filt)
nsamples(new_subset) ##61 samples 

temp2 <- as.data.frame(sample_data(new_subset))
length(rownames(temp2[temp2$status_fixed=='H',])) #26 healthy, 35 UTI



subset_bray2 <- phyloseq::distance(new_subset, method = "bray")

adonis2(subset_bray2 ~ as.factor(Study)+as.factor(status_fixed), data = data.frame(sample_data(new_subset)))
adonis2(subset_bray2 ~  as.factor(status_fixed), data = data.frame(sample_data(new_subset)))


##PCoA
GP.ord <- ordinate(new_subset, "PCoA", "bray")
p1 = plot_ordination(new_subset, GP.ord, color="status_fixed")
p1+stat_ellipse(type = "norm", linetype=2)+theme_bw()

ps.bray <- phyloseq::distance(new_subset, method='bray')
ps.cap.ord <- ordinate(new_subset, 
                       method = 'CAP',
                       distance = ps.bray,
                       formula = distance ~ Study+status_fixed ,
                       na.action = na.omit)

ps.cap.ord

##Figure S1A
cap_plotb <- plot_ordination(physeq = new_subset, ordination = ps.cap.ord, color = 'Study')+
  geom_point(size=2)+
  stat_ellipse(type="norm",linetype=2)+
  theme_classic()

cap_plotb


##Figure S1B
cap_plot <- plot_ordination(physeq = new_subset, ordination = ps.cap.ord, color = 'status_fixed')+
  geom_point(size=2)+
  stat_ellipse(type="norm",linetype=2)+
  theme_classic()

cap_plot



###########maaslin2
library("Maaslin2")

meta <- as.data.frame(sample_data(combined_2groups))
meta$status_fixed <- as.factor(as.character(meta$status_fixed))

## can replace family.df for species.df or genus.df for different level analyses
fit_data = Maaslin2(
  input_data = family.df, 
  input_metadata = meta, 
  transform = "AST",
  normalization = "NONE",
  standardize = FALSE,
  output = "./Maaslin2_HH_family", 
  fixed_effects = c('status_fixed'), random_effects= c('Study'))


##Plotting maaslin output: Figure 2D
input <- read.csv('genus_plot_input.csv')
n <- length(rownames(input))
length(rownames(genus.df)) #147 
input$relative_abundance <- c('.')
input$shape <- c('.')

for (i in 1:n){
  taxon <- as.character(input$feature[i])
  num <- sum(genus.df[,taxon])/147
  
  input$relative_abundance[i] <- num
  
  
  coef <- as.numeric(as.character(input$coef[i]))
  if (coef >0 ){
    input$shape[i] <-'2'
  }
  else{
    input$shape[i] <-'6'
  }
  
}

input$feature <- gsub('g__', '\\', input$feature)

ggplot(input, aes(x=FDR, y=as.numeric(relative_abundance), label=feature))+geom_point(shape=input$shape, size=3)+geom_text(hjust=0, vjust=0)+
  theme_bw()

input$direction <- input$FDR
input[input$shape=='6',]$direction <- -1 * input[input$shape=='6',]$direction

ggplot(input, aes(group=as.factor(feature)))+
  geom_segment(aes(x = coef,y= reorder(feature, desc(feature)), xend = 0, yend = as.factor(feature)), data = input)+
  xlim(-0.1, 0.22)+theme_bw()


