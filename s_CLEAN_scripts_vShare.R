#####Script for microbiome analyses within cohort
##Usage: Figures S2, 3A-3E


##load packages

library(phyloseq)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(vegan)
library(dplyr)
library(ggpubr)
library(Maaslin2)

##note: 'agecat' in metadata file refers to 3 groups. group A: 18-64; groupB: 65-79; group C: ≥80


###LOAD IN DATA###--------------------------------------------------------------

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


###MAKE PHYLOSEQ OBJECT###------------------------------------------------------

#load metadata
meta <- read.csv('master_metadata.csv')

rownames(meta) <- meta$Sample_ID
meta$status_fixed <- as.factor(as.character(meta$status_fixed))

#make phyloseq object
sub_ps <- phyloseq(otu_table(sub_metaphlan_species, taxa_are_rows = TRUE),
                   sample_data(meta),
                   tax_table(sub_species.names))

##filter by at least 10% presence 
GP = filter_taxa(sub_ps, function(x) sum(x > 0) > (0.1*length(x)), TRUE) ## this makes colSums not 100 because I didn't add in "Others" as taxon 

sub_ps <- GP
sub_ps.df <- psmelt(sub_ps) #melt phyloseq object into a data frame, useful for ggplots
ntaxa(sub_ps)  #143  
nsamples(sub_ps) #644


###EPISODE 1 ANALYSES###------------------------------------------------------
##Filter for E1 only 
sub_ps_e1 = prune_samples(sample_data(sub_ps)$episode=='E1', sub_ps)
sub_ps_e1.df <- psmelt(sub_ps_e1) 

species_e1.df<-dcast(sub_ps_e1.df[c("Sample","Species","Abundance")],Sample~Species,fill=0.0, aggregate = sum)
rownames(species_e1.df) <- species_e1.df$Sample
species_e1.df$Sample <- NULL

nsamples(sub_ps_e1) #480
ntaxa(sub_ps_e1) #143

temp <- as.data.frame(sample_data(sub_ps_e1))
View(temp)

#ALPHA DIVERSITY 
alpha.diversity_E1<-data.frame(
  shannon=diversity(species_e1.df,index="shannon"),
  simpson=diversity(species_e1.df,index="simpson"),
  richness=rowSums(species_e1.df>0.0)
)

#add on E1 metadata to alpha.diversity_E1
meta_E1 <- meta[meta$Sample_ID %in% rownames(species_e1.df),]


alpha.diversity_E1$richness<-as.integer(alpha.diversity_E1$richness)
alpha.diversity_E1$sample<- rownames(alpha.diversity_E1)

n <- length(rownames(alpha.diversity_E1)) #480
alpha.diversity_E1$status <- c('')
alpha.diversity_E1$episode <- c('')
alpha.diversity_E1$Gcol <- c('')
alpha.diversity_E1$Ucol <- c('')
alpha.diversity_E1$timebin <- c('')
alpha.diversity_E1$drug <- c('')
alpha.diversity_E1$coltype <- c('')

for (i in 1:n){
  temp <- as.character(alpha.diversity_E1$sample[i])
  
  if (temp %in% meta_E1$Sample_ID){
    
    tempdf <- meta_E1[meta_E1$Sample_ID==temp,]
    alpha.diversity_E1$status[i] <- as.character(tempdf$status_fixed)
    alpha.diversity_E1$episode[i] <- as.character(tempdf$episode)
    alpha.diversity_E1$Gcol[i] <- as.character(tempdf$Gcol)
    alpha.diversity_E1$Ucol[i] <- as.character(tempdf$Ucol)
    alpha.diversity_E1$timebin[i] <- as.character(tempdf$timebin)
    alpha.diversity_E1$drug[i] <- as.character(tempdf$treatment)
    alpha.diversity_E1$coltype[i] <- as.character(tempdf$coltype)
  }
  
}

alpha.diversity_E1$coltype <- gsub('SU', 'B', alpha.diversity_E1$coltype) ## stool (S) + urine (U)  = both (B)


##plot richness
alpha.diversity_E1$status <- as.factor(alpha.diversity_E1$status)
alpha.diversity_E1$sample <- rownames(alpha.diversity_E1)


##Figure S2A: rUTI/non-rUTI vs richness
ggplot(aes(x=status, y=richness,color=status), 
       data=alpha.diversity_E1)+
  geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_point(position=position_jitterdodge(.5), size=2, alpha=0.5)+
  stat_compare_means(method = "kruskal.test", label =  "p.signif")+theme_bw()


kruskal.test(richness ~ status, data=alpha.diversity_E1) #p-val: 0.3745


##Figure S2B: rUTI status vs shannon diversity
ggplot(aes(x=status, y=shannon,color=status), 
       data=alpha.diversity_E1)+
  geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_point(position=position_jitterdodge(.5), size=2, alpha=0.5)+
  stat_compare_means(method = "kruskal.test", label =  "p.signif")+theme_bw()


kruskal.test(shannon ~ status, data=alpha.diversity_E1) #p-val: 0.2362


###ordination###-----------------

GP.ord <- ordinate(sub_ps_e1, "PCoA", "bray")

##Figure S2C: rUTI vs UTI E1 samples, n=480 
ps = plot_ordination(sub_ps_e1, GP.ord, color="status_fixed")
ps + 
  stat_ellipse(type = "norm", linetype = 2) + #normal distribution 
  theme_bw()



E1_bray <- phyloseq::distance(sub_ps_e1, method = "bray")

## patient ID explains greater variability than rUTI status
adonis2(E1_bray ~ Patient_ID + as.factor(status_fixed), data = data.frame(sample_data(sub_ps_e1))) 
adonis2(E1_bray ~ Patient_ID + as.factor(treatment)+as.factor(status_fixed), data = data.frame(sample_data(sub_ps_e1))) 


########AGGREGATE by patient AVERAGE microbiome profile ---------------
## function to get average values 
temp <-dcast(sub_ps_e1.df[c("Patient_ID","OTU","Abundance")],OTU~Patient_ID,fill=0.0, fun.aggregate = mean)


rownames(temp) <- temp$OTU
temp$OTU <- NULL
temp <- as.matrix(temp)
temp2 <- otu_table(temp, taxa_are_rows=T)

merged_sub_ps= merge_samples(sub_ps_e1, "Patient_ID", fun=mean) ##OTU tables are SUMMED not averaged 
otu_table(merged_sub_ps) <- temp2   ##this is why we define otu table separately 
#write.csv(temp2, 'UTI_metaphlan3_merged_by_patient.csv')

meta_e1 <- meta[meta$episode=='E1',]
meta_pat <- meta_e1[!duplicated(meta_e1[,c('Patient_ID')]),]
length(rownames(meta_pat)) #106
View(meta_pat)

rownames(meta_pat) <- meta_pat$Patient_ID
meta_pat$Sample_ID <- NULL

meta_pat$status_fixed <- as.factor(as.character(meta_pat$status_fixed))
sample_data(merged_sub_ps) <- meta_pat


ntaxa(merged_sub_ps) #143
nsamples(merged_sub_ps) #106 patients


#important for maintaining factor status 
sample_data(merged_sub_ps)$status_fixed <- as.factor(sample_data(merged_sub_ps)$status_fixed)


##rUTI status is non-signif 
merged_bray<- phyloseq::distance(merged_sub_ps, method = "bray")
adonis2(merged_bray ~ agecat+ treatment+ status_fixed, data = data.frame(sample_data(merged_sub_ps)))
adonis2(merged_bray ~ status_fixed, data = data.frame(sample_data(merged_sub_ps)))

##Figure S2D: average profile and rUTI status
GP.ord <- ordinate(merged_sub_ps, "PCoA", "bray")
p1 = plot_ordination(merged_sub_ps, GP.ord, color="status_fixed")
p1+stat_ellipse(type = "norm", linetype=2)+theme_bw()


######## Timebin-specific differences (post-abx microbiome recovery)---------------

alpha.diversity_E1 <- na.omit(alpha.diversity_E1) ## some samples marked NA for timebin due to not fitting into category/temporal overlap

##Supp Figure 3A: richness by timebin
ggplot(aes(x=timebin, y=richness, group=timebin), data=alpha.diversity_E1)+
  geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

compare_means(richness ~ timebin,  data = alpha.diversity_E1, ref.group = "5",
              method = "wilcox", p.adjust.method = "BH")


##Supp Fig 3B: line graph of drugs and timebin 
e1_filt <- alpha.diversity_E1[alpha.diversity_E1$drug!="" ,]

ggplot(e1_filt, aes(x=timebin, y=richness, group=drug, color=drug))+geom_point()+ geom_smooth(method = "lm", fill = NA)
kruskal.test(shannon ~ drug, data=e1_filt)


##Kruskal test at each timepoint to see if differences between abx are significant 
e1_filt1 <- e1_filt[e1_filt$timebin=='1',] ##signif diff at timebin 1
kruskal.test(shannon ~ drug, data=e1_filt1)

e1_filt2 <- e1_filt[e1_filt$timebin=='2',] ## signif at timebin 2
kruskal.test(shannon ~ drug, data=e1_filt2)

e1_filt3 <- e1_filt[e1_filt$timebin=='3',] ## no diff between drugs by timebin 3
kruskal.test(shannon ~ drug, data=e1_filt3)

e1_filt4 <- e1_filt[e1_filt$timebin=='4',]
kruskal.test(shannon ~ drug, data=e1_filt4)

e1_filt5 <- e1_filt[e1_filt$timebin=='5',]
kruskal.test(shannon ~ drug, data=e1_filt5)


##post-hoc Dunn test for kruskal wallis test 

library(FSA)

DT = dunnTest(shannon ~ drug,
              data=e1_filt2,
              method="bh")  

DT2 <- as.data.frame(DT[[2]])


#DT2_filt <- DT2[DT2$P.adj<0.05,]
#View(DT2_filt)

########################################################################
#########aggreagate by patient average within S4 and S5---------

###FILTER BY S4 or S5 ONLY 
select <- c('4', '5')
sub_ps_s45 = prune_samples(sample_data(sub_ps_e1)$timepoint %in% select, sub_ps_e1)

nsamples(sub_ps_s45) #114 
View(sample_data(sub_ps_s45))

melt <- psmelt(sub_ps_s45)

## proper function to get average values per patient
temp <-dcast(melt[c("Patient_ID","OTU","Abundance")],OTU~Patient_ID,fill=0.0, fun.aggregate = mean)

rownames(temp) <- temp$OTU
temp$OTU <- NULL
temp <- as.matrix(temp)
temp2 <- otu_table(temp, taxa_are_rows=T)


merged_sub_ps_s45= merge_samples(sub_ps_s45, "Patient_ID", fun=mean) ##OTU tables are SUMMED not averaged 
otu_table(merged_sub_ps_s45) <- temp2   ##this is why we define otu table separately 
nsamples(merged_sub_ps_s45) #96
sample_names(merged_sub_ps_s45)
temp <- as.data.frame(sample_data(merged_sub_ps_s45))

patnames <- rownames(temp)
meta_s45 <- meta_pat[meta_pat$Patient_ID %in% patnames,]
length(rownames(meta_s45)) #96 

sample_data(merged_sub_ps_s45) <- meta_s45 ##timebin, timepoint are inaccurate here since they were taken from meta_pat

ntaxa(merged_sub_ps_s45) #143
nsamples(merged_sub_ps_s45) #96 patients with S4, S5: 63non- Ucol, 33 Ucol


## optional: filter out males 
#f_merged_sub_ps_s45 = prune_samples(sample_data(merged_sub_ps_s45)$Sex=='F', merged_sub_ps_s45)
#nsamples(f_merged_sub_ps_s45) #89 


#important for maintaining factor status!! 
sample_data(merged_sub_ps_s45)$status_fixed <- as.factor(sample_data(merged_sub_ps_s45)$status_fixed)
sample_data(merged_sub_ps_s45)$Ucol <- as.factor(sample_data(merged_sub_ps_s45)$Ucol)
sample_data(merged_sub_ps_s45)$treatment <- as.factor(sample_data(merged_sub_ps_s45)$treatment)
sample_data(merged_sub_ps_s45)$agecat <- as.factor(sample_data(merged_sub_ps_s45)$agecat)
sample_data(merged_sub_ps_s45)$Sex <- as.factor(sample_data(merged_sub_ps_s45)$Sex)
sample_data(merged_sub_ps_s45)$treatment <- as.factor(sample_data(merged_sub_ps_s45)$treatment)



E1_merged_bray<- phyloseq::distance(merged_sub_ps_s45, method = "bray")
adonis2(E1_merged_bray ~  treatment+ agecat+Ucol, data = data.frame(sample_data(merged_sub_ps_s45)))
adonis2(E1_merged_bray ~ agecat + treatment +Ucol + Sex, data = data.frame(sample_data(merged_sub_ps_s45)))
adonis2(E1_merged_bray ~ treatment + agecat+ Sex +Ucol, data = data.frame(sample_data(merged_sub_ps_s45)))



###CAP###-----------------
ps.bray <- phyloseq::distance(merged_sub_ps_s45, method='bray')

ps.cap.ord <- ordinate(merged_sub_ps_s45, 
                       method = 'CAP',
                       distance = ps.bray,
                       formula = distance ~ agecat+treatment+Ucol ,
                       na.action = na.omit)

ps.cap.ord


#cap: Fig 3A
cap_plot <- plot_ordination(physeq = merged_sub_ps_s45, ordination = ps.cap.ord, color = 'Ucol')+
  geom_point(size=2)+
  stat_ellipse(type="norm",linetype=2)+
  theme_classic()

cap_plot


#Access CAP axis. 
cap.df <- cap_plot$data

cap.anova <- anova(ps.cap.ord, by='term')  #Ucol: p=0.017 
cap.anova.mar <- anova(ps.cap.ord, by='mar') #Ucol: p=0.024
cap.anova.axis <- anova(ps.cap.ord, by='axis')


##############################################################################
##MAASLIN2---------------------- 
##format before maaslin
merged_sub_ps_s45_df <- psmelt(merged_sub_ps_s45)
merged_sub_ps_s45_df<-dcast(merged_sub_ps_s45_df[c("Sample","Species","Abundance")],Sample~Species,fill=0.0, fun.aggregate = mean)
rownames(merged_sub_ps_s45_df) <- merged_sub_ps_s45_df$Sample
merged_sub_ps_s45_df$Sample <- NULL


## Figure 3B, 3C: differential taxa between Ucol and non-Ucol samples

fit_data = Maaslin2(
  input_data = merged_sub_ps_s45_df, 
  input_metadata = meta_s45, 
  transform = "AST",
  normalization = "NONE",
  standardize = FALSE,
  output = "./microbiome/temp", 
  fixed_effects = c('Ucol'), random_effects= c('agecat', 'treatment'))



##### Figure 3E: differential taxa between Ucol-nonrUTI and Ucol-rUTI

meta_s45_Ucol <- meta_s45[meta_s45$Ucol_kim=='1',]
merged_sub_ps_s45_df$Sample <- rownames(merged_sub_ps_s45_df)
temp_df <- merged_sub_ps_s45_df[merged_sub_ps_s45_df$Sample %in%meta_s45_Ucol$Patient_ID, ]
length(rownames(temp_df)) #33

temp_df$Sample <- NULL


fit_data = Maaslin2(
  input_data = temp_df, 
  input_metadata = meta_s45_Ucol, 
  transform = "AST",
  normalization = "NONE",
  standardize = FALSE,
  output = "./microbiome/temp2", 
  fixed_effects = c('status_fixed'), random_effects= c('treatment', 'agecat'))


######### calculate log2 fold change 

ecoli_coef <- 0.192454353
fc<- exp(ecoli_coef) # fold change (1.21)
log2(fc) #0.28


para_coef <- 0.131660693
fc<- exp(para_coef)  #(1.14)
log2(fc) #0.19

bacte_coef <- -0.263045394
fc<- exp(bacte_coef) # 0.77
log2(fc) # -0.4


##################################################
#####Figure 3D: Bubble plots for E. coli relative abundance
ecoli <- as.data.frame(species_e1.df$s__Escherichia_coli)
View(ecoli)


rownames(ecoli) <- rownames(species_e1.df)
colnames(ecoli) <- c('Abundance')

length(rownames(ecoli)) #480
length(rownames(meta_e1)) #480

#meta_e1$Abundance <- NULL


##merge E coli abundance data to metadata df 
meta_e1_ecoli <- merge(meta_e1, ecoli,
                 by = 'row.names', all = TRUE)

rownames(meta_e1_ecoli) <- meta_e1_ecoli$Row.names 


##checking if it worked 
rownames(ecoli)[256]
ecoli$Abundance[256]

meta_e1_ecoli["RH-10-E1-S5",]$Abundance

meta_Ucol_ecoli <- meta_e1_ecoli[meta_e1_ecoli$Ucol=='1',]
S5 <- meta_Ucol_ecoli[meta_Ucol_ecoli$timepoint=='5',]
S4 <- meta_Ucol_ecoli[meta_Ucol_ecoli$timepoint=='4',]
S3 <- meta_Ucol_ecoli[meta_Ucol_ecoli$timepoint=='3',]
S2 <- meta_Ucol_ecoli[meta_Ucol_ecoli$timepoint=='2',]


S4_filt <- S4[S4$Patient_ID %in% S3$Patient_ID & S4$Abundance >0.5,] #1 timepoint before 
S5_filt <-  S5[(S5$Patient_ID %in% S3$Patient_ID | S5$Patient_ID %in% S4$Patient_ID) & S5$Abundance >0.5,] # 1 or 2 timepoints before 


temp45 <- meta_Ucol_ecoli[meta_Ucol_ecoli$Patient_ID %in% S5_filt$Patient_ID | meta_Ucol_ecoli$Patient_ID %in% S4_filt$Patient_ID,]
length(unique(temp45$Patient_ID)) #6 patients 
ggplot(temp45, aes(x = as.factor(timepoint), y = as.factor(Patient_ID),size = Abundance ))+ geom_point(alpha = 0.7)+theme_bw()


ggplot(temp45, aes(x = as.factor(timepoint), y = as.factor(Patient_ID),size = Abundance ))+
  geom_point(alpha = 0.7)+theme_bw()+ scale_size(range = c(1,15))




