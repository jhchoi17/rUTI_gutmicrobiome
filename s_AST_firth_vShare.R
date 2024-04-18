
############################################
##script for interpreting AST data
##Usage: Figures 3G, 3H, 3F

'%not in%' <- function(x,y)!('%in%'(x,y))


ast <- read.csv('AST.csv')
 
length(colnames(ast)) #24 - 1 = 23 drugs tested 
length(rownames(ast)) #938 isolates tested 


##make metadata file for ast samples 

meta_ast <- data.frame(sample=ast$Isolate)

meta_ast$Patient_ID <- substr(meta_ast$sample,1,6)
meta_ast$episode <- substr(meta_ast$sample,8,9)
meta_ast$time <- substr(meta_ast$sample,11,13)
meta_ast$niche <- substr(meta_ast$sample,11,11)

meta_ast$Patient_ID <- gsub('-00', '-', meta_ast$Patient_ID)
meta_ast$Patient_ID <- gsub('-0', '-', meta_ast$Patient_ID)


#only keep those included in microbiome analysis 
meta <- read.csv('master_metadata.csv')
meta_ast_filt <- meta_ast[meta_ast$Patient_ID %in% meta$Patient_ID,]
length(rownames(meta_ast_filt)) #895 isolates 
length(unique(meta_ast_filt$Patient_ID)) #95 patients 

#add colonization data (from Thänert & Choi et al., Cell Host Microbe 2022)
library(data.table)

setDT(meta_ast_filt) #crucial to making below line work
setDT(meta)
meta_ast_filt <- meta_ast_filt[meta, c('Gcol_status', 'Ucol_status') := .(i.Gcol, i.Ucol), on = .(Patient_ID)]

##filter AST table .. 
ast_filt <- ast[ast$Isolate %in% meta_ast_filt$sample,]  


##### measuring AST score for each isolate profile 
#approach 1: quantify. 0=sensitive, 0.5=intermediate, 1=resistant 

for (i in 2:24){
  ast_filt[,i] <- gsub('Sensitive', 0, ast_filt[,i])
  ast_filt[,i] <- gsub('Resistant', 1, ast_filt[,i])
  ast_filt[,i] <- gsub('Intermediate', 0.5, ast_filt[,i])
}


guts <- meta_ast_filt[meta_ast_filt$niche=='S',]$sample
urine <- meta_ast_filt[meta_ast_filt$niche=='U' | meta_ast_filt$niche=='D',]$sample


length(guts) #551 Gut isolates
length(urine) #343

ast_gut <- ast_filt[ast_filt$Isolate %in% meta_gut_filt$sample,] ## only get gut isolates belonging to DxU strains!! , n=484
ast_urine <- ast_filt[ast_filt$Isolate %in% urine,] #343

#Note: #WU-004-E2-S06_E_coli fosfomycin is NA 
ast_gut[is.na(ast_gut)] <- 0   #correcting NA error 

rownames(ast_gut) <- ast_gut$Isolate
ast_gut$Isolate <- NULL

rownames(ast_urine) <- ast_urine$Isolate
ast_urine$Isolate <- NULL


#workaround to rowSums 
gutSums <- rowSums(sapply(ast_gut[, c(1:23)],
                          function(x) as.numeric(as.character(x))))

urSums<- rowSums(sapply(ast_urine[, c(1:23)],
                          function(x) as.numeric(as.character(x))))

#new df with iso names and gutSums/urSums 
meta_gut <- data.frame(sample=rownames(ast_gut), AST=gutSums)
meta_urine <- data.frame(sample=rownames(ast_urine), AST=urSums)

###################################comparison step #####
# are the scores diff between Ucol and non Ucol gut isolates?

## add in Ucol data, patient ID
setDT(meta_gut) #crucial to making below line work
setDT(meta_ast_filt)

##combine metadata to AST rowsums
meta_gut <- meta_gut[meta_ast_filt, c('Gcol_status', 'Ucol_status', 'Patient_ID','time') := .(i.Gcol_status, i.Ucol_status, i.Patient_ID, i.time), on = .(sample)]

setDT(meta_urine) #crucial to making below line work
setDT(meta_ast_filt)
meta_urine <- meta_urine[meta_ast_filt, c('Gcol_status', 'Ucol_status', 'Patient_ID','time') := .(i.Gcol_status, i.Ucol_status, i.Patient_ID, i.time), on = .(sample)]



###################################################
##annotate for lineage (strain) definitions 

straindf <- read.csv('210301_upec_isos_reference.csv')
#View(straindf)

length(unique(straindf$patient)) #93
length(unique(straindf$clusterID)) #95 


meta_gut$clusterID <- c('.')
meta_urine$clusterID <- c('.')

n <- length(rownames(meta_gut))

##annotate cluster ID to metadata

for (i in 1:n){
  iso <- as.character(meta_gut$sample[i])
  if (iso %in% straindf$sample){
    temp <- as.character(straindf[straindf$sample==iso,]$clusterID)
    meta_gut$clusterID[i] <- temp
  }
  
}

for (i in 1:n){
  iso <- as.character(meta_urine$sample[i])
  if (iso %in% straindf$sample){
    temp <- as.character(straindf[straindf$sample==iso,]$clusterID)
    meta_urine$clusterID[i] <- temp
  }
  
}


## filter out isolates with non-DxU lineages 
meta_gut_filt <- meta_gut[meta_gut$clusterID!='.',]
meta_urine_filt <- meta_urine[meta_urine$clusterID!='.',]

length(unique(meta_gut_filt$Patient_ID)) #79 Gut /  83 Urine patients 
length(unique(meta_urine_filt$clusterID)) #80 (WU-30 has 2 lineages! 1 for E1,E2 and another for E3), 85 Urine  (WU-30, PN-29)



##### now let's plot: isolates AST score. Ucol vs noncol

library(ggplot2)
library(ggsignif)
library(ggpubr)


##Figure 3G
ggplot(aes(x=Ucol_status, y=AST,color=Ucol_status), 
       data=meta_gut_filt)+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_point(position=position_jitterdodge(.5), size=2, alpha=0.5)+
  stat_compare_means(method = "kruskal.test")  + theme_bw() 

##Figure 3H
ggplot(aes(x=Ucol_status, y=AST,color=Ucol_status), 
       data=meta_urine_filt)+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_point(position=position_jitterdodge(.5), size=2, alpha=0.5)+
  stat_compare_means(method = "kruskal.test")  + theme_bw() 


##########################################
##testing significance of association for each drug type (Gut isolates)

length(rownames(ast_gut)) #484 gut isolates 
ast_gut$sample <- rownames(ast_gut)


Ucol <- meta_gut_filt[meta_gut_filt$Ucol_status=='1',]$sample ##already filtered for DxU strains only 
ast_gut_Ucol <- ast_gut[ast_gut$sample %in% Ucol,]

nonUcol <- meta_gut_filt[meta_gut_filt$Ucol_status=='.',]$sample
ast_gut_noncol <- ast_gut[ast_gut$sample %in% nonUcol,]

ast_gut_Ucol$sample <- NULL #234
ast_gut_noncol$sample <- NULL #250 

ast_gut_Ucol$Ucol <- c('1') #234
ast_gut_noncol$Ucol <- c('.') #250 


drugs <- colnames(ast_gut_Ucol)

ast_guts_annotated <- rbind(ast_gut_Ucol, ast_gut_noncol)
length(rownames(ast_guts_annotated)) #484



##################################################################################
################## Firth bias reduction method ######################


#install.packages("devtools")
#devtools::install_github("georgheinze/logistf")

library(logistf)


#format data into appropriate format
n <- length(rownames(ast_gut)) #484
ast_gut$Ucol <- ('temp') 


for (i in 1:n){
  temp <- as.character(ast_gut$sample[i])
  Ucolstat <- as.character(meta_ast_filt[meta_ast_filt$sample==temp, ]$Ucol_status)
  if(Ucolstat == '.'){
    temp2 <- 0
  }
  else{
    temp2 <- 1
  }
  ast_gut$Ucol[i] <- temp2
}


##replicate df and replace "intermediate" 0.5  with "non-resistant" 0
ast_gut2 <- ast_gut
ast_gut2$Isolate <- NULL


unique(factor(ast_gut2$Amikacin))

ast_gut2[ast_gut2==0.5] <- 0  
ast_gut2$Amikacin <- gsub("1 ", "1", ast_gut2$Amikacin) #typo 

unique(ast_gut2$Ucol)


#write.csv(ast_gut2, './Firth_input.csv')


# list of drugs to test:
#Ampicillin + Cefazolin + Cefotetan + Ceftriaxone + Ceftazidime + Cefepime
#             +Meropenem+Pipercillin.Tazobactam+Ceftolozane.Tazobactam+Ceftazidime.Avibactam+Ampicillin.Sulbactam
#              +Ciprofloxacin+Levofloxacin+Gentamicin+Amikacin+Trimethoprim.sulfa+Fosfomycin
#              +Aztreonam+Doxycycline+Minocycline+Tigecycline+Nitrofurantoin

### testing each drug independently 

ast_gut2$test <- ast_gut2$Imipenem

lf <- logistf(formula=as.integer(test) ~ as.integer(Ucol), ast_gut2) ## important to make sure these are integers
summary(lf)

#######################################################
##plot results : Figure 3F

library(ggplot2)

drugs <- read.csv('firth_results.csv')
drugs$drug <- as.factor(drugs$drug)

typeof(drugs$coef)


p1 <-ggplot(drugs, aes(x=coef, y=as.factor(pseudo)))+geom_point(size=3)+ylab('')+
  xlim(-10,10)+geom_segment(aes(x = lower.95., y = pseudo, xend = upper.95., yend = pseudo))+
  theme_bw()+  geom_vline(xintercept=0)+xlab('log(odds)') + scale_y_discrete(labels=drugs$drug)



################## import count data ######################

input <- read.csv('./AST_counts.csv')
n <- length(rownames(input)) #23
input <- as.data.frame((input))


### calculate % Ucol/ noncol and create bar graphs 

input$Uratio <- (input$Ucol_R/(input$Ucol_S+input$Ucol_R)) *100
input$nonUratio <- (input$noncol_R/(input$noncol_S+input$noncol_R))*100

input_filt <- input[input$drug %in% drugs$drug,]

# add 'pseudo' to df 
setDT(input_filt) #crucial to making below line work
setDT(drugs)
input_filt2 <- input_filt[drugs, c('pseudo') := .(i.pseudo), on = .(drug)]

input_filt$Ucol <- c(1) ## Ucol only 
input_filt$nonUratio <- NULL
input_filt2$Ucol <- c(0)
input_filt2$Uratio <- NULL

p2_data <- rbind(input_filt, input_filt2, use.names=FALSE)

p2_data$Ucol <- as.factor(p2_data$Ucol)

p2 <- ggplot(p2_data, aes(fill=Ucol, group=Ucol, y=as.factor(pseudo), x=Uratio))+geom_bar(position="dodge", stat="identity", width=.5)+ylab('')+
  theme_bw()

##combine p1 and p2: final Figure 3F
library(ggpubr)

ggarrange(p1,p2, widths=c(6,2))


