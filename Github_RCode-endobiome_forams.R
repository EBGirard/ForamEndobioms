#Statistical analyses and visualisation of the diatoms and prokaryote endobiomes datasets
#Author: Elsa B. Girard

options(rstudio.help.showDataPreview = FALSE)

rm(list=ls())
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
setSessionTimeLimit(cpu = Inf, elapsed = Inf)

#package needed
library(stringr)
library(pheatmap)
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(car)
library(patchwork)
library(ggvenn)
library(ggbreak)
library(scales)
library(cowplot)
library(ggpubr)
library(gghalves)
library(tidyverse)
library(rstatix)
library(decontam)
library(phyloseq)
library(viridis)  
library(cooccur)
library(visNetwork)
library(igraph)
library(RColorBrewer) 
library(heatmap3)
library(hrbrthemes)
library(indicspecies)
library(zetadiv)
library(stringi)
library(UpSetR)
library(VennDiagram)
library(ggpicrust2)
library(BiodiversityR) # also loads vegan
library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)

'%not%' <- Negate('%in%') #to create a 'not in' sign, the opposite of %in%
'is.not.na' <- Negate('is.na')
#to install packages: install.packages('package_name', lib='C:/Temp/R/R-4.1.2/library')
#if (!require('BiocManager', quietly = TRUE) + install.packages('BiocManager')
#BiocManager::install('tax4fun')

#________________________________

#--------DIATOMS: prep data---------------------
diat <- read.csv('~downloads/Datasets/DIATOM_apscale_ESV_table.csv', row.names = 1)
diat_tax <- read.csv('~downloads/Datasets/diatoms_all_taxo.csv', row.names = 1)
diat_met <- read.csv('~downloads/Datasets/Metadata_complete.csv', row.names = 1)

#create phyloseq object (make row names for metadata, taxonomy and esv table)
diat <- as.matrix(diat[,-(ncol(diat))])
diat_tax <- as.matrix(diat_tax)

OTU = otu_table(diat, taxa_are_rows = TRUE)
TAX = tax_table(diat_tax)
samples = sample_data(diat_met)

d_data <- phyloseq(OTU, TAX, samples)

#________________decontam package to filter from tag switching_____________________

d_neg <- read.csv('~downloads/Datasets/All_negatives.csv')

#sample_data(d_data)$Sample_or_Control <- NA
#sample_data(d_data)$Sample_or_Control[sample_data(d_data)$sample_cat == 'negative'] <- 'Control sample'
#sample_data(d_data)$Sample_or_Control[sample_data(d_data)$sample_cat != 'negative'] <- 'True sample'
sample_data(d_data)$total_reads <- NA
sample_data(d_data)$total_reads <- colSums(otu_table(d_data))

#NEGATIVES decontamination extractions first, then pcrs (column orders) in a loop

for (i in 1:ncol(d_neg)) {
  
  subset1 <- subset_samples(d_data, registr_code %in% d_neg[,i])
  sample_data(subset1)$is.neg <- sample_data(subset1)$Sample_or_Control == 'Control sample'
  contamdf.prev2 <- isContaminant(subset1, method='prevalence', neg='is.neg',threshold=0.2)
  table(contamdf.prev2$contaminant)
  contam_subset1 <- taxa_names(subset1)[which(contamdf.prev2$contaminant)]
  otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_subset1, colnames(otu_table(d_data)) %in% rownames(sample_data(subset1))] <- 0
  
}

diat <- data.frame(otu_table(d_data))
diat_met <- data.frame(sample_data(d_data))
diat_tax <- data.frame(tax_table(d_data))


#________________Filtering 0.1% + assignment________________________

#df[row,column]
df1 <- diat[rowSums(diat) > 2,]

#transformation switching rows and columns
df1 <- t(df1)
df2 <- as.data.frame(df1)

otudf <- df2
sum(otudf)

#transform raw values into percentages 
df2 <- t(apply(df2, 1, function(x) x/sum(x)))
#make sure row sum equal 1
rowSums(df2)

#remove asv < 0.1% of reads of the total reads but keeping reads number
otudf[df2 < 0.001] <- 0
#keep ESVs that still have reads 0.1%
otudf <- as.data.frame(otudf[,colSums(otudf) > 0]) #from 9760 esvs to 4296 esvs
sum(otudf)
#number of reads retained after 0.1% quality filter = 99.06%
230557715/232755051*100




#transform back data (switch rows and columns)
datafilt <- as.data.frame(t(otudf))

#create a column with ESVs
datafilt$ESVs <- row.names(datafilt) #filtered data
diat_tax$ESVs <- row.names(diat_tax)

#combine datasets
datafilt1 <- merge(datafilt,diat_tax, by = 'ESVs', all.x = TRUE)

colnames(datafilt1[,1000:ncol(datafilt1)])

#rearrange column order
datafilt2 <- datafilt1[, c(1,(ncol(datafilt1)-8):(ncol(datafilt1)-1),2:(ncol(datafilt1)-8))]

write.csv(datafilt2, '~downloads/Datasets/filtered_diat_esv_taxo.csv')

colnames(datafilt2)


#________________Add metadata________________________

#melt dataset to creat baseclear column
df_m <- melt(datafilt2[,c(1,2,5,10:ncol(datafilt2))], id.vars = c('ESVs', 'X.Subject', 'X.Identity.percentage'), 
             value.name = 'reads', variable.name = 'Baseclear_ID')

diat_met$Baseclear_ID <- as.character(row.names(diat_met))
df_m$Baseclear_ID <- as.character(as.factor(df_m$Baseclear_ID))

df2 <- merge(df_m, diat_met, by = 'Baseclear_ID', all = TRUE)

df2 <- df2[, -c(6,29)]


write.csv(df2, '~downloads/Datasets/filtered_diat_esv_taxo_meta.csv', row.names=FALSE)
df2<- read.csv('~downloads/Datasets/filtered_diat_esv_taxo_meta.csv')

#________________taxonomy________________________from R-syst::diatom website
#cite diatom r-syst:
#http:/dx.doi.org/10.5281/zenodo.31137
#http:/database.oxfordjournals.org/content/2016/baw016.full?keytype=ref&ijkey=H324uA95JzzEomz

df2 <- df2[df2$reads > 0,]

splitspecies <- str_split_fixed(df2$X.Subject, '/', 7)
colnames(splitspecies) <- c('Kindgom', 'Phylum', 'Class','Order', 'Family', 'Genus', 'Species')
df2 <- cbind(df2, splitspecies)

write.csv(df2, '~downloads/Datasets/filtered_diat_esv_taxo_meta_species.csv', row.names=FALSE)

#________________filter esvs that are in at least 2 replicates_________________

diat <- read.csv('~downloads/Datasets/filtered_diat_esv_taxo_meta_species.csv')

df1 <- diat
df1$pa[df1$reads > 0] <- 1
df1$pa[df1$reads == 0] <- 0

df1$sample_code <- substr(df1$registr_code, 1, nchar(df1$registr_code)-1)

df1$sample_code[df1$taxon == 'seawater'] <- NA

df2 <- df1 %>% group_by(ESVs, island, water_depth._m, sample_code, category, taxon, sequencin_run) %>% summarize(counts = sum(pa)) 

df3 <- df2[!(df2$counts  < 2), ]

df6 <- merge(df3, df1, by = c('ESVs', 'island', 'water_depth._m', 'sample_code', 'category', 'taxon', 'sequencin_run'))
df7 <- df6[,1:(ncol(df6)-1)]

new_df <- df7

new_df$dist_coast_km <- NA

#correct distances in km
new_df$dist_coast_km[new_df$island == 'Langkai'] <- 42
new_df$dist_coast_km[new_df$island == 'Samalona'] <- 7
new_df$dist_coast_km[new_df$island == 'Badi'] <- 23
new_df$dist_coast_km[new_df$island == 'Kodingareng Keke'] <- 15
new_df$dist_coast_km[new_df$island == 'Langkadea'] <- 14
new_df$dist_coast_km[new_df$island == 'Karanrang'] <- 14
new_df$dist_coast_km[new_df$island == 'Barang Baringan'] <- 5
new_df$dist_coast_km[new_df$island == 'Lae Lae'] <- 1
new_df$dist_coast_km[new_df$island == 'Lumu-lumu'] <- 31
new_df$dist_coast_km[new_df$island == 'Polewali'] <- 11
new_df$dist_coast_km[new_df$island == 'Pajenekang'] <- 19
new_df$dist_coast_km[new_df$island == 'Kapoposang'] <- 62

new_df <- new_df[new_df$usable == 'yes',]

write.csv(new_df, '~downloads/Datasets/filtered_diat_esv_taxo_meta_species_replifilt.csv')

#_________merge replicates___________________________

df <- read.csv('~downloads/Datasets/filtered_diat_esv_taxo_meta_species_replifilt.csv')

colnames(df)

df1 <- df %>% 
  group_by(category, taxon, island, sample_code, water_depth._m, dist_coast_km, long, lat, substrate_type, sequencin_run,
           ESVs, Kindgom,Phylum, Class, Order,Family,Genus,Species, X.Identity.percentage) %>%
  summarise(reads = sum(reads))


df2 <- df1[df1$Phylum == 'Bacillariophyta',] #assigned only to diatoms
df2$pa <- 1

write.csv(df2, '~downloads/Datasets/merged_DIATOM_data_all.csv', row.names = FALSE)

#check how much reads are kept after arranging and merging dataset: 94.95%
sum(df2$reads)

220992068/232755051*100 


#--------DIATOMS: prep working dataset------------------------
df <- read.csv('~downloads/Datasets/merged_DIATOM_data_all.csv')

df$sample_code[is.na(df$sample_code)] <- paste(df$taxon[is.na(df$sample_code)],df$island[is.na(df$sample_code)], df$water_depth._m[is.na(df$sample_code)],sep = '_')

#rarefy to the sample (none negative) with the lowest read number)
test <- df %>% 
  group_by(sample_code) %>% 
  summarise(N_sum=sum(reads))

dfx <- df %>% 
  group_by(sample_code) %>% 
  mutate(N=reads/sum(reads)*5927)

df2 <- dfx

df2a <- df2[df2$category != 'negative',]

df2a$site_code <- paste(df2a$island, df2a$water_depth._m, sep = '_')

#check ESVs kept, to make sure no important ones dropped. First 50 ESVs should be in.

splitESVs <- as.data.frame(str_split_fixed(df2a$ESVs, '_', 2))
colnames(splitESVs) <- c('ESV','ESV_code')
splitESVs$ESV_code <- as.numeric(splitESVs$ESV_code)
splitESVs <- unique(splitESVs)

write.csv(df2a,'~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv', row.names=FALSE)

#--------16S: prep data---------------------
prok <- read.csv('~downloads/Datasets/Proka_apscale_ESV_table.csv', row.names = 1)
prok_tax <- read.csv('~downloads/Datasets/prok_all_taxo.csv', row.names = 1)
prok_met <- read.csv('~downloads/Datasets/Metadata_complete.csv', row.names = 2)

#create phyloseq object (make row names for metadata, taxonomy and esv table)
prok <- as.matrix(prok[,-(ncol(prok))])
prok_tax <- as.matrix(prok_tax)

OTU = otu_table(prok, taxa_are_rows = TRUE)
TAX = tax_table(prok_tax)
samples = sample_data(prok_met)

d_data <- phyloseq(OTU, TAX, samples)

#________________decontam package to filter from tag switching_____________________

d_neg <- read.csv('~downloads/Datasets/All_negatives.csv')

#sample_data(d_data)$Sample_or_Control <- NA
#sample_data(d_data)$Sample_or_Control[sample_data(d_data)$sample_cat == 'negative'] <- 'Control sample'
#sample_data(d_data)$Sample_or_Control[sample_data(d_data)$sample_cat != 'negative'] <- 'True sample'
sample_data(d_data)$total_reads <- NA
sample_data(d_data)$total_reads <- colSums(otu_table(d_data))

#NEGATIVES decontamination extractions first, then pcrs (column orders) in a loop

for (i in 1:ncol(d_neg)) {
  
  subset1 <- subset_samples(d_data, registr_code %in% d_neg[,i])
  sample_data(subset1)$is.neg <- sample_data(subset1)$Sample_or_Control == 'Control sample'
  contamdf.prev2 <- isContaminant(subset1, method='prevalence', neg='is.neg',threshold=0.2)
  table(contamdf.prev2$contaminant)
  contam_subset1 <- taxa_names(subset1)[which(contamdf.prev2$contaminant)]
  otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_subset1, colnames(otu_table(d_data)) %in% rownames(sample_data(subset1))] <- 0
  
}

prok <- data.frame(otu_table(d_data))
prok_met <- data.frame(sample_data(d_data))
prok_tax <- data.frame(tax_table(d_data))

#________________Filtering 0.1% or 0.01% + assignment________________________

#df[row,column]
df1 <- prok[rowSums(prok) > 2,]

#transformation switching rows and columns
df1 <- t(df1)
df2 <- as.data.frame(df1)

otudf <- df2
sum(otudf)

#transform raw values into percentages 
df3 <- t(apply(df2, 1, function(x) x/sum(x)))
#make sure row sum equal 1
rowSums(df3)

#remove asv < 0.1% of reads of the total reads but keeping reads number
otudf[df3 < 0.001] <- 0
#keep ESVs that still have reads 0.1%
otudf <- as.data.frame(otudf[,colSums(otudf) > 0]) #from 61677 esvs to 14745 esvs
sum(otudf)
#number of reads retained after 0.1% quality filter = 92.43%
163577133/176978674*100


#transform back data (switch rows and columns)
datafilt <- as.data.frame(t(otudf))

#create a column with ESVs
datafilt$ESVs <- row.names(datafilt) #filtered data
prok_tax$ESVs <- row.names(prok_tax)

#combine datasets
datafilt1 <- merge(datafilt,prok_tax, by = 'ESVs', all.x = TRUE)

colnames(datafilt1[,1000:ncol(datafilt1)])

#rearrange column order
datafilt2 <- datafilt1[, c(1,(ncol(datafilt1)-7):(ncol(datafilt1)), (ncol(datafilt1)-14), (ncol(datafilt1)-12),2:(ncol(datafilt1)-16))]

write.csv(datafilt2, '~downloads/Datasets/filtered_prok_esv_taxo.csv', row.names = FALSE)


#________________Add metadata________________________

datafilt2 <- read.csv('~downloads/Datasets/filtered_prok_esv_taxo.csv')


#melt dataset to create baseclear column
df_m <- melt(datafilt2[,], id.vars = c('ESVs', 'X.Taxonomy', 'Kindgom', 'Phylum', 'Class','Order', 'Family', 'Genus', 'Species', 'X.Subject.accession', 'X.Identity.percentage'), 
             value.name = 'reads', variable.name = 'real_ID')

df_m$real_ID <- as.character(df_m$real_ID)
prok_met$real_ID <- as.character(row.names(prok_met))

df2 <- merge(df_m, prok_met, by = 'real_ID', all = TRUE)

df2 <- df2[df2$reads > 0,]

write.csv(df2, '~downloads/Datasets/filtered_prok_esv_taxo_meta.csv', row.names=FALSE)


#________________taxonomy________________________

df2 <- read.csv('~downloads/Datasets/filtered_prok_esv_taxo_meta.csv')

#splitspecies <- str_split_fixed(df2$X.Taxonomy, ' / ', 7)
#colnames(splitspecies) <- c('Kindgom', 'Phylum', 'Class','Order', 'Family', 'Genus', 'Species')
#df2 <- cbind(df2, splitspecies)

write.csv(df2, '~downloads/Datasets/filtered_prok_esv_taxo_meta_species.csv', row.names=FALSE)

#________________filter esvs that are in at least 2 replicates_________________

prok <- read.csv('~downloads/Datasets/filtered_prok_esv_taxo_meta_species.csv')

df1 <- prok
df1$pa[df1$reads > 0] <- 1
df1$pa[df1$reads == 0] <- 0

df1$sample_code <- substr(df1$registr_code, 1, nchar(df1$registr_code)-1)

df1$sample_code[df1$taxon == 'seawater'] <- NA

df2 <- df1 %>% group_by(ESVs, island, water_depth._m, sample_code, category, taxon, sequencin_run) %>% summarize(counts = sum(pa)) 

df3 <- df2[!(df2$counts  < 2), ]

df6 <- merge(df3, df1, by = c('ESVs', 'island', 'water_depth._m', 'sample_code', 'category', 'taxon', 'sequencin_run'))
df7 <- df6[,1:(ncol(df6)-1)]

new_df <- df7

new_df$dist_coast_km <- NA

#correct distances in km
new_df$dist_coast_km[new_df$island == 'Langkai'] <- 42
new_df$dist_coast_km[new_df$island == 'Samalona'] <- 7
new_df$dist_coast_km[new_df$island == 'Badi'] <- 23
new_df$dist_coast_km[new_df$island == 'Kodingareng Keke'] <- 15
new_df$dist_coast_km[new_df$island == 'Langkadea'] <- 14
new_df$dist_coast_km[new_df$island == 'Karanrang'] <- 14
new_df$dist_coast_km[new_df$island == 'Barang Baringan'] <- 5
new_df$dist_coast_km[new_df$island == 'Lae Lae'] <- 1
new_df$dist_coast_km[new_df$island == 'Lumu-lumu'] <- 31
new_df$dist_coast_km[new_df$island == 'Polewali'] <- 11
new_df$dist_coast_km[new_df$island == 'Pajenekang'] <- 19
new_df$dist_coast_km[new_df$island == 'Kapoposang'] <- 62

new_df <- new_df[new_df$usable == 'yes',]

write.csv(new_df, '~downloads/Datasets/filtered_prok_esv_taxo_meta_species_replifilt.csv', row.names = FALSE)

#_________remove chloroplast + mitochondria and organelle___________________

df <- read.csv('~downloads/Datasets/filtered_prok_esv_taxo_meta_species_replifilt.csv')

new_df <- df[df$Order != 'Chloroplast',]
new_df <- new_df[new_df$Family != 'Mitochondria',]
new_df <- new_df[new_df$Kindgom != 'Archaea',]
new_df <- new_df[!grepl('human', new_df$Species),] #human related contaminants
new_df <- new_df[!grepl('Cutibacterium acnes', new_df$Species),] #human related contaminants
new_df <- new_df[!grepl('Streptococcus', new_df$Species),] #human related contaminants
new_df <- new_df[!grepl('Streptococcus', new_df$Genus),] #human related contaminants
new_df <- new_df[!grepl('Burkholderia pseudomultivorans', new_df$Species),] #human related contaminants


sort(unique(new_df$Species))

write.csv(new_df, '~downloads/Datasets/filtered_prok_esv_taxo_meta_species_replifilt.csv', row.names = FALSE)

#_________merge replicates___________________________

df <- read.csv('~downloads/Datasets/filtered_prok_esv_taxo_meta_species_replifilt.csv')

colnames(df)

df1 <- df %>% 
  group_by(category, taxon, island, sample_code, water_depth._m, dist_coast_km, long, lat, substrate_type, sequencin_run,
           ESVs, Kindgom,Phylum, Class, Order,Family,Genus,Species,X.Subject.accession, X.Identity.percentage) %>%
  summarise(reads = sum(reads))

df1$pa <- 1


write.csv(df1, '~downloads/Datasets/merged_16S_data.csv', row.names = FALSE)



#--------16S: prep working dataset------------------------
df <- read.csv('~downloads/Datasets/merged_16S_data.csv')

df$sample_code[is.na(df$sample_code)] <- paste(df$taxon[is.na(df$sample_code)],df$island[is.na(df$sample_code)], df$water_depth._m[is.na(df$sample_code)],sep = '_')

#rarefy to the sample (none negative) with the lowest read number)
test <- df %>% 
  group_by(sample_code) %>% 
  summarise(N_sum=sum(reads))

list <- test$sample_code[test$N_sum < 1000]

df <- df[df$sample_code %not% list,]

dfx <- df %>% 
  group_by(sample_code) %>% 
  mutate(N=reads/sum(reads)*1130)

df2 <- dfx

df2a <- df2[df2$category != 'negative',]

df2a$site_code <- paste(df2a$island, df2a$water_depth._m, sep = '_')

#check ESVs kept, to make sure no important ones dropped. First 50 ESVs should be in.

splitESVs <- as.data.frame(str_split_fixed(df2a$ESVs, '_', 2))
colnames(splitESVs) <- c('ESV','ESV_code')
splitESVs$ESV_code <- as.numeric(splitESVs$ESV_code)
splitESVs <- unique(splitESVs)


write.csv(df2a,'~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv', row.names=FALSE)

#--------Figure 2 -- RDA analysis all----------------------
df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

dfx <- df[df$category == 'foram',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2b <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

df3 <- dcast(df2b, sample_code + depth_cat + island + dist_coast_km + taxon + substrate_type ~ ESVs, 
             fun.aggregate = sum, value.var = "percent_N")

df3[df3 == "NaN"] <- 0

df3$depth_cat <- as.factor(df3$depth_cat)
df3$island <- as.factor(df3$island)
df3$gaxon <- as.factor(df3$taxon)
df3$substrate_type <- as.factor(df3$substrate_type)

df_NMDS_vegan <- dcast(df2b, sample_code ~ ESVs, fun.aggregate = sum, value.var = "percent_N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$sample_code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:6]
attach(env)
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns


#RDA analysis
#https://rpubs.com/Roeland-KINDT/694016

Ordination.diatoms <- rda(cca_vegan ~ taxon , data=env, dist="bray", scaling="species")
summary(Ordination.diatoms)
anova(Ordination.diatoms) # overall test of the significant of the analysis
anova(Ordination.diatoms, by="axis", perm.max=500) # test axes for significance
anova(Ordination.diatoms, by="terms", permu=200) # test for sign. environ. variables

plot2 <- ordiplot(Ordination.diatoms, choices=c(1,2))
sites.long2 <- sites.long(plot2, env.data=env)
head(sites.long2)
species.long2 <- species.long(plot2)
species.long2
axis.long2 <- axis.long(Ordination.diatoms, choices=c(1, 2))
axis.long2

plot2 <- ordiplot(Ordination.diatoms, choices=c(1,2))
taxon.ellipses <- ordiellipse(plot2, groups=taxon, display="sites", kind="sd")
taxon.ellipses.long2 <- ordiellipse.long(taxon.ellipses, grouping.name="Taxon")

col_tax <- c('darkred', 'darkolivegreen3', 'skyblue1', 'turquoise','bisque3', 'darkgoldenrod1')

p1 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  
  geom_mark_ellipse(data=sites.long2,aes(x=axis1, y=axis2, colour=taxon, fill=after_scale(alpha(colour, 0.2))), expand=0, size=0.2, show.legend=FALSE) +
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=taxon, shape=taxon),  size=4) +
  scale_color_manual(values = col_tax)+
  theme_pubr()+
  coord_fixed(ratio=1)

Ordination.diatoms.all <- rda(cca_vegan ~ island + substrate_type + depth_cat, data=env, dist="bray", scaling="species")
summary(Ordination.diatoms.all)
anova(Ordination.diatoms.all) # overall test of the significant of the analysis
anova(Ordination.diatoms.all, by="axis", perm.max=500) # test axes for significance
anova(Ordination.diatoms.all, by="terms", permu=200) # test for sign. environ. variables



#_________________________________________________________________________________________________________________________
df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

dfx <- df[df$category == 'foram',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2b <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

df3 <- dcast(df2b, sample_code + depth_cat + island + dist_coast_km + taxon + substrate_type ~ ESVs, 
             fun.aggregate = sum, value.var = "percent_N")

df3[df3 == "NaN"] <- 0

df3$depth_cat <- as.factor(df3$depth_cat)
df3$island <- as.factor(df3$island)
df3$gaxon <- as.factor(df3$taxon)
df3$substrate_type <- as.factor(df3$substrate_type)

df_NMDS_vegan <- dcast(df2b, sample_code ~ ESVs, fun.aggregate = sum, value.var = "percent_N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$sample_code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:6]
attach(env)
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns


#RDA analysis
#https://rpubs.com/Roeland-KINDT/694016

Ordination.bacteria <- rda(cca_vegan ~ taxon , data=env, dist="bray", scaling="species")
summary(Ordination.bacteria)
anova(Ordination.bacteria) # overall test of the significant of the analysis
anova(Ordination.bacteria, by="axis", perm.max=500) # test axes for significance
anova(Ordination.bacteria, by="terms", permu=200) # test for sign. environ. variables


plot2 <- ordiplot(Ordination.bacteria, choices=c(1,2))
sites.long2 <- sites.long(plot2, env.data=env)
head(sites.long2)
species.long2 <- species.long(plot2)
species.long2
axis.long2 <- axis.long(Ordination.bacteria, choices=c(1, 2))
axis.long2

plot2 <- ordiplot(Ordination.bacteria, choices=c(1,2))
taxon.ellipses <- ordiellipse(plot2, groups=taxon, display="sites", kind="sd")
taxon.ellipses.long2 <- ordiellipse.long(taxon.ellipses, grouping.name="Taxon")

p2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  
  geom_mark_ellipse(data=sites.long2,aes(x=axis1, y=axis2, colour=taxon, fill=after_scale(alpha(colour, 0.2))), expand=0, size=0.2, show.legend=FALSE) +
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=taxon, shape=taxon),  size=4) +
  scale_color_manual(values = col_tax)+
  theme_pubr()+
  coord_fixed(ratio=1)

layout <- '
AB
'
#for ESV shared within a species
p1 + p2 + plot_layout(design = layout)


Ordination.bacteria.all <- rda(cca_vegan ~ island + substrate_type + depth_cat , data=env, dist="bray", scaling="species")
summary(Ordination.bacteria.all)
anova(Ordination.bacteria.all) # overall test of the significant of the analysis
anova(Ordination.bacteria.all, by="axis", perm.max=500) # test axes for significance
anova(Ordination.bacteria.all, by="terms", permu=200) # test for sign. environ. variables







#--------Figure 2 -- RDA analysis all per species----------------------
df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

#change host species name here: 
dfx <- df[df$taxon == 'Neorotalia calcar',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2b <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

df3 <- dcast(df2b, sample_code + depth_cat + island + dist_coast_km + taxon + substrate_type ~ ESVs, 
             fun.aggregate = sum, value.var = "percent_N")

df3[df3 == "NaN"] <- 0

df3$depth_cat <- as.factor(df3$depth_cat)
df3$island <- as.factor(df3$island)
df3$gaxon <- as.factor(df3$taxon)
df3$substrate_type <- as.factor(df3$substrate_type)

df_NMDS_vegan <- dcast(df2b, sample_code ~ ESVs, fun.aggregate = sum, value.var = "percent_N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$sample_code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:6]
attach(env)
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

Ordination.diatoms <- rda(cca_vegan ~ island, data=env, dist="bray", scaling="species")
summary(Ordination.diatoms)

plot2 <- ordiplot(Ordination.diatoms, choices=c(1,2))
sites.long2 <- sites.long(plot2, env.data=env)
head(sites.long2)
species.long2 <- species.long(plot2)
species.long2
axis.long2 <- axis.long(Ordination.diatoms, choices=c(1, 2))
axis.long2

plot2 <- ordiplot(Ordination.diatoms, choices=c(1,2))
taxon.ellipses <- ordiellipse(plot2, groups=island, display="sites", kind="sd")
taxon.ellipses.long2 <- ordiellipse.long(taxon.ellipses, grouping.name="Island")

col_tax <- c('darkred', 'darkolivegreen3', 'skyblue1', 'turquoise' ,'bisque3', 'darkgoldenrod1')

p1 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  
  geom_mark_ellipse(data=sites.long2,aes(x=axis1, y=axis2, colour=island, fill=after_scale(alpha(colour, 0.2))), expand=0, size=0.2, show.legend=FALSE) +
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=island, shape=island),  size=4) +
  scale_color_manual(values = col_tax)+
  theme_pubr()+
  coord_fixed(ratio=1)

#RDA analysis
#https://rpubs.com/Roeland-KINDT/694016

Ordination.diatoms.all <- rda(cca_vegan ~ island + substrate_type + depth_cat, data=env, dist="bray", scaling="species")
summary(Ordination.diatoms.all)
anova(Ordination.diatoms.all) # overall test of the significant of the analysis
anova(Ordination.diatoms.all, by="axis", perm.max=500) # test axes for significance
anova(Ordination.diatoms.all, by="terms", permu=200) # test for sign. environ. variables




#change host species name here: 
dfx <- df[df$taxon == 'Amphistegina lobifera',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2b <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

df3 <- dcast(df2b, sample_code + depth_cat + island + dist_coast_km + taxon + substrate_type ~ ESVs, 
             fun.aggregate = sum, value.var = "percent_N")

df3[df3 == "NaN"] <- 0

df3$depth_cat <- as.factor(df3$depth_cat)
df3$island <- as.factor(df3$island)
df3$gaxon <- as.factor(df3$taxon)
df3$substrate_type <- as.factor(df3$substrate_type)

df_NMDS_vegan <- dcast(df2b, sample_code ~ ESVs, fun.aggregate = sum, value.var = "percent_N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$sample_code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:6]
attach(env)
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

Ordination.diatoms <- rda(cca_vegan ~ island, data=env, dist="bray", scaling="species")
summary(Ordination.diatoms)

plot2 <- ordiplot(Ordination.diatoms, choices=c(1,2))
sites.long2 <- sites.long(plot2, env.data=env)
head(sites.long2)
species.long2 <- species.long(plot2)
species.long2
axis.long2 <- axis.long(Ordination.diatoms, choices=c(1, 2))
axis.long2

plot2 <- ordiplot(Ordination.diatoms, choices=c(1,2))
taxon.ellipses <- ordiellipse(plot2, groups=island, display="sites", kind="sd")
taxon.ellipses.long2 <- ordiellipse.long(taxon.ellipses, grouping.name="Island")

col_tax <- c('darkred', 'darkolivegreen3', 'skyblue1', 'turquoise' ,'bisque3', 'darkgoldenrod1')

p2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  
  geom_mark_ellipse(data=sites.long2,aes(x=axis1, y=axis2, colour=island, fill=after_scale(alpha(colour, 0.2))), expand=0, size=0.2, show.legend=FALSE) +
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=island, shape=island),  size=4) +
  scale_color_manual(values = col_tax)+
  theme_pubr()+
  coord_fixed(ratio=1)



#change host species name here: 
dfx <- df[df$taxon == 'Amphistegina radiata',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2b <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

df3 <- dcast(df2b, sample_code + depth_cat + island + dist_coast_km + taxon + substrate_type ~ ESVs, 
             fun.aggregate = sum, value.var = "percent_N")

df3[df3 == "NaN"] <- 0

df3$depth_cat <- as.factor(df3$depth_cat)
df3$island <- as.factor(df3$island)
df3$gaxon <- as.factor(df3$taxon)
df3$substrate_type <- as.factor(df3$substrate_type)

df_NMDS_vegan <- dcast(df2b, sample_code ~ ESVs, fun.aggregate = sum, value.var = "percent_N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$sample_code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:6]
attach(env)
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

Ordination.diatoms <- rda(cca_vegan ~ island, data=env, dist="bray", scaling="species")
summary(Ordination.diatoms)

plot2 <- ordiplot(Ordination.diatoms, choices=c(1,2))
sites.long2 <- sites.long(plot2, env.data=env)
head(sites.long2)
species.long2 <- species.long(plot2)
species.long2
axis.long2 <- axis.long(Ordination.diatoms, choices=c(1, 2))
axis.long2

plot2 <- ordiplot(Ordination.diatoms, choices=c(1,2))
taxon.ellipses <- ordiellipse(plot2, groups=island, display="sites", kind="sd")
taxon.ellipses.long2 <- ordiellipse.long(taxon.ellipses, grouping.name="Island")

col_tax <- c('darkred', 'darkolivegreen3', 'skyblue1', 'turquoise' ,'bisque3', 'darkgoldenrod1')

p3 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  
  geom_mark_ellipse(data=sites.long2,aes(x=axis1, y=axis2, colour=island, fill=after_scale(alpha(colour, 0.2))), expand=0, size=0.2, show.legend=FALSE) +
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=island, shape=island),  size=4) +
  scale_color_manual(values = col_tax)+
  theme_pubr()+
  coord_fixed(ratio=1)

layout <- '
ABC
'
#for ESV shared within a species
p1 + p2 + p3 + plot_layout(design = layout)

#_________________________________________________________________________________________________________________________
df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

#change host species name here: 
dfx <- df[df$taxon == 'Amphistegina lobifera',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2b <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

df3 <- dcast(df2b, sample_code + depth_cat + island + dist_coast_km + taxon + substrate_type ~ ESVs, 
             fun.aggregate = sum, value.var = "percent_N")

df3[df3 == "NaN"] <- 0

df3$depth_cat <- as.factor(df3$depth_cat)
df3$island <- as.factor(df3$island)
df3$gaxon <- as.factor(df3$taxon)
df3$substrate_type <- as.factor(df3$substrate_type)

df_NMDS_vegan <- dcast(df2b, sample_code ~ ESVs, fun.aggregate = sum, value.var = "percent_N")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$sample_code
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:6]
attach(env)
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns


#RDA analysis
#https://rpubs.com/Roeland-KINDT/694016

Ordination.bacteria.all <- rda(cca_vegan ~ island + substrate_type + depth_cat , data=env, dist="bray", scaling="species")
summary(Ordination.bacteria.all)
anova(Ordination.bacteria.all) # overall test of the significant of the analysis
anova(Ordination.bacteria.all, by="axis", perm.max=500) # test axes for significance
anova(Ordination.bacteria.all, by="terms", permu=200) # test for sign. environ. variables







#--------Dunn's test----------------

#_____________Diatoms_____________________________
df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$code <- paste(df$site_code,df$sample_code, sep = '_')

df_N <- dcast(df, code + taxon + category ~ ESVs, 
              fun.aggregate = sum, value.var = 'N')
df_N <- melt(df_N , id.vars = c('code', "category", 'taxon'), 
             value.name = 'N', variable.name = 'ESVs')
df_N$N[is.na(df_N$N)] <- 0

n_res <- dunn_test(df_N, N ~ taxon, p.adjust.method = "holm", detailed = FALSE)

df_pa <- dcast(df, code + taxon + category ~ ESVs, 
               fun.aggregate = sum, value.var = 'pa')
df_pa <- melt(df_pa , id.vars = c('code', "category", 'taxon'), 
              value.name = 'pa', variable.name = 'ESVs')
df_pa$pa[is.na(df_pa$pa)] <- 0

pa_res <- dunn_test(df_pa, pa ~ taxon, p.adjust.method = "holm", detailed = FALSE)

dunnRBCL <- rbind(pa_res, n_res)

#___________________Bacteria_____________________________

df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$code <- paste(df$site_code,df$sample_code, sep = '_')

df_N <- dcast(df, code + taxon + category ~ ESVs, 
              fun.aggregate = sum, value.var = 'N')
df_N <- melt(df_N , id.vars = c('code', "category", 'taxon'), 
             value.name = 'N', variable.name = 'ESVs')
df_N$N[is.na(df_N$N)] <- 0

n_res <- dunn_test(df_N, N ~ taxon, p.adjust.method = "holm", detailed = FALSE)

df_pa <- dcast(df, code + taxon + category ~ ESVs, 
               fun.aggregate = sum, value.var = 'pa')
df_pa <- melt(df_pa , id.vars = c('code', "category", 'taxon'), 
              value.name = 'pa', variable.name = 'ESVs')
df_pa$pa[is.na(df_pa$pa)] <- 0

pa_res <- dunn_test(df_pa, pa ~ taxon, p.adjust.method = "holm", detailed = FALSE)

dunn16S <- rbind(pa_res, n_res)



#--------Figure 3 -- bar chart all----------------------------------
df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

dfx <- df[df$category == 'foram',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df3 <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

#20 proportionally most abundant ESVs in the foram dataset
selected <- c('ESV_1','ESV_2','ESV_3','ESV_4','ESV_5','ESV_6','ESV_7',
              'ESV_8','ESV_9','ESV_10','ESV_11','ESV_12','ESV_13','ESV_14',
              'ESV_15','ESV_16','ESV_17','ESV_18','ESV_20','ESV_21', 'other ESVs')

df3$ESV_plot <- NA
df3$ESV_plot[df3$ESVs %in% selected] <- df3$ESVs[df3$ESVs %in% selected]
df3$ESV_plot[df3$ESVs %not% selected] <- 'other ESVs'
df3$spec_plot <- NA
df3$spec_plot[df3$ESVs %in% selected] <- df3$Species[df3$ESVs %in% selected]
df3$spec_plot[df3$ESVs %not% selected] <- 'other species'

df4 <- df3 %>% group_by(sample_code) %>% slice_max(percent_N, n = 1) %>% arrange(factor(ESV_plot, selected))
#df5 <- df3 %>% group_by(ESVs) %>% summarize(sum = sum(percent_N))
#df11 <- df3 %>% group_by(taxon, sample_code) %>% summarise(n = n())
#df12 <- df11 %>% group_by(taxon) %>% summarise(n = n())

p1 <- ggplot(df3, aes(x = factor(sample_code, df4$sample_code), y = percent_N, fill = factor(ESV_plot, selected), col = factor(ESV_plot, selected))) +
  geom_bar(stat='identity', position = 'fill', width=1) +
  facet_grid(~taxon, scale = 'free_x', space = "free_x")+
  scale_fill_manual(limits = selected, values = c(brewer.pal(10, "BrBG"), brewer.pal(10, "RdYlBu"), 'grey'))+
  scale_color_manual(limits = selected, values = c(brewer.pal(10, "BrBG"), brewer.pal(10, "RdYlBu"), 'grey'))+
  theme_pubr() + 
  theme(axis.text.x = element_blank(), legend.position = 'right', legend.title = element_blank())+
  xlab('Specimens') +
  ylab('Diatom RA (%) ESVs per sample')+
  labs(tag = 'A')

spec_ord1 <- c("Bacillaria_paxillifer" , "Navicula_sp.", "Nitzschia_aurariae", "Nitzschia_cf._microcephala", "Nitzschia_inconspicua",
               "Serratifera_corallina" ,  "Serratifera_sp.","Talaroneis_posidoniae" ,  "Thalassionema_frauenfeldii",  "other species" )

p2 <- ggplot(df3, aes(factor(sample_code, df4$sample_code), y = percent_N, fill = factor(spec_plot, spec_ord1), col = factor(spec_plot, spec_ord1))) +
  geom_bar(stat='identity', position = 'fill', width=1) +
  facet_grid(~taxon, scale = 'free_x', space = "free_x")+
  scale_fill_manual(limits = spec_ord1, values = c(brewer.pal(9, "BrBG"), 'grey'))+
  scale_color_manual(limits = spec_ord1, values = c(brewer.pal(9, "BrBG"), 'grey'))+
  theme_pubr() + 
  theme(axis.text.x = element_blank(), legend.position = 'right', legend.title = element_blank())+
  xlab('Specimens') +
  ylab('Diatom RA (%) species per sample') +
  labs(tag = 'B')


#_____________________________________________________________________
df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

dfx <- df[df$category == 'foram',]

df5 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df6 <- df5 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

#20 proportionally most abundant ESVs in the foram dataset
selected2 <- c('ESV_5','ESV_15','ESV_16',
               'ESV_17','ESV_18','ESV_19','ESV_20','ESV_21','ESV_25','ESV_28',
               'ESV_31','ESV_32','ESV_34','ESV_36','ESV_38','ESV_39','ESV_41',
               'ESV_43','ESV_46','ESV_151', 'other ESVs')

df6$ESV_plot <- NA
df6$ESV_plot[df6$ESVs %in% selected2] <- df6$ESVs[df6$ESVs %in% selected2]
df6$ESV_plot[df6$ESVs %not% selected2] <- 'other ESVs'
df6$spec_plot <- NA
df6$spec_plot[df6$ESVs %in% selected2] <- df6$Family[df6$ESVs %in% selected2]
df6$spec_plot[df6$ESVs %not% selected2] <- 'other families'

df7 <- df6 %>% group_by(sample_code) %>% slice_max(percent_N, n = 1) %>% arrange(factor(ESV_plot, selected2))
#df8 <- df6 %>% group_by(ESVs) %>% summarize(sum = sum(percent_N))
#df9 <- df6 %>% group_by(taxon, sample_code) %>% summarise(n = n())
#df10 <- df9 %>% group_by(taxon) %>% summarise(n = n())

p3 <- ggplot(df6, aes(x = factor(sample_code, df7$sample_code), y = percent_N, fill = factor(ESV_plot, selected2), col = factor(ESV_plot, selected2))) +
  geom_bar(stat='identity', position = 'fill', width=1) +
  facet_grid(~taxon, scale = 'free_x', space = "free_x")+
  scale_fill_manual(limits = selected2, values = c(brewer.pal(10, "BrBG"), brewer.pal(10, "RdYlBu"), 'grey'))+
  scale_color_manual(limits = selected2, values = c(brewer.pal(10, "BrBG"), brewer.pal(10, "RdYlBu"), 'grey'))+
  theme_pubr() + 
  theme(axis.text.x = element_blank(), legend.position = 'right', legend.title = element_blank())+
  xlab('Specimens') +
  ylab('Bacteria RA (%) ESVs per sample')+
  labs(tag = 'C')


spec_ord <- c("Alteromonadaceae","Burkholderiaceae", "Comamonadaceae",  "Enterobacteriaceae", "Moraxellaceae",
              "Pseudoalteromonadaceae", "Pseudomonadaceae",  "Rhizobiaceae", "Rhodobacteraceae",  "Stappiaceae"  , "Vibrionaceae",  "Xanthomonadaceae" , "unknownfamily", "other families")

p4 <- ggplot(df6, aes(factor(sample_code, df7$sample_code), y = percent_N, fill = factor(spec_plot, spec_ord), col = factor(spec_plot, spec_ord))) +
  geom_bar(stat='identity', position = 'fill', width=1) +
  facet_grid(~taxon, scale = 'free_x', space = "free_x")+
  scale_fill_manual(limits = spec_ord, values = c(brewer.pal(9, "BrBG"), brewer.pal(9, "RdYlBu")[c(1,2,3)], 'black','grey'))+
  scale_color_manual(limits = spec_ord, values = c(brewer.pal(9, "BrBG"), brewer.pal(9, "RdYlBu")[c(1,2,3)], 'black','grey'))+
  theme_pubr() + 
  theme(axis.text.x = element_blank(), legend.position = 'right', legend.title = element_blank())+
  xlab('Specimens') +
  ylab('Bacteria RA (%) families per sample') +
  labs(tag = 'D')

layout <- '
A
B
C
D
'
#for ESV shared within a species
p1 + p2 + p3 + p4 + plot_layout(design = layout)


#--------Figure 5 -- heatmaps all ----------------------------------
df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

#dfx <- df[df$category == 'foram',]

df2 <- df %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df3 <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

#20 proportionally most abundant ESVs in the foram dataset
selected <- c('ESV_1','ESV_2','ESV_3','ESV_4','ESV_5','ESV_6','ESV_7',
              'ESV_8','ESV_9','ESV_10','ESV_11','ESV_12','ESV_13','ESV_14',
              'ESV_15','ESV_16','ESV_17','ESV_18','ESV_20','ESV_21')

df3 <- df3[df3$ESVs %in% selected,]
df3$ESVs_spec <- paste(df3$Species, ':',df3$ESVs , sep= ' ')
df3$taxon[df3$taxon %in% c('seawater', 'substrate')] <- paste(df3$taxon[df3$taxon %in% c('seawater', 'substrate')],
                                                              df3$depth_cat[df3$taxon %in% c('seawater', 'substrate')], sep = ': ')

df4 <- df3 %>% group_by(category, taxon, ESVs_spec) %>% summarise(median_percent = round(median(percent_N),2), max_percent = max(percent_N), min_percent = min(percent_N))
df4 <- df4 %>% group_by(ESVs_spec) %>% mutate(total_max = max(median_percent))

df5 <- dcast(df4, taxon ~ ESVs_spec , value.var = 'median_percent')
row.names(df5) <- df5[,1]
df5 <- df5[,-1]
df5[is.na(df5)] <- 0

rowclus <- hclust(dist(df5)) #dendogram for sample types
colclus <- hclust(dist(t(df5))) #dendogram for bacteria ESVs

order_name <- data.frame(colnames(df5), 1:ncol(df5))   
order_name <- order_name %>% mutate(X1.ncol.df5.[order(colclus$order)])
ord <- rev(colnames(df5)[order(order_name$`X1.ncol.df5.[order(colclus$order)]`)])

order <- rev(c('seawater: flat','seawater: mid slope','seawater: deep slope','substrate: flat','substrate: mid slope','substrate: deep slope', 
               'Neorotalia calcar', 'Amphistegina lobifera', 'Amphistegina lessonii' ,'Amphistegina radiata','Calcarina spengleri', 'Heterostegina depressa' ))

p1 <- ggplot(df4, aes(y =factor(taxon, order), x=ESVs_spec)) +
  geom_tile(aes(fill = log(median_percent)), col = 'black') +
  facet_grid(.~factor(category, c('foram', 'substrate', 'seawater')), scales = 'free', space="free") +
  scale_fill_gradient2(low = 'floralwhite', mid = 'floralwhite', high = 'darkred') +
  scale_size(range = c(1, 10)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = 'bottom')+
  xlab('Twenty relatively most abundant ESVs and their taxonomic association') +
  ylab('Sample type') +
  ggtitle('Diatom endobiomes and microbiomes') +
  coord_flip() +
  labs(tag = 'A')




#_____________________________________________________________________
df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

#dfx <- df[df$category == 'foram',]

df6 <- df %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df7 <- df6 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

#20 proportionally most abundant ESVs in the foram dataset
selected2 <- c('ESV_5','ESV_15','ESV_16',
               'ESV_17','ESV_18','ESV_19','ESV_20','ESV_21','ESV_25','ESV_28',
               'ESV_31','ESV_32','ESV_34','ESV_36','ESV_38','ESV_39','ESV_41',
               'ESV_43','ESV_46','ESV_151')


df7 <- df7[df7$ESVs %in% selected2,]
df7$ESVs_spec <- paste(df7$Family, ':',df7$ESVs , sep= ' ')
df7$taxon[df7$taxon %in% c('seawater', 'substrate')] <- paste(df7$taxon[df7$taxon %in% c('seawater', 'substrate')],
                                                              df7$depth_cat[df7$taxon %in% c('seawater', 'substrate')], sep = ': ')

df8 <- df7 %>% group_by(category, taxon, ESVs_spec) %>% summarise(median_percent = round(median(percent_N),2), max_percent = max(percent_N), min_percent = min(percent_N))
df8 <- df8 %>% group_by(ESVs_spec) %>% mutate(total_max = max(median_percent))

order <- rev(c('seawater: flat','seawater: mid slope','seawater: deep slope','substrate: flat','substrate: mid slope','substrate: deep slope', 
               'Neorotalia calcar', 'Amphistegina lobifera', 'Amphistegina lessonii' ,'Amphistegina radiata','Calcarina spengleri', 'Heterostegina depressa' ))

p2 <- ggplot(df8, aes(y =factor(taxon, order), x=ESVs_spec)) +
  geom_tile(aes(fill = log(median_percent)), col = 'black') +
  facet_grid(.~factor(category, c('foram', 'substrate', 'seawater')), scales = 'free', space="free") +
  scale_fill_gradient2(low = 'floralwhite', mid = 'floralwhite', high = 'darkred') +
  scale_size(range = c(1, 10)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = 'bottom')+
  xlab('Twenty relatively most abundant ESVs and their taxonomic association') +
  ylab('Sample type') +
  ggtitle('Bacterial endobiomes and microbiomes') +
  coord_flip() +
  labs(tag = 'B')


layout <- '
AB
'
#for ESV shared within a species
p1 + p2 + plot_layout(design = layout)


#--------Figure 5 -- DIATOMS: venn diagram----------------------------

df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df_NMDS_vegan <- dcast(df, ESVs ~ taxon, fun.aggregate = sum, value.var = 'pa')
rownames(df_NMDS_vegan) <- df_NMDS_vegan$ESVs
df_NMDS_vegan <- df_NMDS_vegan[,-1]

#keep ESVs that are found it at least 2 specimens/samples
df2 <- df_NMDS_vegan
df2[df2 < 2] <- 0
df2[df2 > 0] <- 1
df2 <- df2[!rowSums(df2) == 0,]

df4 <- df[df$ESVs %in% unique(rownames(df2)),]

head(df2,5)

colnames(df2)

at_least_1 <- list(foram = df %>% filter(category=="foram") %>% select(ESVs) %>% unlist() , 
                   seawater = df %>% filter(category=="seawater") %>% select(ESVs) %>% unlist() , 
                   substrate = df %>% filter(category=="substrate") %>% select(ESVs) %>% unlist())

at_least_2 <- list(foram = df4 %>% filter(category=="foram") %>% select(ESVs) %>% unlist() , 
                   seawater = df4 %>% filter(category=="seawater") %>% select(ESVs) %>% unlist() , 
                   substrate = df4 %>% filter(category=="substrate") %>% select(ESVs) %>% unlist())

#to identify which ESVs are share between all three groups
names(which(colSums(table(stack(at_least_1)[2:1])> 0) > 2))

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(at_least_1, 
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("pink", 'blue3','bisque3'),
             # Numbers
             cex = 2,
             # Set names
             cat.cex = 2,
             cat.fontface = "bold",
             cat.col = c("pink", 'blue3','bisque3'),
             cat.default.pos = "outer")

display_venn(at_least_2, 
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("pink", 'blue3','bisque3'),
             # Numbers
             cex = 2,
             # Set names
             cat.cex = 2,
             cat.fontface = "bold",
             cat.col = c("pink", 'blue3','bisque3'),
             cat.default.pos = "outer")



#--------Figure 5 -- 16S: venn diagram----------------------------

#_____________foram vs substrate vs seawater___________________________

df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df_NMDS_vegan <- dcast(df, ESVs ~ category, fun.aggregate = sum, value.var = 'pa')
rownames(df_NMDS_vegan) <- df_NMDS_vegan$ESVs
df_NMDS_vegan <- df_NMDS_vegan[,-1]

#keep ESVs that are found it at least 2 specimens/samples
df2 <- df_NMDS_vegan
df2[df2 < 2] <- 0
df2[df2 > 0] <- 1
df2 <- df2[!rowSums(df2) == 0,]

df3 <- df[df$ESVs %in% unique(rownames(df2)),]

head(df2,5)

colnames(df2)


at_least_1 <- list(foram = df %>% filter(category=="foram") %>% select(ESVs) %>% unlist() , 
                   seawater = df %>% filter(category=="seawater") %>% select(ESVs) %>% unlist() , 
                   substrate = df %>% filter(category=="substrate") %>% select(ESVs) %>% unlist())

at_least_2 <- list(foram = df3 %>% filter(category=="foram") %>% select(ESVs) %>% unlist() , 
                   seawater = df3 %>% filter(category=="seawater") %>% select(ESVs) %>% unlist() , 
                   substrate = df3 %>% filter(category=="substrate") %>% select(ESVs) %>% unlist())

venn.diagram(at_least_1, filename = "C:/Temp/Chapter7-holobiont/Figure/venn-16s.png")

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(at_least_1, 
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("pink", 'blue3','bisque3'),
             # Numbers
             cex = 2,
             # Set names
             cat.cex = 2,
             cat.fontface = "bold",
             cat.col = c("pink", 'blue3','bisque3'),
             cat.default.pos = "outer")

display_venn(at_least_2, 
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("pink", 'blue3','bisque3'),
             # Numbers
             cex = 2,
             # Set names
             cat.cex = 2,
             cat.fontface = "bold",
             cat.col = c("pink", 'blue3','bisque3'),
             cat.default.pos = "outer")


#--------Figure S1 -- diatom: plot percent vs relative abundance-------------------------

df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

dfx <- df[df$category == 'foram',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)


#because the sequencing runs are comparable, we average the relative abundance for the samples that were in both sequencing run
df3 <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

#df3 <- df3[df3$ESVs %in% list_indic_ESVs,]

df3a <- df3 %>% 
  group_by(sample_code, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species,  ESVs, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(sum_percent=sum(percent_N))

df5 <- df3a %>% group_by(category, taxon, ESVs) %>% summarise(mean = median(sum_percent), sd = sd(sum_percent), counts = sum(pa))

df5$tot_nm_specimen <- NA
df5$tot_nm_specimen[df5$taxon == 'Amphistegina lobifera'] <- 35
df5$tot_nm_specimen[df5$taxon == 'Amphistegina lessonii'] <- 72
df5$tot_nm_specimen[df5$taxon == 'Amphistegina radiata'] <- 49
df5$tot_nm_specimen[df5$taxon == 'Calcarina spengleri'] <- 10
df5$tot_nm_specimen[df5$taxon == 'Heterostegina depressa'] <- 55
df5$tot_nm_specimen[df5$taxon == 'Neorotalia calcar'] <- 29

df5 <- df5 %>% group_by(taxon) %>% mutate(counts_percent_all = counts/tot_nm_specimen*100)
df5$upper <- df5$mean+df5$sd
df5$upper[df5$upper > 100] <- 100
df5$lower <- df5$mean-df5$sd
df5$lower[df5$lower < 0] <- 0
df5$color <- ifelse(df5$counts_percent_all > 75, df5$taxon, NA_character_)
df5$label <- ifelse(df5$counts_percent_all > 75, df5$ESVs, NA_character_)

p1 <- ggplot(df5, aes(x = counts_percent_all,y = mean, label = label, color = color)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 30) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9), col = 'black') +
  #facet_grid(taxon~., scales = 'free', space="free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom', legend.title = element_blank())+
  ylab('Mean relative abundance (%) +- sd') +
  xlab('Percentage presence in foraminifera specimens within the same host') +
  ggtitle('Diatoms') +
  coord_cartesian(ylim = c(0,100)) +
  labs(tag = 'A')

df6 <- df3a %>% group_by(category, ESVs) %>% summarise(mean = mean(sum_percent), sd = sd(sum_percent), counts = sum(pa))
df6 <- df6 %>% group_by(category) %>% mutate(counts_percent_all = counts/243*100)
df6$upper <- df6$mean+df6$sd
df6$upper[df6$upper > 100] <- 100
df6$lower <- df6$mean-df6$sd
df6$lower[df6$lower < 0] <- 0
df6$color <- ifelse(df6$counts_percent_all > 50 & df6$mean > 50, 'color', NA_character_)
df6$label <- ifelse(df6$counts_percent_all > 50, df6$ESVs, NA_character_)

p2 <- ggplot(df6, aes(x = counts_percent_all ,y = mean, label = label)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 30) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9), col = 'black') +
  #facet_grid(taxon~., scales = 'free', space="free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom', legend.title = element_blank())+
  ylab('Mean relative abundance (%) +- sd') +
  xlab('Percentage presence in foraminifera specimens overall') +
  ggtitle('Diatoms') +
  coord_cartesian(ylim = c(0,100)) +
  labs(tag = 'B')

layout <- '
AB
'
#for ESV shared within a species
p1 + p2 + plot_layout(design = layout)

#--------Figure S1 -- 16S: plot percent vs relative abundance-------------------------

df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')
dfx <- df[df$category == 'foram',]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

#_______________transform ESVs into species___________________________
df2$Species[df2$X.Identity.percentage < 94.5] <- 'unknown Bacteria'
df2$Genus[df2$X.Identity.percentage < 94.5] <- 'unknown Bacteria'
df2$Genus[df2$Genus == 'unknowngenus'] <- 'unknown Bacteria'
df2$Family[df2$X.Identity.percentage < 86.5] <- 'unknown Bacteria'
df2$Family[df2$Family == 'unknownfamily'] <- 'unknown Bacteria'
df2$Order[df2$Order == 'unknownorder'] <- 'unknown Bacteria'
df2$Class[df2$Class == 'unknownclass'] <- 'unknown Bacteria'


#because the sequencing runs are comparable, we average the relative abundance for the samples that were in both sequencing run
df3 <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

#df3 <- df3[df3$ESVs %in% list_indic_ESVs,]

df3a <- df3 %>% 
  group_by(sample_code, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species,  ESVs, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(sum_percent=sum(percent_N))



df5 <- df3a %>% group_by(category, taxon, ESVs) %>% summarise(mean = mean(sum_percent), sd = sd(sum_percent), counts = sum(pa))

df5$tot_nm_specimen <- NA
df5$tot_nm_specimen[df5$taxon == 'Amphistegina lobifera'] <- 35
df5$tot_nm_specimen[df5$taxon == 'Amphistegina lessonii'] <- 72
df5$tot_nm_specimen[df5$taxon == 'Amphistegina radiata'] <- 49
df5$tot_nm_specimen[df5$taxon == 'Calcarina spengleri'] <- 10
df5$tot_nm_specimen[df5$taxon == 'Heterostegina depressa'] <- 55
df5$tot_nm_specimen[df5$taxon == 'Neorotalia calcar'] <- 29

df5 <- df5 %>% group_by(taxon) %>% mutate(counts_percent_all = counts/tot_nm_specimen*100)
df5$upper <- df5$mean+df5$sd
df5$upper[df5$upper > 100] <- 100
df5$lower <- df5$mean-df5$sd
df5$lower[df5$lower < 0] <- 0
df5$color <- ifelse(df5$counts_percent_all > 75, df5$taxon, NA_character_)
df5$label <- ifelse(df5$counts_percent_all > 75, df5$ESVs, NA_character_)

p3 <- ggplot(df5, aes(x = counts_percent_all,y = mean, label = label, color = color)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 30) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9), col = 'black') +
  #facet_grid(taxon~., scales = 'free', space="free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom', legend.title = element_blank())+
  ylab('Mean relative abundance (%) +- sd') +
  xlab('Percentage presence in foraminifera specimens within the same host') +
  ggtitle('Bacteria') +
  coord_cartesian(ylim = c(0,100)) +
  labs(tag = 'C')

df6 <- df3a %>% group_by(category, ESVs) %>% summarise(mean = mean(sum_percent), sd = sd(sum_percent), counts = sum(pa))
df6 <- df6 %>% group_by(category) %>% mutate(counts_percent_all = counts/243*100)
df6$upper <- df6$mean+df6$sd
df6$upper[df6$upper > 100] <- 100
df6$lower <- df6$mean-df6$sd
df6$lower[df6$lower < 0] <- 0
df6$color <- ifelse(df6$counts_percent_all > 50 & df6$mean > 50, 'color', NA_character_)
df6$label <- ifelse(df6$counts_percent_all > 50, df6$ESVs, NA_character_)

p4 <- ggplot(df6, aes(x = counts_percent_all ,y = mean, label = label)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 30) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9), col = 'black') +
  #facet_grid(taxon~., scales = 'free', space="free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom', legend.title = element_blank())+
  ylab('Mean relative abundance (%) +- sd') +
  xlab('Percentage presence in foraminifera specimens overall') +
  ggtitle('Bacteria') +
  coord_cartesian(ylim = c(0,100)) +
  labs(tag = 'D')


layout <- '
AC
BD
'
#for ESV shared within a species
p3 + p4+ plot_layout(design = layout)




#--------Figure S2 -- diversity index all-----------------------

#_______________DIATOMS______________________


df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

dfx <- df[df$taxon %not% c('seawater', 'substrate'),]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2$X.Identity.percentage <- df2$X.Identity.percentage/1000

#_______________transform ESVs into species___________________________
df2$species_95 <- NA
df2$species_95[df2$X.Identity.percentage >= 95] <- df2$Species[df2$X.Identity.percentage >= 95]
df2$species_95[df2$X.Identity.percentage < 95] <- paste(df2$Species[df2$X.Identity.percentage < 95], '< 95% ID', sep = '_')

#because the sequencing runs are comparable, we average the relative abundance for the samples that were in both sequencing run
df3 <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, species_95, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

#df3 <- df3[df3$ESVs %in% list_indic_ESVs,]
df3$code <- paste(df3$island, df3$depth_cat ,df3$sample_code, sep = '_')

df_diat <- dcast(df3, code ~ ESVs, fun.aggregate = sum, value.var = 'percent_N')
rownames(df_diat) <- df_diat$code
df_diat[is.na(df_diat)] <- 0


#_______________PROKARYOTES______________________


df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

dfx <- df[df$taxon %not% c('seawater', 'substrate'),]

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

#_______________transform ESVs into species___________________________
df2$Species[df2$X.Identity.percentage < 94.5] <- 'unknown Bacteria'
df2$Genus[df2$X.Identity.percentage < 94.5] <- 'unknown Bacteria'
df2$Genus[df2$Genus == 'unknowngenus'] <- 'unknown Bacteria'
df2$Family[df2$X.Identity.percentage < 86.5] <- 'unknown Bacteria'
df2$Family[df2$Family == 'unknownfamily'] <- 'unknown Bacteria'
df2$Order[df2$Order == 'unknownorder'] <- 'unknown Bacteria'
df2$Class[df2$Class == 'unknownclass'] <- 'unknown Bacteria'


#because the sequencing runs are comparable, we average the relative abundance for the samples that were in both sequencing run
df3 <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, X.Subject.accession,X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

#df3 <- df3[df3$ESVs %in% list_indic_ESVs,]
df3$code <- paste(df3$island, df3$depth_cat ,df3$sample_code, sep = '_')

df_bac <- dcast(df3, code ~ ESVs, fun.aggregate = sum, value.var = 'percent_N')
rownames(df_bac) <- df_bac$code
df_bac[is.na(df_bac)] <- 0


#_____________________diversity indexes______________________

#shannon
bac_shannon <- vegan::diversity(df_bac[-1], index="shannon")
diat_shannon <- vegan::diversity(df_diat[-1], index="shannon")

#dominance
bac_dominance <- vegan::diversity(df_bac[-1], index="simpson")
diat_dominance <- vegan::diversity(df_diat[-1], index="simpson")

#evenness
bac_S <- apply(df_bac[,-1] > 0, 1, sum)
bac_evenness <- vegan::diversity(df_bac[-1], index="simpson")/log(bac_S)

diat_S <- apply(df_diat[,-1] > 0, 1, sum)
diat_evenness <- vegan::diversity(df_diat[-1], index="simpson")/log(diat_S)

l <- tibble::lst(diat_evenness, bac_evenness, diat_dominance, bac_dominance, diat_shannon, bac_shannon, diat_S, bac_S)
res_df <- data.frame(lapply(l, `length<-`, max(lengths(l))))
res_df$code <- row.names(res_df)

res_df <- melt(res_df, id.vars = c('code'), 
               value.name = 'diversity_value', variable.name = 'diversity_index')

res_df$endobiont <- NA
res_df$endobiont[res_df$diversity_index %in% c('diat_evenness','diat_dominance','diat_shannon', 'diat_S')] <- 'diatoms'
res_df$endobiont[res_df$diversity_index %in% c('bac_evenness','bac_dominance','bac_shannon', 'bac_S')] <- 'bacteria'

res_df$index <- NA
res_df$index[res_df$diversity_index %in% c('diat_evenness','bac_evenness')] <- 'Pilou evenness'
res_df$index[res_df$diversity_index %in% c('bac_dominance','diat_dominance')] <- 'Simpson Index (1-D)'
res_df$index[res_df$diversity_index %in% c('diat_shannon','bac_shannon')] <- 'Shannon-Wiener Index'
res_df$index[res_df$diversity_index %in% c('diat_S','bac_S')] <- 'Species richness'

ggplot(res_df, aes(x = endobiont, y = diversity_value, fill = endobiont)) +
  geom_boxplot(alpha = 0.5) +
  geom_signif(comparisons = list(c('bacteria', 'diatoms')), test = 't.test',
              map_signif_level=TRUE, tip_length = 0.01, vjust = .1, col = 'black') +
  facet_wrap(~index, scales = 'free', ncol = 4) +
  scale_fill_manual(values = c('darkred','floralwhite')) +
  theme_classic() +
  ggtitle('Considering all ESVs')



#--------Figure S3 -- DIATOMS: - occupancy model---------------------

library(camtrapR)
library(purrr)
library(DT)
library(knitr)
library(sf)
library(rjags)
#JAGS_ROOT="C:/Users/elsa.girard/AppData/Roaming/Microsoft/Windows/Start Menu/Programs/JAGS"


#_____PA - Transform dataframes into a list (https:/cran.r-project.org/web/packages/cit/vignettes/camtrapr5.html)____________

df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

#df$Genus[df$Species == 'Serratifera_sp.'] <- 'Serratifera'

dfx <- df

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2$X.Identity.percentage <- df2$X.Identity.percentage/1000

#_______________transform ESVs into species___________________________
df2$species_95 <- NA
df2$species_95[df2$X.Identity.percentage >= 95] <- df2$Species[df2$X.Identity.percentage >= 95]
df2$species_95[df2$X.Identity.percentage < 95] <- paste(df2$Species[df2$X.Identity.percentage < 95], '< 95% ID', sep = '_')

#because the sequencing runs are comparable, we average the relative abundance for the samples that were in both sequencing run
df3 <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, species_95, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

df3a <- df3 %>% 
  group_by(sample_code, species_95, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(sum_percent=sum(percent_N))


#find list of top 2 ESVs per sample
df4_2 <- df3a %>% arrange(desc(sum_percent)) %>% group_by(island, depth_cat, sample_code) %>% slice(1:2)
#find list of top 1 ESVs per sample
df4_1 <- df3a %>% arrange(desc(sum_percent)) %>% group_by(island, depth_cat, sample_code) %>% slice(1:1)
df4_1$dominant_diatom <- NA

list_samples <- unique(df4_1$sample_code)

for (i in 1:length(list_samples)) {
  
  subset_df <- df4_2[df4_2$sample_code == list_samples[[i]],]
  
  if (df4_2$sum_percent[1]/df4_2$sum_percent[2] >= 2) {
    
    df4_1$dominant_diatom[df4_1$sample_code == list_samples[[i]]] <- 'yes'
    
  } else{
    
    df4_1$dominant_diatom[df4_1$sample_code == list_samples[[i]]] <- 'no'
    
  }
  
}


df5 <- df4_1 %>% group_by(taxon, species_95) %>% summarise(mean = mean(sum_percent), sd = sd(sum_percent), counts = sum(pa))
df5 <- df5 %>% group_by(taxon) %>% mutate(counts_percent = round(counts/sum(counts)*100,0))

list_dominant <- unique(df5$species_95[df5$taxon %not% c('seawater', 'substrate') & df5$counts > 1])
df6 <- df3a[df3a$species_95 %in% list_dominant,]
df6$sample_code[df6$taxon %in% c('seawater', 'substrate')] <- paste(df6$taxon[df6$taxon %in% c('seawater', 'substrate')],df6$island[df6$taxon %in% c('seawater', 'substrate')],
                                                                    df6$depth_cat[df6$taxon %in% c('seawater', 'substrate')], sep = ': ')


df_model <- dcast(df6, species_95 ~ sample_code, 
                  fun.aggregate = sum, value.var = 'sum_percent')
rownames(df_model) <- df_model[,1]
df_model <- df_model[,-1]


#from here, we need a nmds style matrix

OTU_Index_pa <- vegan::decostand(df_model, method = "pa")
OTU_Index1 <- as.matrix(t(OTU_Index_pa))

# Initialize an empty list to store rows
ylist <- list()

# Iterate over each row of the matrix
for (i in 1:ncol(OTU_Index1)) {
  # Extract each row and add it to the list
  ylist[[i]] <- OTU_Index1[, i, drop = FALSE]
}

names(ylist) <- colnames(OTU_Index1)

metadata_reordered <- data.frame(unique(df6[order(df6$sample_code),c('sample_code', 'island', 'taxon','category', 'depth_cat', 'substrate_type', 'dist_coast_km')]))
metadata_reordered$taxon[metadata_reordered$taxon %in% c('substrate','seawater')] <- 'environment'
metadata_reordered$sample_code <- as.factor(metadata_reordered$sample_code)
metadata_reordered$taxon <-  as.factor(as.character(metadata_reordered$taxon))
metadata_reordered$taxon <- factor(metadata_reordered$taxon, 
                                   levels = c( 'environment','Neorotalia calcar','Amphistegina lobifera','Amphistegina lessonii','Amphistegina radiata' ,
                                               'Calcarina spengleri','Heterostegina depressa'))
metadata_reordered$depth_cat <-  as.factor(metadata_reordered$depth_cat)
metadata_reordered$category <- as.factor(metadata_reordered$category)
metadata_reordered$island <- as.factor(metadata_reordered$island)
metadata_reordered$dist_coast_km <- as.numeric(metadata_reordered$dist_coast_km)
metadata_reordered$substrate_type <- as.factor(metadata_reordered$substrate_type)


unique(metadata_reordered$taxon)

rownames(metadata_reordered) <- NULL

effort <- matrix(1, nrow =  nrow(OTU_Index1))

input_data <- list(ylist    = ylist,
                   siteCovs = metadata_reordered,
                   obsCovs  = list(effort = effort))  # is identical for all species

modelFile_jags_categ1 <- tempfile(fileext = ".txt")

model_jags_categ1 <- communityModel(
  occuCovs = list(ranef = c('taxon')),
  detCovsObservation = list(ranef = "effort"),
  intercepts = list(det = "ranef", occu = "ranef"),
  data_list = input_data,
  modelFile = modelFile_jags_categ1)

# short model run for demonstration
out_ahm_jags_categ1 <- fit(model_jags_categ1, 
                           n.iter = 1000, 
                           n.burnin = 250,
                           thin = 5,
                           chains = 3,
                           quiet = T
)




#histogram/bar chart
plot_effects(object = model_jags_categ1,
             mcmc.list = out_ahm_jags_categ1)

#per category all species
plot <- plot_coef(object = model_jags_categ1,
                  mcmc.list = out_ahm_jags_categ1,
                  ordered = FALSE,
                  level = c(outer = 0.999, inner = 0.95))

results <- as.data.frame(plot$taxon$data)
results1 <- results[results$species != "community",]
results1$significance[results1$significance == "no"] <- "p > 0.05"
results1$significance[results1$significance == "outer"] <- "p < 0.001"
results1$significance[results1$significance == "inner"] <- "p < 0.05"

taxon_lab <- c("Environment", "A. lessonii",  "A. lobifera",  "A. radiata",  "C. spengleri",  "H. depressa",  "N. calcar")

results1$covariate <- factor(results1$covariate, levels = c("environment", "Amphistegina lessonii",  "Amphistegina lobifera",  
                                                            "Amphistegina radiata",  "Calcarina spengleri",  "Heterostegina depressa",  
                                                            "Neorotalia calcar"),labels = taxon_lab)

ggplot(results1, aes(y = species, shape = significance, col = significance)) +
  geom_vline(aes(xintercept = 0)) +
  geom_vline(aes(xintercept = -10), col = "gray") +
  geom_vline(aes(xintercept = 10), col = "gray") +
  geom_linerange(aes(xmin=lower_inner, xmax=upper_inner), linewidth = 1) +
  geom_linerange(aes(xmin=lower_outer, xmax=upper_outer), linewidth = 0.5) +
  geom_point(aes(x = median), size = 3) +
  #geom_point(aes(x = median), size = 3, col = "black") +
  scale_shape_manual(values = c(17, 15, 19))+
  scale_color_manual(values = c("pink", "black", "gray"))+
  facet_grid(.~covariate) + 
  theme_classic() +
  xlab("Effect size")+
  ylab("Dominant diatom per host foraminifera")



summary(model_jags_categ1)
info <- summary(out_ahm_jags_categ1)
#--------Figure S4 -- DIATOMS: community composition most abundant diatom species (barchart and heatmap)---------------
df <- read.csv('~downloads/Datasets/workingdataset_DIATOM_all_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'
df$category[df$taxon == 'seawater'] <- 'seawater'
df$category[df$taxon == 'substrate'] <- 'substrate'

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')

dfx <- df

df2 <- dfx %>% 
  group_by(sample_code, island, depth_cat, sequencin_run) %>% 
  mutate(percent_N=N/sum(N)*100)

df2$X.Identity.percentage <- df2$X.Identity.percentage/1000

#_______________transform ESVs into species___________________________
df2$species_95 <- NA
df2$species_95[df2$X.Identity.percentage >= 95] <- df2$Species[df2$X.Identity.percentage >= 95]
df2$species_95[df2$X.Identity.percentage < 95] <- paste(df2$Species[df2$X.Identity.percentage < 95], '< 95% ID', sep = '_')

#because the sequencing runs are comparable, we average the relative abundance for the samples that were in both sequencing run
df3 <- df2 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, species_95, X.Identity.percentage, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(percent_N=mean(percent_N))

df3a <- df3 %>% 
  group_by(sample_code, ESVs, category, taxon, island, depth_cat, Kindgom, Class, Order, Family, Genus, Species, pa, 
           substrate_type, dist_coast_km, water_depth._m) %>% 
  summarize(sum_percent=sum(percent_N))

df3a$ord_gen <- paste(substr(df3a$Order, 1,4),df3a$Genus, sep = ': ')

df4_1 <- df3a %>% arrange(desc(sum_percent)) %>% group_by(island, depth_cat, sample_code) %>% slice(1:1)

df4_1$ord_spe <- paste(substr(df4_1$Order, 1,4),df4_1$Species, df4_1$ESVs, sep = ': ')

df5 <- df4_1 %>% group_by(taxon, Order, ord_spe) %>% summarise(mean = mean(sum_percent), sd = sd(sum_percent), counts = sum(pa))


col3 <- c('darkred','red2','blue4', 'bisque3','bisque1', 'turquoise') #for diatom in > 1 species

ggplot(df5[df5$taxon %not% c('seawater', 'substrate'),], aes(x = ord_spe ,y = mean, fill = factor(Order, sort(unique(Order)))), labels = paste(counts, '%', sep = '')) +
  geom_col(position = 'dodge', col = 'black', alpha = 0.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean), width=.2, position=position_dodge(.9), col = 'black') +
  facet_grid(~taxon, scales = 'free', space="free") + 
  geom_text(aes(label = counts), y = 2) +
  scale_fill_manual(values = col3)+
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom', legend.title = element_blank())+
  ylab('Relative abundance (%)') +
  xlab('Most abundant ESV per specimen')



#_______________Check that the two sequencing runs are comparable___________________________
#find list of top 1 ESVs per sample
df4 <- df2 %>% arrange(desc(percent_N)) %>% group_by(island, depth_cat, sample_code, sequencin_run) %>% slice(1:10)
df4$ord_spe <- paste(substr(df4$Order, 1,4),df4$Species, df4$ESVs, round(df4$X.Identity.percentage/1000,2), sep = ': ')
df4$ord_gen <- paste(substr(df4$Order, 1,4),df4$Genus, sep = ': ')

df_1 <- df4[df4$sequencin_run == 1, c('sample_code', 'sequencin_run', 'percent_N', 'ESVs')]
df_2 <- df4[df4$sequencin_run == 2,c('sample_code', 'sequencin_run', 'percent_N', 'ESVs')]

df_12 <- merge(df_1, df_2, c('sample_code'))

df_test <- df4[df4$sample_code %in% unique(df_12$sample_code),]

df_test$sequencin_run <- as.factor(df_test$sequencin_run)

#check similarity between sequencing runs
ggplot(df_test, aes(x = sequencin_run, y = percent_N, fill = ESVs)) +
  geom_bar(stat='identity', position = 'stack') +
  facet_grid(~sample_code, scale = 'free_x')+
  xlab('Sequencing run') +
  ylab('Relative abundance (%) of the top 10 ESVs per sample, per run')

#__________________________________________________________________________________________





#--------Figure S5 -- 16S: number of shared esvs 0.1%--------------------------------------------------

df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df$sample_code[df$taxon == 'seawater'] <- 'seawater'
df$sample_code[df$taxon == 'substrate'] <- 'substrate'


#resulting dataframe
site_code <- c(1:2000)
taxon <- NA
V1 <- NA
V2 <- NA
nmbr_shared_esvs <- NA
nmbr_all_esvs <- NA
percent_shared_esvs <- NA
results <- data.frame(site_code, taxon, V1, V2, nmbr_shared_esvs, nmbr_all_esvs, percent_shared_esvs)

microhab <- unique(df$site_code[!is.na(df$site_code)])

r <- 0

#make for loop here
for (i in 1:length(microhab)) { #going through different sites
  
  df1 <- df[df$site_code == microhab[[i]],]
  
  foram <- unique(df1$taxon[df1$taxon %not% c('seawater', 'substrate')])
  
  for (x in 1:length(foram)) { #going through different species in that site
    
    sp <- unique(df1$sample_code[df1$taxon == foram[[x]] & !is.na(df1$taxon)])
    sp_subs_water <- c(sp, 'substrate', 'seawater')
    
    combi <- as.data.frame(t(combn(sp_subs_water,2))) #create all combinations
    
    for (y in 1:nrow(combi)) { #doing all combinations for that species in that site
      
      r <- r + 1
      
      dfv1 <- subset(df1, df1$sample_code == combi[y,1])
      dfv2 <- subset(df1, df1$sample_code == combi[y,2])
      
      dfcom <- unique(merge(dfv1[,c('site_code','ESVs')], dfv2[,c('site_code','ESVs')], by = c('site_code','ESVs')))
      dfall <- unique(merge(dfv1[,c('site_code','ESVs')], dfv2[,c('site_code','ESVs')], by = c('site_code','ESVs'), all = TRUE))
      
      results$site_code[r] <- microhab[[i]]
      results$taxon[r] <- foram[[x]]
      results$V1[r] <- combi[y,1]
      results$V2[r] <- combi[y,2]
      results$nmbr_shared_esvs[r] <- nrow(dfcom)
      results$nmbr_all_esvs[r] <- nrow(dfall)
      results$percent_shared_esvs[r] <- nrow(dfcom)/nrow(dfall) * 100
    } 
  }
}

results <- results[!is.na(results$taxon),]


#resulting dataframe shared between same species across sites
site_code <- c(1:10000)
taxon <- NA
V1 <- NA
V2 <- NA
nmbr_shared_esvs <- NA
nmbr_all_esvs <- NA
percent_shared_esvs <- NA
results2 <- data.frame(site_code, taxon, V1, V2, nmbr_shared_esvs, nmbr_all_esvs, percent_shared_esvs)

r <- 0

foram <- unique(df$taxon[df$taxon %not% c('seawater', 'substrate')])

for (x in 1:length(foram)) { #going through different species overall
  
  sp <- unique(df$sample_code[df$taxon == foram[[x]] & !is.na(df$taxon)])
  
  combi <- as.data.frame(t(combn(sp,2))) #create all combinations
  
  for (y in 1:nrow(combi)) { #doing all combinations for that species overall
    
    r <- r + 1
    
    dfv1 <- subset(df, df$sample_code == combi[y,1])
    dfv2 <- subset(df, df$sample_code == combi[y,2])
    
    dfcom <- unique(merge(dfv1[,c('taxon','ESVs')], dfv2[,c('taxon','ESVs')], by = c('taxon','ESVs')))
    dfall <- unique(merge(dfv1[,c('taxon','ESVs')], dfv2[,c('taxon','ESVs')], by = c('taxon','ESVs'), all = TRUE))
    
    results2$taxon[r] <- foram[[x]]
    results2$V1[r] <- combi[y,1]
    results2$V2[r] <- combi[y,2]
    results2$nmbr_shared_esvs[r] <- nrow(dfcom)
    results2$nmbr_all_esvs[r] <- nrow(dfall)
    results2$percent_shared_esvs[r] <- nrow(dfcom)/nrow(dfall) * 100
  } 
}


results2 <- results2[!is.na(results2$taxon),]

results3 <- rbind(results[,c(2:7)], results2[,c(2:7)])

results3 <- unique(results3[!is.na(results3$V1),])

results3$category <- NA
results3$category[results3$V2 == 'seawater'] <- 'shared foram - seawater same microhabitat'
results3$category[results3$V2 == 'substrate'] <- 'shared foram - substrate same microhabitat'
results3$category[results3$V1 %in% c('seawater', 'substrate') & results3$V2 %in% c('seawater', 'substrate')] <- 'shared seawater - substrate same microhabitat'

for (i in 1:nrow(results3)) { 
  
  if (results3$V1[i] %not% c('seawater', 'substrate') & results3$V2[i] %not% c('seawater', 'substrate')) {
    
    if (unique(df$site_code[df$sample_code == results3$V1[i]]) == unique(df$site_code[df$sample_code == results3$V2[i]])) {
      
      results3$category[i] <- 'shared foram sp1 - foram sp1 same microhabitat'
      
    } else if (unique(df$site_code[df$sample_code == results3$V1[i]]) != unique(df$site_code[df$sample_code == results3$V2[i]])) {
      
      results3$category[i] <- 'shared foram sp1 - foram sp1 different microhabitat'
      
    } else {
      
      results3$category[i] <- results3$category[i]
    }
  }
}


unique(results3$category)

#resulting dataframe shared between different species from same site
site_code <- c(1:10000)
taxon <- NA
V1 <- NA
V2 <- NA
nmbr_shared_esvs <- NA
nmbr_all_esvs <- NA
percent_shared_esvs <- NA
results4 <- data.frame(site_code, taxon, V1, V2, nmbr_shared_esvs, nmbr_all_esvs, percent_shared_esvs)

microhab <- unique(df$site_code[!is.na(df$site_code)])

r <- 0

for (i in 1:length(microhab)) { #going through different sites
  
  df1 <- df[df$site_code == microhab[[i]],]
  
  foram <- unique(df1$sample_code[df1$taxon %not% c('seawater', 'substrate')])
  
  combi <- as.data.frame(t(combn(foram,2))) #create all combinations
  
  for (y in 1:nrow(combi)) { #doing all combinations for that species in that site
    
    r <- r + 1
    
    dfv1 <- subset(df1, df1$sample_code == combi[y,1])
    dfv2 <- subset(df1, df1$sample_code == combi[y,2])
    
    dfcom <- unique(merge(dfv1[,c('site_code','ESVs')], dfv2[,c('site_code','ESVs')], by = c('site_code','ESVs')))
    dfall <- unique(merge(dfv1[,c('site_code','ESVs')], dfv2[,c('site_code','ESVs')], by = c('site_code','ESVs'), all = TRUE))
    
    results4$site_code[r] <- microhab[[i]]
    results4$V1[r] <- combi[y,1]
    results4$V2[r] <- combi[y,2]
    results4$nmbr_shared_esvs[r] <- nrow(dfcom)
    results4$nmbr_all_esvs[r] <- nrow(dfall)
    results4$percent_shared_esvs[r] <- nrow(dfcom)/nrow(dfall) * 100
  } 
}

results4 <- results4[!is.na(results4$V1),]

results4$category <- NA

for (i in 1:nrow(results4)) { 
  
  if (unique(df$taxon[df$sample_code == results4$V1[i]]) == unique(df$taxon[df$sample_code == results4$V2[i]])) {
    
    results4$category[i] <- 'shared foram sp1 - foram sp1 same microhabitat'
    
  } else if (unique(df$taxon[df$sample_code == results4$V1[i]]) != unique(df$taxon[df$sample_code == results4$V2[i]])) {
    
    results4$category[i] <- 'shared foram sp1 - foram sp2 same microhabitat'
    
  }
}

col1 <- c('black','black','black','black', 'black','black' ,'blue3','bisque3')
df5 <- df %>% group_by(site_code, taxon, sample_code) %>% summarise(sum_pa = sum(pa)) 

df5 <- df5 %>% group_by(taxon) %>% mutate(mean_pa = mean(sum_pa))

p1 <- ggplot(df5, aes(x = taxon, y = sum_pa, col = taxon, fill = taxon)) +
  geom_jitter() +
  geom_boxplot(alpha = 0.3) +
  #stat_summary(aes(y = mean_pa), linewidth= 0.3, geom = 'crossbar') + 
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  theme_pubr() +
  theme(axis.title.x = element_blank(), legend.position = 'none') +
  ylab('Number of ESVs') +
  ggtitle('Number of bacteria ESVs found in foraminifera, substrate and seawater samples')


col2 <- c('pink','darkred','bisque3', 'blue3','turquoise')
results3 <- results3 %>% group_by(taxon, category) %>% mutate(mean_percent = mean(percent_shared_esvs),mean_nmbr = round(mean(nmbr_shared_esvs)), mean_all = round(mean(nmbr_all_esvs)))

p2 <- ggplot(results3, aes(x = category, y = percent_shared_esvs, col = category, fill = category)) +
  geom_violin(alpha = 0.3) +
  facet_grid(~taxon) +
  stat_summary(aes(y = mean_percent), linewidth= 0.3, geom = 'crossbar') + 
  geom_text(aes(label = mean_nmbr), y = -10) +
  geom_text(aes(label = mean_all), y = -20) +
  scale_color_manual(values = col2) +
  scale_fill_manual(values = col2) +
  geom_signif(comparisons = list(c('shared foram sp1 - foram sp1 different microhabitat', 'shared foram sp1 - foram sp1 same microhabitat')), test = 't.test',
              map_signif_level=TRUE,y_position = 90, tip_length = 0.01, vjust = .1, col = 'black') + 
  geom_signif(comparisons = list(c('shared foram - substrate same microhabitat', 'shared foram sp1 - foram sp1 same microhabitat')), test = 't.test',
              map_signif_level=TRUE,y_position = 80, tip_length = 0.01, vjust = .1, col = 'black') + 
  geom_signif(comparisons = list(c('shared foram - seawater same microhabitat', 'shared foram sp1 - foram sp1 same microhabitat')), test = 't.test',
              map_signif_level=TRUE,y_position = 70, tip_length = 0.01, vjust = .1, col = 'black') + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, -10, -20),labels = c('0%', '25%', '50%', '75%', '100%', 'mean # ESVs shared', 'mean # ESVs total')) +
  coord_cartesian( ylim = c(-20, 100), clip = 'off') +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  ylab('Proportion of shared ESVs') +
  ggtitle('Propotion of bacteria ESVs shared between 2 specimens of the same species')

order_group <- c("shared foram sp1 - foram sp1 same microhabitat","shared foram sp1 - foram sp1 different microhabitat","shared foram - substrate same microhabitat",
                 "shared foram - seawater same microhabitat","shared seawater - substrate same microhabitat")

p2a <- ggplot(results3, aes(x = factor(category, order_group), y = nmbr_shared_esvs, col = factor(category, order_group), fill = factor(category, order_group))) +
  geom_violin(alpha = 0.3) +
  facet_grid(~taxon) +
  stat_summary(aes(y = mean_nmbr), linewidth= 0.3, geom = 'crossbar') + 
  geom_text(aes(label = mean_all), y = -5) +
  scale_color_manual(values = col2) +
  scale_fill_manual(values = col2) +
  geom_signif(comparisons = list(c('shared foram sp1 - foram sp1 different microhabitat', 'shared foram sp1 - foram sp1 same microhabitat')), test = 't.test',
              map_signif_level=TRUE,y_position = 35, tip_length = 0.01, vjust = .1, col = 'black') + 
  geom_signif(comparisons = list(c('shared foram - substrate same microhabitat', 'shared foram sp1 - foram sp1 same microhabitat')), test = 't.test',
              map_signif_level=TRUE,y_position = 30, tip_length = 0.01, vjust = .1, col = 'black') + 
  geom_signif(comparisons = list(c('shared foram - seawater same microhabitat', 'shared foram sp1 - foram sp1 same microhabitat')), test = 't.test',
              map_signif_level=TRUE,y_position = 25, tip_length = 0.01, vjust = .1, col = 'black') + 
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, -5),labels = c('0','10', '20', '30', '40' , 'mean # ESVs total')) +
  coord_cartesian( ylim = c(-5, 40), clip = 'off') +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  ylab('Proportion of shared ESVs') +
  xlab('Category of shared ESVs') +
  ggtitle('Propotion of bacteria ESVs shared between 2 specimens of the same species')

results3a <- results3[results3$category %in% c("shared foram sp1 - foram sp1 same microhabitat","shared foram sp1 - foram sp1 different microhabitat"),]
results3b <- results3a %>% group_by(taxon, category) %>% summarize(mean_shared = mean(nmbr_shared_esvs), sd_shared = sd(nmbr_shared_esvs))

col4 <- c('pink','darkred')

p2b <- ggplot(results3b, aes(x = factor(category, order_group), y = mean_shared, col = factor(category, order_group), fill = factor(category, order_group))) +
  geom_jitter(data = results3a, aes(x = factor(category, order_group), y = nmbr_shared_esvs), shape = 1) +
  geom_bar(stat = 'identity', alpha = 0.3, col = 'black') +
  facet_grid(~taxon) +
  geom_errorbar(aes(ymin=mean_shared-sd_shared, ymax=mean_shared+sd_shared), width=.2, position=position_dodge(.9), col = 'black')+
  geom_signif(data = results3a, aes(x = category, y = nmbr_shared_esvs), comparisons = list(c("shared foram sp1 - foram sp1 same microhabitat","shared foram sp1 - foram sp1 different microhabitat")), test = 't.test',
              map_signif_level=TRUE, y_position = 22, tip_length = 0.01, vjust = .1, col = 'black') +
  scale_color_manual(values = col4) +
  scale_fill_manual(values = col4) +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  ylab('Proportion of shared ESVs') +
  ggtitle('Propotion of bacteria ESVs shared between 2 specimens of the same species')


col3 <- c('pink','red2')
results5 <- results4 %>% group_by(site_code, category) %>% mutate(mean_percent = mean(percent_shared_esvs),mean_nmbr = round(mean(nmbr_shared_esvs)), mean_all = round(mean(nmbr_all_esvs)))

splitloc <- as.data.frame(str_split_fixed(results5$site_code, '_', 2))
colnames(splitloc) <- c('island', 'depth_m')
splitloc$depth_m <- as.numeric(splitloc$depth_m)
results5  <- cbind(results5,splitloc)
results5$depth_cat <- NA
results5$depth_cat[results5$depth_m == 1] <- 'flat'
results5$depth_cat[results5$depth_m > 15] <- 'deep slope'
results5$depth_cat[is.na(results5$depth_cat)] <- 'mid slope'
depth <- c('flat', 'mid slope', 'deep slope')
island_ord <- c('Samalona', 'Kodingareng Keke', 'Pajenekang', 'Badi', 'Langkai', 'Kapoposang')

p3 <- ggplot(results5, aes(x = category, y = percent_shared_esvs, col = category, fill = category)) +
  geom_violin(alpha = 0.3) +
  facet_grid(~site_code) +
  stat_summary(aes(y = mean_percent), linewidth= 0.3, geom = 'crossbar') + 
  geom_text(aes(label = mean_nmbr), y = -10) +
  geom_text(aes(label = mean_all), y = -20) +
  geom_signif(comparisons = list(c('shared foram sp1 - foram sp2 same microhabitat', 'shared foram sp1 - foram sp1 same microhabitat')), test = 't.test',
              map_signif_level=TRUE, y_position = 90, tip_length = 0.01, vjust = .1, col = 'black') +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, -10),labels = c('0%', '25%', '50%', '75%', '100%', 'mean # ESVs shared','mean # ESVs total')) +
  coord_cartesian( ylim = c(-20, 100), clip = 'off') +
  scale_color_manual(values = col3) +
  scale_fill_manual(values = col3) +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  ylab('Proportion of shared ESVs') +
  ggtitle('Proportion of bacteria ESVs shared between 2 specimens of the same microhabitat')

results6 <- results5 %>% group_by(island, depth_cat, category) %>% summarize(mean_shared = mean(nmbr_shared_esvs), sd_shared = sd(nmbr_shared_esvs))

p3a <- ggplot(results6, aes(x = category, y = mean_shared, col = category, fill = category)) +
  geom_jitter(data = results5, aes(x = category, y = nmbr_shared_esvs), shape = 1) +
  geom_bar(stat = 'identity', alpha = 0.3, col = 'black') +
  facet_grid(factor(depth_cat, depth)~factor(island, island_ord)) +
  geom_errorbar(aes(ymin=mean_shared-sd_shared, ymax=mean_shared+sd_shared), width=.2, position=position_dodge(.9), col = 'black')+
  geom_signif(data = results5, aes(x = category, y = nmbr_shared_esvs), comparisons = list(c('shared foram sp1 - foram sp2 same microhabitat', 'shared foram sp1 - foram sp1 same microhabitat')), test = 't.test',
              map_signif_level=TRUE, y_position = 22, tip_length = 0.01, vjust = .1, col = 'black') +
  scale_color_manual(values = col3) +
  scale_fill_manual(values = col3) +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  ylab('Proportion of shared ESVs') +
  ggtitle('Proportion of bateria ESVs shared between 2 specimens of the same microhabitat')

#finish
#--------Figure S5 -- 16S: ploting shared esv result analysis 0.1%----------------------
layout <- '
A
B
B
B
'
#for ESV shared within a species
p2b + p3a + plot_layout(design = layout)

layout <- '
A
B
C
'

#for ESV shared within a species
p1 + p2a + p3a + plot_layout(design = layout)


layout <- '
A
B
C
'

#for ESV shared within a species
p1 + p2 + p4 + plot_layout(design = layout)

#for ESV shared within a microhabitat
p3 + p5 + plot_layout(design = layout)



#--------Figure S6 -- 16S: METAGENassist-------------------------------
df <- read.csv('~downloads/Datasets/workingdataset_16S_withoutneg_01percent.csv')

df <- df[df$taxon %not% c('substrate', 'seawater'),]

df$depth_cat <- NA
df$depth_cat[df$water_depth._m == 1] <- 'flat'
df$depth_cat[df$water_depth._m > 15] <- 'deep slope'
df$depth_cat[is.na(df$depth_cat)] <- 'mid slope'

df$code <- paste(df$site_code,df$sample_code, sep = '_')
df$code2 <- paste(df$taxon,df$depth_cat, df$island, sep = '_')



#arrange taxonomy based on Yarza et al. 2014

df$Species[df$X.Identity.percentage < 94.5] <- 'unknownspecies'
df$Genus[df$X.Identity.percentage < 94.5] <- 'unknowngenus'
df$Family[df$X.Identity.percentage < 86.5] <- 'unknownfamily'
df$Order[df$X.Identity.percentage < 82] <- 'unknownorder'

df$taxonomy_all <- paste(df$Kindgom, df$Phylum, df$Class, df$Order, df$Family, df$Genus, df$Species, sep = ';')

#__________________all samples__________________________________________

df$SampleID <- df$code
df$value <- round(df$N)

df_ESV <- dcast(df, SampleID ~ taxonomy_all, 
                fun.aggregate = sum, value.var = 'value')

write.csv(df_ESV, '~downloads/Datasets/METAGENassist_taxonomic_profile_16S_all.csv', row.names = FALSE)

df_metadata <- unique(df[, c('SampleID','island','depth_cat','substrate_type','taxon','category')])
order_sample <- df_ESV$SampleID
df_metadata <- df_metadata %>% arrange(factor(SampleID, levels = order_sample))

write.csv(df_metadata, '~downloads/Datasets/METAGENassist_metadata_16S_all.csv', row.names = FALSE)

#__________________Taxon_________________________________________

df$SampleID <- df$code2
df$value <- round(df$N)

df_ESV <- dcast(df, SampleID ~ taxonomy_all, 
                fun.aggregate = sum, value.var = 'value')

write.csv(df_ESV, '~downloads/Datasets/METAGENassist_taxonomic_profile_16S_taxon_all.csv', row.names = FALSE)

df_metadata <- unique(df[, c('SampleID','island','depth_cat','substrate_type','taxon')])
order_sample <- df_ESV$SampleID
df_metadata <- df_metadata %>% arrange(factor(SampleID, levels = order_sample))

write.csv(df_metadata, '~downloads/Datasets/METAGENassist_metadata_16S_taxon_all.csv', row.names = FALSE)

#____________________Plot metabolism heatmap______________________________
res <- read.csv('~downloads/Datasets/METAGENassist_metabolism_results_all.csv')

#res <- res[,-7] #remove unknown category

res_melt <- melt(res, id.vars = c('sample_id'), 
                 value.name = 'abundance', variable.name = 'metabolism')

colnames(res_melt) <- c('SampleID', 'InferedMetabolism', 'Abundance')

res_final <- merge(df_metadata, res_melt, by = 'SampleID', all = TRUE)

res_final <- res_final %>% 
  group_by(SampleID) %>% 
  mutate(percent=Abundance/sum(Abundance)*100)

res15 <- res_final %>% arrange(desc(percent)) %>% group_by(taxon) %>% slice(1:15)

res_final15 <- res_final[res_final$InferedMetabolism %in% unique(res15$InferedMetabolism),]



ggplot(res_final15, aes(y = InferedMetabolism, x=paste(res_final15$depth_cat, res_final15$island, sep = ' : '))) +
  geom_tile(aes(fill = percent), color = 'black') +
  #geom_point(aes(col = log(max_percent)), size = 7) +
  #geom_point(aes(col = log(min_percent)), size = 3) +
  #geom_text(aes(label = mean_percent))+
  facet_grid(~taxon, scales = 'free', space = 'free') +
  scale_fill_gradient2(low = 'blue4', mid = 'floralwhite', high = 'darkred') +
  scale_size(range = c(1, 10)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = 'bottom')+
  xlab('Infered taxonomy-based functional roles of bacteria') +
  ylab('Samples') +
  ggtitle('Infered taxonomy-based functional roles of bacteria in different foram host species')

#--------Figure S7 -- DIATOMS: define identity threshold from similarity matrix based on rbcl database---------------------------

dfm <- read.csv("~downloads/Datasets/diatom_rcbl_rsyst_v8-idmatrix.csv", row.names = 1)


dfm$taxonomy1 <- rownames(dfm)

df1 <- melt(dfm, id.vars = c('taxonomy1'), 
            value.name = 'percent_ID', variable.name = 'taxonomy2' )

df1$percent_ID[is.na(df1$percent_ID)] <- 100

split <- str_split_fixed(df1$taxonomy1, "_", 7)
colnames(split) <- c('kingdom1', 'phylum1','class1','order1','family1','genus1','species1')
df2 <- cbind(df1, split)

split <- str_split_fixed(df1$taxonomy2, "_", 7)
colnames(split) <- c('kingdom2', 'phylum2','class2','order2','family2','genus2','species2')
df2 <- cbind(df2, split)

split <- str_split_fixed(df2$species1, " ", 2)
colnames(split) <- c('species1a', 'other')
df2 <- cbind(df2, split[,1])

split <- str_split_fixed(df2$species2, " ", 2)
colnames(split) <- c('species2a', 'other')
df2 <- cbind(df2, split[,1])

colnames(df2) <- c("taxonomy1" , "taxonomy2" , "percent_ID" ,"kingdom1"  , "phylum1" ,   "class1"  ,   "order1"    , "family1"  ,  "genus1"    , "species1" ,  "kingdom2" ,  "phylum2" ,   "class2",    
                   "order2"  ,   "family2"  ,  "genus2"   ,  "species2"  , "species"  ,  "genera"   ,  "families",   "orders"  ,   "classes"  ,  "species1a", "species2a")

df2$species <- NA
df2$species[df2$species1a == df2$species2a] <- "same species"
df2$species[df2$species1a != df2$species2a] <- "different species"
df2$genera <- NA
df2$genera[df2$genus1 == df2$genus2] <- "same genus"
df2$genera[df2$genus1 != df2$genus2] <- "different genera"
df2$families <- NA
df2$families[df2$family1 == df2$family2] <- "same family"
df2$families[df2$family1 != df2$family2] <- "different families"
df2$orders <- NA
df2$orders[df2$order1 == df2$order2] <- "same order"
df2$orders[df2$order1 != df2$order2] <- "different orders"
df2$classes <- NA
df2$classes[df2$class1 == df2$class2] <- "same order"
df2$classes[df2$class1 != df2$class2] <- "different orders"

df2$percent_ID[df2$percent_ID == 0] <- NA

df3 <- df2 %>% group_by(genera) %>% summarize(mean_ID = mean(percent_ID, na.rm = TRUE),
                                              min_ID = min(percent_ID, na.rm = TRUE), 
                                              max_ID = max(percent_ID, na.rm = TRUE),
                                              sd_ID = sd(percent_ID, na.rm = TRUE))
colnames(df3) <- c('category', 'mean_ID', 'min_ID', 'max_ID', 'sd_ID')
df3$group <- 'genus'
df4 <- df2 %>% group_by(families) %>% summarize(mean_ID = mean(percent_ID, na.rm = TRUE),
                                                min_ID = min(percent_ID, na.rm = TRUE), 
                                                max_ID = max(percent_ID, na.rm = TRUE),
                                                sd_ID = sd(percent_ID, na.rm = TRUE))
colnames(df4) <- c('category', 'mean_ID', 'min_ID', 'max_ID', 'sd_ID')
df4$group <- 'family'
df5 <- df2 %>% group_by(orders) %>% summarize(mean_ID = mean(percent_ID, na.rm = TRUE),
                                              min_ID = min(percent_ID, na.rm = TRUE), 
                                              max_ID = max(percent_ID, na.rm = TRUE),
                                              sd_ID = sd(percent_ID, na.rm = TRUE))
colnames(df5) <- c('category', 'mean_ID', 'min_ID', 'max_ID', 'sd_ID')
df5$group <- 'order'
df6 <- df2 %>% group_by(classes) %>% summarize(mean_ID = mean(percent_ID, na.rm = TRUE),
                                               min_ID = min(percent_ID, na.rm = TRUE), 
                                               max_ID = max(percent_ID, na.rm = TRUE),
                                               sd_ID = sd(percent_ID, na.rm = TRUE))
colnames(df6) <- c('category', 'mean_ID', 'min_ID', 'max_ID', 'sd_ID')
df6$group <- 'class'
df7 <- df2 %>% group_by(species) %>% summarize(mean_ID = mean(percent_ID, na.rm = TRUE),
                                               min_ID = min(percent_ID, na.rm = TRUE), 
                                               max_ID = max(percent_ID, na.rm = TRUE),
                                               sd_ID = sd(percent_ID, na.rm = TRUE))
colnames(df7) <- c('category', 'mean_ID', 'min_ID', 'max_ID', 'sd_ID')
df7$group <- 'species'

df8 <- rbind(df3, df4, df5, df6, df7) 

groups <- c('class', 'order', 'family', 'genus', 'species')

ggplot(df8, aes(x = category)) + 
  geom_point(aes(y = mean_ID)) +
  geom_hline(yintercept = 95, linetype = 'dashed', color = 'darkred') +
  #geom_text(aes(x = 0.7, y = 95.5, label = "Species"),stat = "unique", color = 'darkred') +
  #geom_point(aes(y = min_ID), color = "turquoise", shape = 8, size = 4) +
  #geom_hline(yintercept = 93, linetype = 'dashed', color = 'darkgreen') +
  #geom_text(aes(x = 0.7, y = 93, label = "Family"),stat = "unique", color = 'darkgreen') +
  #geom_hline(yintercept = 92, linetype = 'dashed', color = 'darkblue') +
  #geom_text(aes(x = 0.7, y = 92, label = "Order"),stat = "unique", color = 'darkblue') +
  #geom_hline(yintercept = 91, linetype = 'dashed', color = 'darkorange3') +
  #geom_text(aes(x = 0.7, y = 91, label = "Class"),stat = "unique", color = 'darkorange3') +
  geom_errorbar(aes(ymin=mean_ID-sd_ID, ymax=mean_ID+sd_ID), width=.2, position=position_dodge(.9), col = 'black') +
  facet_grid(~factor(group, groups), scales = 'free_x') +
  scale_y_continuous(breaks = seq(80, 100, by = 1)) +
  coord_cartesian(ylim = c(80, 100), clip = 'off') +
  theme_bw() +
  ylab('Percentage identity (%) between two sequences')

species <- unique(df2$species1a)
genus <- unique(df2$genus1)
family <- unique(df2$family1)
order <- unique(df2$order1)
class <- unique(df2$class1)






