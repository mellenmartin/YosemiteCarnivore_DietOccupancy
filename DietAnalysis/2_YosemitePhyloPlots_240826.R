###################
##=============================
# Diet data - species level plots
#=============================


#### Scats collected by Rogue Detection Teams (2019-2020)
#### Funded by  OSU, YOSE NP, Yosemite Conservancy - Matthews, Stock, Green, Levi 
#### Last modified 230118, Marie Martin
#### Modified from Brent Barry PEPE Diet Prelim  # BRB 29042022   
###################

#library(plyr)
library(dplyr)
library(ggplot2)
remove(list=ls()) 
gc()
### Step 1, Load from csv
setwd("~/YosemiteDiet/")
diet <- read.csv("ProcessedYosemite_DietGeoDataSimplified_240414.csv") 
names(diet) = tolower(names(diet))

##### phylo plots #####
require(dplyr)
require(reshape2)
require(metacoder)
require(ggplot2)
require(ggpubr)
require(lubridate)
require(plyr)

#create functions for metacoder
plots <- function(season, meta.obj, phylum, tree_color_range, seednum, nobs = 0)
{
  set.seed(seednum)
  temp <- meta.obj %>%
    filter_taxa(n_obs > nobs, reassign_obs = F)
  ht <- heat_tree(temp,
                  #node_size=n_obs,
                  node_size = temp$data$carnivore_occ[["total"]],
                  node_size_range = c(0.001, 0.01),
                  node_color = temp$data$carnivore_occ[[season]],
                  node_color_range = tree_color_range,
                  edge_size_range = c(0.005, 0.015),
                  tree_color_range = tree_color_range,
                  layout="davidson-harel",
                  title = phylum,
                  node_color_axis_label="focc",
                  node_size_axis_label="total focc")
  return(ht)
}

key_plot <- function(meta.obj, phylum, tree_color_range, seednum, nobs = 0)
{
  set.seed(seednum)
  temp <- meta.obj %>%
    filter_taxa(n_obs > nobs, reassign_obs = F)
  ht <- heat_tree(temp,
                  node_label = taxon_names, 
                  node_label_size_range = c(0.02, 0.05),
                  #node_size=n_obs,
                  node_size = temp$data$carnivore_occ[["total"]],
                  node_size_range = c(0.005, 0.015),
                  node_color = temp$data$carnivore_occ[["total"]],
                  node_color_range = tree_color_range,
                  edge_size_range = c(0.005, 0.01),
                  tree_color_range = tree_color_range,
                  overlap_avoidance = 2,
                  margin_size = c(0.01, 0.01, 0, 0),
                  title = phylum,
                  layout="davidson-harel",
                  node_size_axis_label = " ",
                  node_color_axis_label="No. samples",
                  repel_labels = TRUE)
  return(ht)
}

############
#load data
library(dplyr)
detach(package:plyr)
#all data
#data.all <- read.table("Data_Raw/dataall_meta_man_2022-04-18.txt", sep=",") #data with metabarcoding and manual sorting
data.all2 <- diet %>% dplyr::group_by(labid.y, speciesid_osu, speciessimpleidgroomed.y, species, genus, family, order, class, phylum) %>% 
  dplyr::filter(!speciesid_osu == "Buteo spp.",
               !speciesid_osu == "No Amp",
               !speciesid_osu == "Unknown") %>% 
  dplyr::summarize(reps = sum(replicates),
                      numreads = sum(rep_reads),
                      totalreads = mean(total_reads),
                      pidmatch = mean(best_match_pid.mean)) %>% ungroup
  
#sk.final <- read.table(file="Data_Raw/spilogale_gracilis_metabarcoding_manual_byphylum_2022-01-25.txt", sep=",", header=T) #NA values in rows for MT904 and F37-193, samples only had clupea
#sk.final <- sk.final[!sk.final$ScatID %in% c("MT904","F37-193"),] #remove MT904 and F37-193 from sk.final, samples only had clupea

#give data.all new season information
#data.all <- merge(data.all, sk.final[,c("ScatID","CollectionType","CollectionDate","scatseason","logged","e","e.m","YrsSinceDi")], by.x="sample", by.y="ScatID", all.x=T)

data.all2$lineage <- paste(data.all2$phylum, data.all2$class, data.all2$order, data.all2$family, data.all2$genus, data.all2$species, sep=",") #make lineages line up across metabarcoding data and manual sort data
data.all2$lineage <-  gsub(data.all2$lineage, pattern=",{2,}", replacement="")
data.all2$scatseason <- rep("summer", length(data.all2$labid.y))

write.csv(data.all2, "ProcessedYosemite_PhyloSummaryData_240417.csv", row.names = FALSE)
##### marten #####
#vertebrate data
data.all <- data.all2
vertdata2 <- data.all[data.all$phylum == "Chordata",]

library(plyr)
library(ddply)
#summarize number of verts per scat
vertcount2 <- ddply(vertdata2, .(labid.y), nrow)
#
samples <- unique(vertdata2[,c("labid.y")])

sk2 <- vertdata2 %>% dplyr::select(labid.y, lineage)
sk2 <- melt(sk2, id=c("labid.y", "lineage"))
sk2 <- dcast(data=sk2, lineage~labid.y, fun.aggregate=length)

obj <- parse_tax_data(sk2, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata2$labid.y))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata2$labid.y), groups=unique(vertdata2[,c("labid.y","speciesid_osu")])$speciesid_osu)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$`Canis latrans`+ 
  obj$data$carnivore_occ$`Martes caurina`+ 
  obj$data$carnivore_occ$`Puma concolor`+
  obj$data$carnivore_occ$`Lynx rufus`+
  obj$data$carnivore_occ$`Urocyon cinereoargenteus`+
  obj$data$carnivore_occ$`Ursus americanus`+
  obj$data$carnivore_occ$`Vulpes vulpes`+
  obj$data$carnivore_occ$`Pekania pennanti`

v_color <- c("#36728E", "#498EA4", "#67A9B6","#88C3C8", "#D2EEEA")

marten <- plots("Martes caurina", meta.obj=obj, phylum="", tree_color_range=v_color, seednum=2)
marten
vkey <- key_plot(meta.obj=obj, phylum="Yosemite prey community", tree_color_range=v_color, seednum=2)
vkey
ggsave("prelimresults/YosemiteDiet_AllPhylo_240417.jpg", width=20,height=14, units = 'in', dpi= 300)


##### puma #####
#vertebrate data
vertdata2 <- data.all2[data.all2$speciesid_osu == "Puma concolor" & data.all2$phylum == "Chordata",]
#remove MT906 and MT907 because probably clupea bones via manual sort
#vertdata <- vertdata[!vertdata$sample %in% c("MT906","MT907"),]
library(plyr)
#summarize number of verts per scat
vertcount2 <- ddply(vertdata2, .(labid.y), nrow)
#
martsamples <- unique(vertdata2[,c("labid.y")])

sk2 <- vertdata2 %>% select(labid.y, lineage)
sk2 <- melt(sk2, id=c("labid.y", "lineage"))
sk2 <- dcast(data=sk2, lineage~labid.y, fun.aggregate=length)

obj <- parse_tax_data(sk2, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata2$labid.y))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata2$labid.y), groups=unique(vertdata2[,c("labid.y","scatseason")])$scatseason)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$summer 

v_color <- c("azure3" , "sienna3" , "darkcyan")

f <- plots("summer", meta.obj=obj, phylum="Cougar (n = 98)", tree_color_range=v_color, seednum=2)
vkey <- key_plot(meta.obj=obj, phylum="Cougar (n = 98)", tree_color_range=v_color, seednum=2)
vkey
ggsave("prelimresults/YosemitePUCO_prelimphylo_230208.png", width=11,height=11, units = 'in', dpi= 150)

##### coyote #####
#vertebrate data
vertdata2 <- data.all2[data.all2$speciesid_osu == "Canis latrans" & data.all2$phylum == "Chordata",]
#remove MT906 and MT907 because probably clupea bones via manual sort
#vertdata <- vertdata[!vertdata$sample %in% c("MT906","MT907"),]
library(plyr)
#summarize number of verts per scat
vertcount2 <- ddply(vertdata2, .(labid.y), nrow)
#
martsamples <- unique(vertdata2[,c("labid.y")])

sk2 <- vertdata2 %>% select(labid.y, lineage)
sk2 <- melt(sk2, id=c("labid.y", "lineage"))
sk2 <- dcast(data=sk2, lineage~labid.y, fun.aggregate=length)

obj <- parse_tax_data(sk2, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata2$labid.y))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata2$labid.y), groups=unique(vertdata2[,c("labid.y","scatseason")])$scatseason)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$summer 

v_color <- c("azure3" , "sienna3" , "darkcyan")

f <- plots("summer", meta.obj=obj, phylum="Coyote (n = 707)", tree_color_range=v_color, seednum=2)
vkey <- key_plot(meta.obj=obj, phylum="Coyote (n = 707)", tree_color_range=v_color, seednum=2)
vkey
ggsave("prelimresults/YosemiteCALA_prelimphylo_230208.png", width=11,height=11, units = 'in', dpi= 150)

##### bobcat #####
#vertebrate data
vertdata2 <- data.all2[data.all2$speciesid_osu == "Lynx rufus" & data.all2$phylum == "Chordata",]
#remove MT906 and MT907 because probably clupea bones via manual sort
#vertdata <- vertdata[!vertdata$sample %in% c("MT906","MT907"),]
library(plyr)
#summarize number of verts per scat
vertcount2 <- ddply(vertdata2, .(labid.y), nrow)
#
martsamples <- unique(vertdata2[,c("labid.y")])

sk2 <- vertdata2 %>% select(labid.y, lineage)
sk2 <- melt(sk2, id=c("labid.y", "lineage"))
sk2 <- dcast(data=sk2, lineage~labid.y, fun.aggregate=length)

obj <- parse_tax_data(sk2, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata2$labid.y))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata2$labid.y), groups=unique(vertdata2[,c("labid.y","scatseason")])$scatseason)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$summer 

v_color <- c("azure3" , "sienna3" , "darkcyan")

f <- plots("summer", meta.obj=obj, phylum="Bobcat (n = 101)", tree_color_range=v_color, seednum=2)
vkey <- key_plot(meta.obj=obj, phylum="Bobcat (n = 101)", tree_color_range=v_color, seednum=2)
vkey
ggsave("prelimresults/YosemiteLYRU_prelimphylo_230208.png", width=11,height=11, units = 'in', dpi= 150)

##### grey fox #####
#vertebrate data
vertdata2 <- data.all2[data.all2$speciesid_osu == "Urocyon cinereoargenteus" & data.all2$phylum == "Chordata",]
#remove MT906 and MT907 because probably clupea bones via manual sort
#vertdata <- vertdata[!vertdata$sample %in% c("MT906","MT907"),]
library(plyr)
#summarize number of verts per scat
vertcount2 <- ddply(vertdata2, .(labid.y), nrow)
#
martsamples <- unique(vertdata2[,c("labid.y")])

sk2 <- vertdata2 %>% select(labid.y, lineage)
sk2 <- melt(sk2, id=c("labid.y", "lineage"))
sk2 <- dcast(data=sk2, lineage~labid.y, fun.aggregate=length)

obj <- parse_tax_data(sk2, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata2$labid.y))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata2$labid.y), groups=unique(vertdata2[,c("labid.y","scatseason")])$scatseason)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$summer 

v_color <- c("azure3" , "sienna3" , "darkcyan")

f <- plots("summer", meta.obj=obj, phylum="Grey fox (n = 43)", tree_color_range=v_color, seednum=2)
vkey <- key_plot(meta.obj=obj, phylum="Grey fox (n = 43)", tree_color_range=v_color, seednum=2)
vkey
ggsave("prelimresults/YosemiteURCI_prelimphylo_230208.png", width=11,height=11, units = 'in', dpi= 150)

##### fisher #####
#vertebrate data
vertdata2 <- data.all2[data.all2$speciesid_osu == "Pekania pennanti" & data.all2$phylum == "Chordata",]
#remove MT906 and MT907 because probably clupea bones via manual sort
#vertdata <- vertdata[!vertdata$sample %in% c("MT906","MT907"),]
library(plyr)
#summarize number of verts per scat
vertcount2 <- ddply(vertdata2, .(labid.y), nrow)
#
martsamples <- unique(vertdata2[,c("labid.y")])

sk2 <- vertdata2 %>% select(labid.y, lineage)
sk2 <- melt(sk2, id=c("labid.y", "lineage"))
sk2 <- dcast(data=sk2, lineage~labid.y, fun.aggregate=length)

obj <- parse_tax_data(sk2, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata2$labid.y))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata2$labid.y), groups=unique(vertdata2[,c("labid.y","scatseason")])$scatseason)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$summer 

v_color <- c("azure3" , "sienna3" , "darkcyan")

f <- plots("summer", meta.obj=obj, phylum="Fisher (n = 12)", tree_color_range=v_color, seednum=2)
vkey <- key_plot(meta.obj=obj, phylum="Fisher (n = 12)", tree_color_range=v_color, seednum=2)
vkey
ggsave("prelimresults/YosemitePEPE_prelimphylo_230208.png", width=11,height=11, units = 'in', dpi= 150)

##### red fox #####
#vertebrate data
vertdata2 <- data.all2[data.all2$speciesid_osu == "Vulpes vulpes" & data.all2$phylum == "Chordata",]
#remove MT906 and MT907 because probably clupea bones via manual sort
#vertdata <- vertdata[!vertdata$sample %in% c("MT906","MT907"),]
library(plyr)
#summarize number of verts per scat
vertcount2 <- ddply(vertdata2, .(labid.y), nrow)
#
martsamples <- unique(vertdata2[,c("labid.y")])

sk2 <- vertdata2 %>% select(labid.y, lineage)
sk2 <- melt(sk2, id=c("labid.y", "lineage"))
sk2 <- dcast(data=sk2, lineage~labid.y, fun.aggregate=length)

obj <- parse_tax_data(sk2, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata2$labid.y))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata2$labid.y), groups=unique(vertdata2[,c("labid.y","scatseason")])$scatseason)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$summer 

v_color <- c("#36728E", "#498EA4", "#67A9B6","#88C3C8", "#D2EEEA")

f <- plots("summer", meta.obj=obj, phylum="Red fox (n = 4)", tree_color_range=v_color, seednum=2)
vkey <- key_plot(meta.obj=obj, phylum="Red fox (n = 4)", tree_color_range=v_color, seednum=2)
vkey
ggsave("prelimresults/YosemiteVUVU_prelimphylo_230208.png", width=11,height=11, units = 'in', dpi= 150)

