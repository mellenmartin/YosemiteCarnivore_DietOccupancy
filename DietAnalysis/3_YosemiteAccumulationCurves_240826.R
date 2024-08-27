## produce rarefraction curves for carnivore scat species diversity

#### Scats collected by Rogue Detection Teams (2019-2020)
#### Funded by OSU, YOSE NP, Yosemite Conservancy - Sean Matthews, Sarah Stock, David Green, Ben Sacks, Taal Levi 
#### Modified by Marie Martin
#### Modified from code written by Marie Tosa

# help with package from:
# https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html

rm(list=ls());gc() 
#####
#load packages
require(iNEXT)
require(vegan)
require(ggplot2)
require(tidyr) #for spread and gather

# set working drive
setwd("~/YosemiteDiet/")
diet <- read.csv("ProcessedYosemite_DietGeoDataSimplified_240414.csv") 
names(diet) = tolower(names(diet))
data.all <- diet %>% 
  dplyr::group_by(labid.y, speciesid_osu, speciessimpleidgroomed.y, species, genus, family, order, class, phylum) %>% 
  dplyr::summarize(reps = sum(replicates),
            numreads = sum(rep_reads),
            totalreads = mean(total_reads),
            pidmatch = mean(best_match_pid.mean)) %>% 
  ungroup
data.all$lineage <- paste(data.all$phylum, data.all$class, data.all$order, data.all$family, data.all$genus, data.all$species, sep=",") #make lineages line up across metabarcoding data and manual sort data
data.all$lineage <-  gsub(data.all$lineage, pattern=",{2,}", replacement="")

#sk.final <- read.table(file="Data_Raw/spilogale_gracilis_metabarcoding_manual_byphylum_2022-01-25.txt", sep=",", header=T) #NA values in rows for MT904 and F37-193

require(plyr)
#separate data by taxon
### bobcat
bobdata <- data.all[data.all$speciesid_osu == "Lynx rufus" & data.all$phylum == "Chordata",]
bobdata <- bobdata[!duplicated(bobdata[,c("lineage","labid.y")]),]
bob.wide <- spread(bobdata[,c("labid.y","lineage","numreads")], lineage, numreads, fill=0) #dims: 75 rows, 32 species, should be 128 rows
bob.wide[is.na(bob.wide)] <- 0 #replace NAs with 0
bob.wide[,2:33][bob.wide[,2:33] > 0] <- 1

### cougar
cougdata <- data.all[data.all$speciesid_osu == "Puma concolor" & data.all$phylum == "Chordata",]
cougdata <- cougdata[!duplicated(cougdata[,c("lineage","labid.y")]),]
coug.wide <- spread(cougdata[,c("labid.y","lineage","numreads")], lineage, numreads, fill=0) #dims: 75 rows, 32 species, should be 128 rows
coug.wide[is.na(coug.wide)] <- 0 #replace NAs with 0
coug.wide[,2:15][coug.wide[,2:15] > 0] <- 1

### coyote
coydata <- data.all[data.all$speciesid_osu == "Canis latrans" & data.all$phylum == "Chordata",]
coydata <- coydata[!duplicated(coydata[,c("lineage","labid.y")]),]
coy.wide <- spread(coydata[,c("labid.y","lineage","numreads")], lineage, numreads, fill=0) #dims: 75 rows, 32 species, should be 128 rows
coy.wide[is.na(coy.wide)] <- 0 #replace NAs with 0
coy.wide[,2:54][coy.wide[,2:54] > 0] <- 1

### fisher
fishdata <- data.all[data.all$speciesid_osu == "Pekania pennanti" & data.all$phylum == "Chordata",]
fishdata <- fishdata[!duplicated(fishdata[,c("lineage","labid.y")]),]
fish.wide <- spread(fishdata[,c("labid.y","lineage","numreads")], lineage, numreads, fill=0) #dims: 75 rows, 32 species, should be 128 rows
fish.wide[is.na(fish.wide)] <- 0 #replace NAs with 0
fish.wide[,2:11][fish.wide[,2:11] > 0] <- 1

### grey fox
gfoxdata <- data.all[data.all$speciesid_osu == "Urocyon cinereoargenteus" & data.all$phylum == "Chordata",]
gfoxdata <- gfoxdata[!duplicated(gfoxdata[,c("lineage","labid.y")]),]
gfox.wide <- spread(gfoxdata[,c("labid.y","lineage","numreads")], lineage, numreads, fill=0) #dims: 75 rows, 32 species, should be 128 rows
gfox.wide[is.na(gfox.wide)] <- 0 #replace NAs with 0
gfox.wide[,2:20][gfox.wide[,2:20] > 0] <- 1

### marten
martdata <- data.all[data.all$speciesid_osu == "Martes caurina" & data.all$phylum == "Chordata",]
martdata <- martdata[!duplicated(martdata[,c("lineage","labid.y")]),]
mart.wide <- spread(martdata[,c("labid.y","lineage","numreads")], lineage, numreads, fill=0) #dims: 75 rows, 32 species, should be 128 rows
mart.wide[is.na(mart.wide)] <- 0 #replace NAs with 0
mart.wide[,2:34][mart.wide[,2:34] > 0] <- 1

### red fox
rfoxdata <- data.all[data.all$speciesid_osu == "Vulpes vulpes" & data.all$phylum == "Chordata",]
rfoxdata <- rfoxdata[!duplicated(rfoxdata[,c("lineage","labid.y")]),]
rfox.wide <- spread(rfoxdata[,c("labid.y","lineage","numreads")], lineage, numreads, fill=0) #dims: 75 rows, 32 species, should be 128 rows
rfox.wide[is.na(rfox.wide)] <- 0 #replace NAs with 0
rfox.wide[,2:5][rfox.wide[,2:5] > 0] <- 1
### 
###############
#create and plot curves
data.inext <- list(
  bobcat = t(bob.wide[,!names(bob.wide) %in% c("labid.y"),]),
  cougar = t(coug.wide[,!names(coug.wide) %in% c("labid.y"),]),
  coyote = t(coy.wide[,!names(coy.wide) %in% c("labid.y"),]),
  fisher = t(fish.wide[,!names(fish.wide) %in% c("labid.y"),]),
  greyfox = t(gfox.wide[,!names(gfox.wide) %in% c("labid.y"),]),
  marten = t(mart.wide[,!names(mart.wide) %in% c("labid.y"),]),
  redfox = t(rfox.wide[,!names(rfox.wide) %in% c("labid.y"),]))
bob.inext <- list(bobcat = t(bob.wide[,!names(bob.wide) %in% c("labid.y"),]))
coug.inext <- list(cougar = t(coug.wide[,!names(coug.wide) %in% c("labid.y"),]))
coy.inext <- list(coyote = t(coy.wide[,!names(coy.wide) %in% c("labid.y"),]))
fish.inext <- list(fisher = t(fish.wide[,!names(fish.wide) %in% c("labid.y"),]))
grey.inext <- list(greyfox = t(gfox.wide[,!names(gfox.wide) %in% c("labid.y"),]))
mart.inext <- list(marten = t(mart.wide[,!names(mart.wide) %in% c("labid.y"),]))
red.inext <- list(redfox = t(rfox.wide[,!names(rfox.wide) %in% c("labid.y"),]))

out.inext <- iNEXT(x=data.inext, datatype="incidence_raw", endpoint=200) #data needs to be rows = species, columns = samples

bob.inext <- iNEXT(x=bob.inext, datatype="incidence_raw", endpoint=200) #data needs to be rows = species, columns = samples
coug.inext <- iNEXT(x=coug.inext, datatype="incidence_raw", endpoint=200) #data needs to be rows = species, columns = samples
coy.inext <- iNEXT(x=coy.inext, datatype="incidence_raw", endpoint=800) #data needs to be rows = species, columns = samples
fish.inext <- iNEXT(x=fish.inext, datatype="incidence_raw", endpoint=100) #data needs to be rows = species, columns = samples
grey.inext <- iNEXT(x=grey.inext, datatype="incidence_raw", endpoint=150) #data needs to be rows = species, columns = samples
mart.inext <- iNEXT(x=mart.inext, datatype="incidence_raw", endpoint=200) #data needs to be rows = species, columns = samples
red.inext <- iNEXT(x=red.inext, datatype="incidence_raw", endpoint=100) #data needs to be rows = species, columns = samples

p <- ggiNEXT(out.inext, facet.var = "Assemblage", color.var = "Assemblage",type=1) + 
  xlab("Number of scat samples") + 
  ylab("Taxonomic richness") +
  scale_color_manual(values=c("#0288D1"))+
  scale_fill_manual(values = c("#0288D1"))+
  #theme(axis.title = element_text(size = 16))+
  #theme(legend.position="right")+
  #theme(legend.title = element_text(colour="black", size=16))+
  #theme(legend.text = element_text(colour="black", size = 16))+
  #theme(axis.text=element_text(size=10))
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))
p
ggsave("prelimresults/YosemiteComm_prelimaccu_230208.png", width=20,height=6, units = 'in', dpi= 150)
# ggsave(p, filename="Figures/Rarefaction_Spilogale_gracilis.tiff", height=6, width=8, units="in", dpi=300, compression="lzw")

hcl.colors(7, palette = "ag_GrnYl")
"#255668" "#007178" "#008C80" "#17A77E" "#61C074" "#A4D764" "#EDEF5C"


bob<- ggiNEXT(bob.inext, color.var = "Order.q", type=1) + 
  xlab("") + 
  ylab("Taxonomic richness") +
  xlim(0,200)+
  scale_color_manual(values=c("#4B0055"), guide = "none")+
  scale_fill_manual(values = c("#4B0055"), guide = "none")+
  ggtitle("Bobcat")+
  #theme(axis.title = element_text(size = 16))+
  #theme(legend.position="right")+
  #theme(legend.title = element_text(colour="black", size=16))+
  #theme(legend.text = element_text(colour="black", size = 16))+
  #theme(axis.text=element_text(size=10))
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))


coug<- ggiNEXT(coug.inext, color.var = "Order.q", type=1) + 
  xlab("") + 
  ylab("") +
  xlim(0,200)+
  scale_color_manual(values=c("#00588B"), guide = "none")+
  scale_fill_manual(values = c("#00588B"), guide = "none")+
  ggtitle("Cougar")+
  #theme(axis.title = element_text(size = 16))+
  #theme(legend.position="right")+
  #theme(legend.title = element_text(colour="black", size=16))+
  #theme(legend.text = element_text(colour="black", size = 16))+
  #theme(axis.text=element_text(size=10))
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))

coy<- ggiNEXT(coy.inext, color.var = "Order.q", type=1) + 
  xlab("") + 
  ylab("") +
  xlim(0,800)+
  scale_color_manual(values=c("#009B95"), guide = "none")+
  scale_fill_manual(values = c("#009B95"), guide = "none")+
  ggtitle("Coyote")+
  #theme(axis.title = element_text(size = 16))+
  #theme(legend.position="right")+
  #theme(legend.title = element_text(colour="black", size=16))+
  #theme(legend.text = element_text(colour="black", size = 16))+
  #theme(axis.text=element_text(size=10))
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))

grey<- ggiNEXT(grey.inext, color.var = "Order.q", type=1) + 
  xlab("Number of samples") + 
  ylab("Taxonomic richness") +
  xlim(0,150)+
  scale_color_manual(values=c("#53cc67"), guide = "none")+
  scale_fill_manual(values = c("#53cc67"), guide = "none")+
  ggtitle("Grey fox")+
  #theme(axis.title = element_text(size = 16))+
  #theme(legend.position="right")+
  #theme(legend.title = element_text(colour="black", size=16))+
  #theme(legend.text = element_text(colour="black", size = 16))+
  #theme(axis.text=element_text(size=10))
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))

mart<- ggiNEXT(mart.inext, color.var = "Order.q", type=1) + 
  xlab("Number of samples") + 
  ylab(" ") +
  xlim(0,200)+
  scale_color_manual(values=c("#FDE333"), guide = "none")+
  scale_fill_manual(values = c("#FDE333"), guide = "none")+
  ggtitle("Marten")+
  #theme(axis.title = element_text(size = 16))+
  #theme(legend.position="right")+
  #theme(legend.title = element_text(colour="black", size=16))+
  #theme(legend.text = element_text(colour="black", size = 16))+
  #theme(axis.text=element_text(size=10))
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))
mart
ggsave(ggarrange(bob,
                 coug,
                 coy,
                 grey,
                 mart,
                 nrow=2, ncol = 3), 
       filename=paste("prelimresults/YosemiteDiet_AccuCurves_", Sys.Date(), ".jpg", sep=""), 
       height=6, width=12, units="in", dpi=300)
