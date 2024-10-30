
library(dplyr)
library(tidyr)
library(tidyverse)
dir()
#######################
### Step 1, Load from csv
setwd("~/YosemiteDiet/")
sealD <- read.csv("ProcessedYosemite_DietGeoDataSimplified_241002.csv")
head(sealD)
sealD <- sealD %>% filter(!defecator == "Neogale frenata",
                          !defecator == "Buteo",
                          !defecator == "Ursus americanus",
                          !defecator == "Vulpes vulpes",
                          !defecator == "Marmota flaviventris",
                          !defecator == "Pekania pennanti",
                          !defecator == "Odocoileus",
                          !idlevel=="Family",
                          !idlevel == "Order",
                          !idlevel=="Class") 

evalq(tapply(labid.y,defecator,function(x) length(unique(x))),sealD)

#  Number of prey species
length(unique(sealD$species))

# Per sample
n.tps=evalq(tapply(speciessimpleidgroomed.y,labid.y,function(x) length(unique(x))),sealD)
hist(n.tps)
mean(n.tps)

seq.depth=evalq(tapply(rep_reads,labid.y,sum),sealD)
summary(seq.depth)

sum(seq.depth<50)
names(seq.depth[seq.depth<50])
hist(seq.depth)
# which samples have less than 50 reads
sealD$sample %in% names(seq.depth[seq.depth<50])

#sum(seq.depth50<50)
sealD50 <- sealD

sealD50$defcommonname <- ifelse(sealD50$defecator=="Martes caurina", "Marten",
                              ifelse(sealD50$defecator == "Canis latrans", "Coyote",
                                               ifelse(sealD50$defecator == "Puma concolor", "Cougar",
                                                      ifelse(sealD50$defecator == "Lynx rufus", "Bobcat", "Grayfox"))))

#14 Strata - removing 2 with low samples
Marten = evalq(defcommonname =="Marten", sealD50)
Cougar = evalq(defcommonname=="Cougar", sealD50)
Bobcat = evalq(defcommonname=="Bobcat", sealD50)
Coyote = evalq(defcommonname=="Coyote", sealD50)
Grayfox = evalq(defcommonname=="Grayfox", sealD50)

IndexNAme= c(unique(sealD50$defcommonname))

length(IndexNAme) # This is how many subsets of data we are using
# Create subsets
t=0
for( i in IndexNAme){
  t=t+1
  nam =paste("Subset",formatC(t,width = 2, format = "d", flag = "0"),sep="")
  temp=evalq(tapply(rep_reads,list(speciessimpleidgroomed.y,labid.y),sum),sealD50[eval(as.name(i)),])
  temp[is.na(temp)]=0
  temp50=  temp[,apply(temp,2,sum)>50]
  assign(nam,temp50 )
}
## Look at one subset in detail
Subset01 #first subset
dim(Subset01) # number of prey and number of samples

#### marten ####
#Look at the counts of prey identified in this subset
Marten=evalq(tapply(rep_reads,list(speciessimpleidgroomed.y,labid.y),sum),sealD50[Marten,])
Marten[is.na(Marten)]=0
apply(Marten,1,sum)
# Convert to proportions
Martenprop=prop.table(Marten,2)
dim(Martenprop)
apply(Marten,1,sum)
# No prey occurences are less than 1%
sum(Martenprop < .01)
sum(Martenprop== 0)
# Make presence/absence dataset
PA_Marten= Marten; PA_Marten[Marten>0]=1
apply(PA_Marten,1,sum)
#POO summary
dim(PA_Marten)# 66 samples
apply(PA_Marten,1,sum) # Number of occurences
POO=apply(PA_Marten,1,sum)/sum(apply(PA_Marten,1,sum))*100
apply(PA_Marten,2,sum)# how many prey per scat
#wPOO summary
dim(PA_Marten)# 66 samples
wPOO=apply(prop.table(PA_Marten,2),1,mean)*100
#RRA summary
RRA=apply(prop.table(as.matrix(Marten),2),1,mean)*100
MartenSum <- as.data.frame(cbind(POO, RRA, wPOO))
MartenSum <- tibble::rownames_to_column(as.data.frame(MartenSum), "PreyTaxa")
MartenSum$Species <- rep("Marten", nrow(MartenSum))
MartenSum$reads <- rowSums(Marten)

#### coyote ####
#Look at the counts of prey identified in this subset
Coyote=evalq(tapply(rep_reads,list(speciessimpleidgroomed.y,labid.y),sum),sealD50[Coyote,])
Coyote[is.na(Coyote)]=0
apply(Coyote,1,sum)
# Convert to proportions
Coyoteprop=prop.table(Coyote,2)
dim(Coyoteprop)
apply(Coyote,1,sum)
# No prey occurences are less than 1%
sum(Coyoteprop < .01)
sum(Coyoteprop== 0)
# Make presence/absence dataset
PA_Coyote= Coyote; PA_Coyote[Coyote>0]=1
apply(PA_Coyote,1,sum)
#POO summary
dim(PA_Coyote)
apply(PA_Coyote,1,sum) # Number of occurences
POO=apply(PA_Coyote,1,sum)/sum(apply(PA_Coyote,1,sum))*100
apply(PA_Coyote,2,sum)# how many prey per scat
#wPOO summary
dim(PA_Coyote)
wPOO=apply(prop.table(PA_Coyote,2),1,mean)*100
#RRA summary
RRA=apply(prop.table(as.matrix(Coyote),2),1,mean)*100
CoyoteSum <- as.data.frame(cbind(POO, RRA, wPOO))
CoyoteSum <- tibble::rownames_to_column(as.data.frame(CoyoteSum), "PreyTaxa")
CoyoteSum$Species <- rep("Coyote", nrow(CoyoteSum))
CoyoteSum$reads <- rowSums(Coyote)

##### Cougar #####
#Look at the counts of prey identified in this subset
Cougar=evalq(tapply(rep_reads,list(speciessimpleidgroomed.y,labid.y),sum),sealD50[Cougar,])
Cougar[is.na(Cougar)]=0
apply(Cougar,1,sum)
# Convert to proportions
Cougarprop=prop.table(Cougar,2)
dim(Cougarprop)
apply(Cougar,1,sum)
# No prey occurences are less than 1%
sum(Cougarprop < .01)
sum(Cougarprop== 0)
# Make presence/absence dataset
PA_Cougar= Cougar; PA_Cougar[Cougar>0]=1
apply(PA_Cougar,1,sum)
#POO summary
dim(PA_Cougar)
apply(PA_Cougar,1,sum) # Number of occurences
POO=apply(PA_Cougar,1,sum)/sum(apply(PA_Cougar,1,sum))*100
apply(PA_Cougar,2,sum)# how many prey per scat
#wPOO summary
dim(PA_Cougar)
wPOO=apply(prop.table(PA_Cougar,2),1,mean)*100
#RRA summary
#RRA summary
RRA=apply(prop.table(as.matrix(Cougar),2),1,mean)*100
CougarSum <- as.data.frame(cbind(POO, RRA, wPOO))
CougarSum <- tibble::rownames_to_column(as.data.frame(CougarSum), "PreyTaxa")
CougarSum$Species <- rep("Cougar", nrow(CougarSum))
CougarSum$reads <- rowSums(Cougar)
  
##### grey fox #####
#Look at the counts of prey identified in this subset
Grayfox=evalq(tapply(rep_reads,list(speciessimpleidgroomed.y,labid.y),sum),sealD50[Grayfox,])
Grayfox[is.na(Grayfox)]=0
apply(Grayfox,1,sum)
# Convert to proportions
Grayfoxprop=prop.table(Grayfox,2)
dim(Grayfoxprop)
apply(Grayfox,1,sum)
# No prey occurences are less than 1%
sum(Grayfoxprop < .01)
sum(Grayfoxprop== 0)
# Make presence/absence dataset
PA_Grayfox= Grayfox; PA_Grayfox[Grayfox>0]=1
apply(PA_Grayfox,1,sum)
#POO summary
dim(PA_Grayfox)
apply(PA_Grayfox,1,sum) # Number of occurences
POO=apply(PA_Grayfox,1,sum)/sum(apply(PA_Grayfox,1,sum))*100
apply(PA_Grayfox,2,sum)# how many prey per scat
#wPOO summary
dim(PA_Grayfox)
wPOO=apply(prop.table(PA_Grayfox,2),1,mean)*100
#RRA summary
RRA=apply(prop.table(as.matrix(Grayfox),2),1,mean)*100
GrayfoxSum <- as.data.frame(cbind(POO, RRA, wPOO))
GrayfoxSum <- tibble::rownames_to_column(as.data.frame(GrayfoxSum), "PreyTaxa")
GrayfoxSum$Species <- rep("Grayfox", nrow(GrayfoxSum))
GrayfoxSum$reads <- rowSums(Grayfox)

##### grey fox #####
#Look at the counts of prey identified in this subset
Bobcat=evalq(tapply(rep_reads,list(speciessimpleidgroomed.y,labid.y),sum),sealD50[Bobcat,])
Bobcat[is.na(Bobcat)]=0
apply(Bobcat,1,sum)
# Convert to proportions
Bobcatprop=prop.table(Bobcat,2)
dim(Bobcatprop)
apply(Bobcat,1,sum)
# No prey occurences are less than 1%
sum(Bobcatprop < .01)
sum(Bobcatprop== 0)
# Make presence/absence dataset
PA_Bobcat= Bobcat; PA_Bobcat[Bobcat>0]=1
apply(PA_Bobcat,1,sum)
#POO summary
dim(PA_Bobcat)
apply(PA_Bobcat,1,sum) # Number of occurences
POO=apply(PA_Bobcat,1,sum)/sum(apply(PA_Bobcat,1,sum))*100
apply(PA_Bobcat,2,sum)# how many prey per scat
#wPOO summary
dim(PA_Bobcat)
wPOO=apply(prop.table(PA_Bobcat,2),1,mean)*100
#RRA summary
RRA=apply(prop.table(as.matrix(Bobcat),2),1,mean)*100
BobcatSum <- as.data.frame(cbind(POO, RRA, wPOO))
BobcatSum <- tibble::rownames_to_column(as.data.frame(BobcatSum), "PreyTaxa")
BobcatSum$Species <- rep("Bobcat", nrow(BobcatSum))
BobcatSum$reads <- rowSums(Bobcat)

AllSum <- bind_rows(BobcatSum,
                CougarSum,
                CoyoteSum,
                GrayfoxSum,
                MartenSum)

write.csv(AllSum, "YosemiteDiet_AllspeciesSummaryFOOwPOO_241003.csv", row.names = FALSE)
AllSum <- read.csv("YosemiteDiet_AllspeciesSummaryFOOwPOO_241003.csv")
#########################################################

Subset01 <- tibble::rownames_to_column(as.data.frame(Subset01), "PreyTaxa")
write.csv(Subset01, "YosemiteDiet_GroomedMartenPreyMatrix_240418.csv", row.names = FALSE)
Subset02 <- tibble::rownames_to_column(as.data.frame(Subset02), "PreyTaxa")
write.csv(Subset02, "YosemiteDiet_GroomedCoyotePreyMatrix_240418.csv", row.names = FALSE)
Subset03 <- tibble::rownames_to_column(as.data.frame(Subset03), "PreyTaxa")
write.csv(Subset03, "YosemiteDiet_GroomedCougarPreyMatrix_240418.csv", row.names = FALSE)
Subset04 <- tibble::rownames_to_column(as.data.frame(Subset04), "PreyTaxa")
write.csv(Subset04, "YosemiteDiet_GroomedGrayfoxPreyMatrix_240418.csv", row.names = FALSE)
Subset05 <- tibble::rownames_to_column(as.data.frame(Subset05), "PreyTaxa")
write.csv(Subset05, "YosemiteDiet_GroomedBobcatPreyMatrix_240418.csv", row.names = FALSE)

#### plotweb ###
sizecats <- read.csv("YOSE_SpeciesSummary_SizeTaxonomy_230703.csv")
allsum_fullinfo <- left_join(AllSum, sizecats, by = c("PreyTaxa" = "speciessimpleidgroomed"))
write.csv(allsum_fullinfo, "YosemiteDiet_AllSpeciesAllInfoSummary_241003.csv", row.names = FALSE)
allsum_fullinfo <- read.csv("YosemiteDiet_AllSpeciesAllInfoSummary_241003.csv")
allsum_pivotclass <- allsum_fullinfo %>% 
  group_by(Species, class) %>% 
  dplyr::summarise(classwPOO = sum(wPOO)) %>% 
  dplyr::select(c("class", "classwPOO", "Species")) %>% 
  pivot_wider(names_from = class, values_from = classwPOO) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames(var="Species")

allsum_class <- allsum_fullinfo %>% 
  dplyr::group_by(Species, class) %>% 
  dplyr::summarise(classwPOO = sum(wPOO)) %>% 
  dplyr::select(c("class", "classwPOO", "Species"))
vert.colors <- c("#4B0055","#353E7C", "#007094","#009B95","#00BE7D")
classwpoo <- ggplot(allsum_class, aes(x=Species, y=classwPOO, fill=class)) + 
  #rra <- ggplot(verttab, aes(x=sciname, y=finalrra_all, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) +
  scale_fill_manual(name="Prey class", values=vert.colors) + 
  ylab("Weighted percent of occurrence") +
  xlab("")+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, size = 18))

classwpoo

allsum_group <- allsum_fullinfo %>% 
  dplyr::group_by(Species, taxasizegroup) %>% 
  dplyr::summarise(groupwPOO = sum(wPOO)) %>% 
  dplyr::select(c("taxasizegroup", "groupwPOO", "Species")) %>% 
  dplyr::filter(!taxasizegroup == "Reptiles")

unique(allsum_group$taxasizegroup)

vert.colors <- hcl.colors(9, palette = "tofino")
groupwpoo <- ggplot(allsum_group, aes(x=Species, y=groupwPOO, fill=taxasizegroup)) + 
  #rra <- ggplot(verttab, aes(x=sciname, y=finalrra_all, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) +
  scale_fill_manual(name="Prey group", values=rev(vert.colors)) + 
  ylab("Weighted percent of occurrence") +
  xlab("")+
  theme(#strip.background = element_blank(),
    axis.text = element_text(size=20),
    #axis.text.y = element_text(angle = 90, size = 18),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22, angle = 90),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    strip.text.x=element_text(size=20, margin = margin(2,2,2,2)),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18))

groupwpoo
ggsave("prelimresults/YosemiteCarn_PreyGroup_240418.jpg", width=8,height=6, units = 'in', dpi= 600)

#### bray curtis and piankas ####
allsum_pivotspecies <- allsum_fullinfo %>% 
  dplyr::select(c("PreyTaxa", "wPOO", "Species")) %>% 
  pivot_wider(names_from = PreyTaxa, values_from = wPOO) %>% 
  replace(is.na(.), 0) 
#%>% 
 # column_to_rownames(var="Species")

carn.dist <- vegdist(allsum_pivotspecies[,2:66], method="bray", diag=TRUE, upper=TRUE)
carn.dist
set.seed(36) #reproducible results
carn.div<-adonis2(carn.dist~Species, data=allsum_pivotspecies, permutations=1000)
carn.div
summary(carn.div)

library(spaa)
library(sjmisc)
allsum_pivotspecieslong <- allsum_pivotspecies %>%
  column_to_rownames(var="Species") %>% rotate_df()
niche.overlap(allsum_pivotspecieslong, method = "pianka")

library(EcoSimR)
allsum_pivotspecies <- allsum_pivotspecies %>% column_to_rownames(var = "Species")
nullmod <- niche_null_model(allsum_pivotspecies,nReps=1000)
nullmod

## Summary and plot info
summary(nullmod)
plot(nullmod,type="niche")
niche.overlap.boot(allsum_pivotspecieslong,method="pianka",times=999,quant=c(0.025,0.975))

#### nmds and hull plots ##### 
allsum_pivotspecies <- allsum_fullinfo %>% 
  dplyr::select(c("PreyTaxa", "wPOO", "Species")) %>% 
  pivot_wider(names_from = PreyTaxa, values_from = wPOO) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames(var="Species")


#make community matrix - extract columns with abundance information
com = allsum_pivotspecies[,1:ncol(allsum_pivotspecies)]
#turn abundance data frame into a matrix
m_com = as.matrix(com)

nmds_results <- metaMDS(comm = m_com, distance = "bray", try = 500) 
plot(nmds_results)

data.scores <- as.data.frame(scores(nmds_results)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
species.scores <- as.data.frame(scores(nmds_results, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

nmds_hulljoin <- left_join(AllSum, species.scores, by = c("PreyTaxa" = "species"))
nmds_hulljoin$family <- ifelse(nmds_hulljoin$Defecator=="Marten", "Mustelid",
                                ifelse(nmds_hulljoin$Defecator == "Coyote", "Canid",
                                       ifelse(nmds_hulljoin$Defecator == "Cougar", "Felid",
                                              ifelse(nmds_hulljoin$Defecator == "Bobcat", "Felid", "Canid"))))

ordiellipse(nmds_results, rownames(data.scores), display = "species", kind = "se", label = T)

#originalresponse
NMDS = data.frame(MDS1 = nmds_hulljoin$NMDS1, MDS2=nmds_hulljoin$NMDS2, Species= nmds_hulljoin$Defecator)
NMDS.mean=aggregate(NMDS[,1:2],list(Species=NMDS$Species),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 50) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in unique(NMDS$Species)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Species==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,Species=g))
}

nmds_ellipse <- ggplot(data = NMDS, aes(MDS1, MDS2)) + 
  geom_point(shape = 3, size = 3) +
  #geom_jitter(aes(color = Species), size = 3, width = 0.01)+
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=Species), size=2, linetype=2)+
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$Species, size = 8)+
  labs(y= "NMDS2", x = "NMDS1") +
  theme(#strip.background = element_blank(),
    axis.text = element_text(size=20),
    #axis.text.y = element_text(angle = 90, size = 18),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22, angle = 90),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text.x=element_text(size=20, margin = margin(2,2,2,2)),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    #strip.background = element_rect(colour = "black", fill = "white"),
    panel.border = element_rect(linetype = "solid", fill = NA),
    axis.line = element_line(colour = "black"),
    panel.spacing.x = unit(20, "mm"))

ggsave("prelimresults/YosemiteCarn_NMDSEllipses_240827.pdf", width=10,height=8, units = 'in', dpi= 300)


library(vegan)
library(ggplot2)
library(ggConvexHull)
library(viridis)
library(viridisLite)
nmds <- ggplot(nmds_hulljoin,aes(x = NMDS1, y = NMDS2, col = Species)) +
  geom_convexhull(alpha = 0.25,aes(fill = Species)) + 
  #facet_wrap(~family)+
  geom_point() +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(#strip.background = element_blank(),
        axis.text = element_text(size=20),
        #axis.text.y = element_text(angle = 90, size = 18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        strip.text.x=element_text(size=20, margin = margin(2,2,2,2)),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(20, "mm"))
nmds
nmdsspecies <- ggplot(nmds_hulljoin,aes(x = NMDS1, y = NMDS2, col = Species)) +
  geom_convexhull(alpha = 0.5,aes(fill = Species)) + 
  facet_wrap(~Species, ncol = 5)+
  geom_point() +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        #axis.text.y = element_text(angle = 90, size = 18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #legend.title = element_text(size = 16),
        #legend.text = element_text(size = 14),
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(10, "mm"))
nmdsspecies
ggsave("prelimresults/YosemiteCarn_NMDSNicheSpace_Separate_240418.png", width=16,height=4.5, units = 'in', dpi= 600)

library(ggpubr)
### save wpoo and niche figure 
ggsave(ggarrange(groupwpoo,
                 nmds_ellipse,
                 widths = c(1.4,1.3),
                 nrow=1, ncol = 2), 
       filename=paste("prelimresults/YosemiteDiet_wPOO_nicheellipse", Sys.Date(), ".jpg", sep=""), 
       height=6, width=17, units="in", dpi=300)

### save wpoo and niche figure 
ggsave(ggarrange(groupwpoo,
                 nmds_ellipse,
                 widths = c(1.4,1.3),
                 nrow=1, ncol = 2), 
       filename=paste("Manuscript/ManuscriptFigures/YosemiteDiet_wPOO_nicheellipse", Sys.Date(), ".pdf", sep=""), 
       height=6, width=17, units="in", dpi=300)

#### reshape data so you have focc of prey items per samples per species ####
library(sjmisc)
Subset01b <- Subset01 %>% column_to_rownames(var="PreyTaxa")%>% rotate_df()
Subset01b$Species <- rep("Marten",nrow(Subset01b))
Subset02b <- Subset02 %>% column_to_rownames(var="PreyTaxa")%>% rotate_df()
Subset02b$Species <- rep("Coyote",nrow(Subset02b))
Subset03b <- Subset03 %>% column_to_rownames(var="PreyTaxa")%>% rotate_df()
Subset03b$Species <- rep("Bobcat",nrow(Subset03b))
Subset04b <- Subset04 %>% column_to_rownames(var="PreyTaxa")%>% rotate_df()
Subset04b$Species <- rep("Grayfox",nrow(Subset04b))
Subset05b <- Subset05 %>% column_to_rownames(var="PreyTaxa")%>% rotate_df()
Subset05b$Species <- rep("Cougar",nrow(Subset05b))
list_of_dataframes <- list(Subset01b,
                           Subset02b,
                           Subset03b,
                           Subset04b,
                           Subset05b)
allsamps <- bind_rows(list_of_dataframes, .id = "column_label") %>% 
  replace(is.na(.), 0) %>% 
  rownames_to_column(.,"Sample") %>% 
  mutate(Species = c(rep("Marten",nrow(Subset01b)),
                     rep("Coyote",nrow(Subset02b)),
                     rep("Bobcat",nrow(Subset03b)),
                     rep("Grayfox",nrow(Subset04b)),
                     rep("Cougar",nrow(Subset05b)))) %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

martshandat <- Subset01 %>% 
  column_to_rownames(var = "PreyTaxa")
coyshandat <- Subset02 %>% 
  column_to_rownames(var = "PreyTaxa")
bobshandat <- Subset03 %>% 
  column_to_rownames(var = "PreyTaxa")
greyshandat <- Subset04 %>% 
  column_to_rownames(var = "PreyTaxa")
cougshandat <- Subset05 %>% 
  column_to_rownames(var = "PreyTaxa")

martshan <- diversity(martshandat)
mean(martshan) # [1] 0.82
coyshan <- diversity(coyshandat)
mean(coyshan) # [1] 1.25
bobshan <- diversity(bobshandat)
mean(bobshan) # [1] 0.76
greyshan <- diversity(greyshandat)
mean(greyshan) # [1] 0.65
cougshan <- diversity(cougshandat)
mean(cougshan) # [1] 0.72

martrichdat <- Subset01 %>% 
  column_to_rownames(var = "PreyTaxa")%>% 
  mutate_if(is.numeric, ~1 * (. > 0))
coyrichdat <- Subset02 %>% 
  column_to_rownames(var = "PreyTaxa")%>% 
  mutate_if(is.numeric, ~1 * (. > 0))
bobrichdat <- Subset03 %>% 
  column_to_rownames(var = "PreyTaxa")%>% 
  mutate_if(is.numeric, ~1 * (. > 0))
greyrichdat <- Subset04 %>% 
  column_to_rownames(var = "PreyTaxa")%>% 
  mutate_if(is.numeric, ~1 * (. > 0))
cougrichdat <- Subset05 %>% 
  column_to_rownames(var = "PreyTaxa")%>% 
  mutate_if(is.numeric, ~1 * (. > 0))

martrich <- colSums(martrichdat)
mean(martrich) # [1] 1.30
coyrich <- colSums(coyrichdat)
mean(coyrich) # [1] 1.45
bobrich <- colSums(bobrichdat)
mean(bobrich) # [1] 1.04
greyrich <- colSums(greyrichdat)
mean(greyrich) # [1]1.56
cougrich <- colSums(cougrichdat)
mean(cougrich) # [1] 1.64

### summary statistics ###
# % sctas w/ mammals
classsum <- sealD %>% group_by(class) %>% mutate(count = length(class))
#### taxa summaries ###
# class
classsum <- sealD %>% 
  dplyr::group_by(labid.y,
                  class) %>% 
  dplyr::summarize(count = nrow(labid.y)) %>% ungroup
# order
ordsum <- sealD %>% 
  dplyr::group_by(labid.y,
                  order) %>% 
  dplyr::summarize(count = nrow(labid.y)) %>% ungroup
# family
famsum <- sealD %>% 
  dplyr::group_by(labid.y,
                  family) %>% 
  dplyr::summarize(count = nrow(labid.y)) %>% ungroup

# genus
gensum <- sealD %>% 
  dplyr::group_by(labid.y,
                  genus) %>% 
  dplyr::summarize(count = nrow(labid.y)) %>% ungroup

length(unique(classsum$labid.y))
# 1079
table(classsum$class)
#Actinopterygii       Amphibia           Aves       Mammalia       Reptilia 
#71                     23            144           1050              1
# fish   71/1079*100
# 8.110 [1] 6.580167
# amphib  23/1079*100
# 2.11981 [1] 2.131603
# birds  114/1079*100
# 13.271 [1] 10.56534
# mammals  1050/1079*100
# 96.77419 [1] 97.31233
# reptiles 1/1079*100
# 0.0921659
length(unique(ordsum$order)) # 16
length(unique(famsum$family)) # 36
length(unique(gensum$genus)) # 57


#################### prey/scat content relationships to landcover and structure ###########
library(dplyr)
library(rjags)
library(abind)
library(jagsUI)
#library(rgeos)
#library(sp)
library(reshape2)
#library(rgdal)
#library(raster)
library(adehabitatHR)
library(readxl)
library(sf)
library(sfheaders)
library(terra)

# set working directory
setwd("~/YosemiteDiet/OccupancyAnalysis")
scat <- read.csv("~/YosemiteDiet/YOSE_ResultsAll_12sDietMetabarcoding_scatmetadata.csv", header = T)
names(scat) <- tolower(names(scat))

# read in occupancy and detection covariates
# read in occupancy and detection covariates
elev <- terra::rast('./YosemiteDietOccupancyCovariates/elevation_YNPclip.tif')
snpck <- terra::rast('./YosemiteDietOccupancyCovariates/snpck1000m.tif')
snpck2019 <- snpck[[9]]
snpck2020 <- snpck[[10]]
snpckstack <- c(snpck2019,
                snpck2020)
snpckmean <- terra::mean(snpckstack, na.rm = TRUE)
cancov <- terra::rast('./YosemiteDietOccupancyCovariates/cc1000m.tif')
cancov2019 <- cancov[[9]]
cancov2020 <- cancov[[10]]
cancovstack <- c(cancov2019,
                   cancov2020)
cancovmean <- terra::mean(cancovstack, na.rm = TRUE)
sdcancov <- terra::rast('./YosemiteDietOccupancyCovariates/sdcc1000m.tif')
sdcancov2019 <- sdcancov[[9]]
sdcancov2020 <- sdcancov[[10]]
sdcancovstack <- c(sdcancov2019,
                   sdcancov2020)
sdcancovmean <- terra::mean(sdcancovstack, na.rm = TRUE)
crs <-"+init=EPSG:3310"
elev <- terra::project(elev, crs)

# scat locations
scats <- st_as_sf(scat, coords = c("utmx","utmy"), crs = 26911) %>% sf::st_transform(., crs = 3310) %>% filter(speciesid_osu == "Canis latrans"|
                                                                                      speciesid_osu == "Martes caurina"|
                                                                                      speciesid_osu == "Lynx rufus"|
                                                                                      speciesid_osu == "Puma concolor"|
                                                                                      speciesid_osu == "Urocyon cinereoargenteus")

# buffer marten scats by 1km
scats_buffer <- st_buffer(scats, dist = 1000)
scats_buffer$ID <- c(1:nrow(scats_buffer))
scats_buffer$area <- st_area(scats_buffer)
scats_buffer$areakm <- units::set_units(scats_buffer$area, km^2)

####
#buff.vect<-terra::vect(scats_buffer) # terra extract works with vector not sf object
scats_buffer$elev <- terra::extract(elev, scats_buffer, fun = mean)
scats_buffer$cancov <- terra::extract(cancovmean, scats_buffer, fun = mean)
scats_buffer$sdcancov <- terra::extract(sdcancovmean, scats_buffer, fun = mean)
scats_buffer$snow <- terra::extract(snpckmean, scats_buffer, fun = mean)
scats_buffer$elev <- scats_buffer$elev$elevation
scats_buffer$cancov <- scats_buffer$cancov$mean
scats_buffer$sdcancov <- scats_buffer$sdcancov$mean
scats_buffer$snow <- scats_buffer$snow$mean

##### manyglm to look at effects of habitat on probability of consuming certain prey ########
#load packages
require(mvabund)
require(dplyr)
require(reshape2)

###############
#load packages
###############
require(rgdal)
require(raster)
require(sf)
require(landscapemetrics)
require(tidyr)
require(mvabund)
require(ggplot2)
require(ggpubr)
library(tidyr)
#summarize number of verts per scat
# set working directory
setwd("~/YosemiteDiet")
diet <- read.csv("ProcessedYosemite_PhyloSummaryData_240417.csv")
sizecats <- read.csv("YOSE_SpeciesSummary_SizeTaxonomy_230703.csv")
diet <- left_join(diet, sizecats, by = c("speciessimpleidgroomed.y" = "speciessimpleidgroomed"))
dietdata <- diet[!duplicated(diet[,c("taxasizegroup","labid.y")]),]
dietdata <- dietdata %>% filter(speciesid_osu == "Canis latrans")
dietdata <- spread(dietdata[,c("labid.y","taxasizegroup","numreads")], taxasizegroup, numreads, fill=0) #dims: 202 rows, 55 species, should be 128 rows
dietdata[is.na(dietdata)] <- 0 #replace NAs with 0
dietdata[,2:9][dietdata[,2:9] > 0] <- 1

joined <- left_join(scats_buffer, dietdata, 
                    by = join_by("labid" == "labid.y")) %>% dplyr::filter(speciesid_osu == "Canis latrans")
joined <- joined %>% filter(!is.na(BirdsSmall))%>% dplyr::filter(!is.nan(cancov)) 
joined <- as.data.frame(joined)
data <- joined %>% dplyr::select(BirdsMedium, BirdsSmall, MammalsLarge, MammalsMedium, MammalsSmall)
varsall <- joined %>% dplyr::select(snow, cancov, sdcancov) %>% dplyr::filter(!is.nan(cancov))
cor(varsall)
vars <- joined %>% dplyr::select(cancov)
###############
set.seed(1234)
data <- mvabund(data)
#1000 m buffer
mod.sl_for <- manyglm(formula = data ~ cancov, data=vars, family="binomial")
plot(mod.sl_for)
anova(mod.sl_for)
anova(mod.sl_for, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(cancov = rep(seq(0,54, by=0.55)))

pred <- predict(mod.sl_for, newdata, interval="prediction")
conf <- predict(mod.sl_for, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("cancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
vert.colors <- hcl.colors(9, palette = "tofino")
vert.colors
"#D6E0FF" "#99A5E0" "#666E9A" "#373D58" "#111111" "#2D4424" "#557A47" "#82B56F" "#C2EFB4"
coyotecancov <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=cancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=cancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=cancov)) +
  xlab("Canopy cover") + ylab("") + 
  ylim(0,1)+
  xlim(0, 55)+
  facet_wrap(~taxa,ncol = 6) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
coyotecancov
ggsave(coyotecancov, filename="CoyoteDiet_PreyProbs_NDVI.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(snow)
###############
mod.sl_shrub <- manyglm(formula = data ~ snow, data=vars, family="binomial")
plot(mod.sl_shrub)
anova(mod.sl_shrub)
anova(mod.sl_shrub, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(snow = rep(seq(138,840, by = 8)))

pred <- predict(mod.sl_shrub, newdata, interval="prediction")
conf <- predict(mod.sl_shrub, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("snow","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

coyotesnowplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=snow, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=snow, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=snow)) +
  xlab("Maximum snowpack (m)") + ylab("") + 
  ylim(0,1)+
  xlim(130,840)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
coyotesnowplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_percentshrub.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(sdcancov)
###############
data <- mvabund(data)
#1000 m buffer
mod.sl_edge <- manyglm(formula = data ~ sdcancov, data=vars, family="binomial")
plot(mod.sl_edge)
anova(mod.sl_edge)
anova(mod.sl_edge, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(sdcancov = rep(seq(0,20, by = 0.2)))

pred <- predict(mod.sl_edge, newdata, interval="prediction")
conf <- predict(mod.sl_edge, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("sdcancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
coyotesdcancovplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=sdcancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=sdcancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=sdcancov)) +
  xlab(bquote("Canopy cover (SD)")) + ylab("Presence in coyote scat\n") + 
  ylim(0,1)+
  xlim(0,9)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
coyotesdcancovplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_edgeden.tiff", height=6, width=12, units="in", dpi=300, compression="lzw")

yotediet <- ggsave(ggarrange(coyotecancov, coyotesdcancovplot, coyotesnowplot, 
                             ncol = 1,
                             nrow = 3), 
                   filename=paste("./prelimresults/CoyoteDiet_PreyProbs_covresponse", Sys.Date(), ".jpeg", sep=""), 
                   height=12, width=20, units="in", dpi=600)


##### bobcat #####
dietdata <- diet[!duplicated(diet[,c("taxasizegroup","labid.y")]),]
dietdata <- dietdata %>% filter(speciesid_osu == "Lynx rufus")
dietdata <- spread(dietdata[,c("labid.y","taxasizegroup","numreads")], taxasizegroup, numreads, fill=0) #dims: 202 rows, 55 species, should be 128 rows
dietdata[is.na(dietdata)] <- 0 #replace NAs with 0
dietdata[,2:6][dietdata[,2:6] > 0] <- 1

joined <- left_join(scats_buffer, dietdata, 
                    by = join_by("labid" == "labid.y")) %>% dplyr::filter(speciesid_osu == "Lynx rufus")
joined <- joined %>% filter(!is.na(BirdsSmall))%>% dplyr::filter(!is.nan(cancov)) 
joined <- as.data.frame(joined)
data <- joined %>% dplyr::select(BirdsMedium, BirdsSmall, MammalsLarge, MammalsMedium, MammalsSmall)
varsall <- joined %>% dplyr::select(snow, cancov, sdcancov) %>% dplyr::filter(!is.nan(cancov))
cor(varsall)
vars <- joined %>% dplyr::select(cancov)
###############
set.seed(1234)
data <- mvabund(data)
#1000 m buffer
mod.sl_for <- manyglm(formula = data ~ cancov, data=vars, family="binomial")
plot(mod.sl_for)
anova(mod.sl_for)
anova(mod.sl_for, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(cancov = rep(seq(0,54, by=0.55)))

pred <- predict(mod.sl_for, newdata, interval="prediction")
conf <- predict(mod.sl_for, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("cancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
vert.colors <- hcl.colors(9, palette = "tofino")
vert.colors
"#D6E0FF" "#99A5E0" "#666E9A" "#373D58" "#111111" "#2D4424" "#557A47" "#82B56F" "#C2EFB4"
bobcatcancov <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=cancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=cancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=cancov)) +
  xlab("Canopy cover") + ylab("") + 
  ylim(0,1)+
  xlim(0, 55)+
  facet_wrap(~taxa,ncol = 6) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
bobcatcancov
ggsave(coyotecancov, filename="CoyoteDiet_PreyProbs_NDVI.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(snow)
###############
mod.sl_shrub <- manyglm(formula = data ~ snow, data=vars, family="binomial")
plot(mod.sl_shrub)
anova(mod.sl_shrub)
anova(mod.sl_shrub, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(snow = rep(seq(138,840, by = 8)))

pred <- predict(mod.sl_shrub, newdata, interval="prediction")
conf <- predict(mod.sl_shrub, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("snow","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

bobcatsnowplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=snow, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=snow, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=snow)) +
  xlab("Maximum snowpack (m)") + ylab("") + 
  ylim(0,1)+
  xlim(130,840)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
bobcatsnowplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_percentshrub.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(sdcancov)
###############
data <- mvabund(data)
#1000 m buffer
mod.sl_edge <- manyglm(formula = data ~ sdcancov, data=vars, family="binomial")
plot(mod.sl_edge)
anova(mod.sl_edge)
anova(mod.sl_edge, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(sdcancov = rep(seq(0,20, by = 0.2)))

pred <- predict(mod.sl_edge, newdata, interval="prediction")
conf <- predict(mod.sl_edge, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("sdcancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
bobcatsdcancovplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=sdcancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=sdcancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=sdcancov)) +
  xlab(bquote("Canopy cover (SD)")) + ylab("Presence in bobcat scat\n") + 
  ylim(0,1)+
  xlim(0,9)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
bobcatsdcancovplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_edgeden.tiff", height=6, width=12, units="in", dpi=300, compression="lzw")

bobdiet <- ggsave(ggarrange(bobcatcancov, bobcatsdcancovplot, bobcatsnowplot, 
                             ncol = 1,
                             nrow = 3), 
                   filename=paste("./prelimresults/BobcatDiet_PreyProbs_covresponse", Sys.Date(), ".jpeg", sep=""), 
                   height=12, width=20, units="in", dpi=600)

##### marten #####
dietdata <- diet[!duplicated(diet[,c("taxasizegroup","labid.y")]),]
dietdata <- dietdata %>% filter(speciesid_osu == "Martes caurina")
dietdata <- spread(dietdata[,c("labid.y","taxasizegroup","numreads")], taxasizegroup, numreads, fill=0) #dims: 202 rows, 55 species, should be 128 rows
dietdata[is.na(dietdata)] <- 0 #replace NAs with 0
dietdata[,2:8][dietdata[,2:8] > 0] <- 1

joined <- left_join(scats_buffer, dietdata, 
                    by = join_by("labid" == "labid.y")) %>% dplyr::filter(speciesid_osu == "Martes caurina")
joined <- joined %>% filter(!is.na(BirdsSmall))%>% dplyr::filter(!is.nan(cancov)) 
joined <- as.data.frame(joined)
data <- joined %>% dplyr::select(BirdsMedium, BirdsSmall, MammalsLarge, MammalsMedium, MammalsSmall)
varsall <- joined %>% dplyr::select(snow, cancov, sdcancov) %>% dplyr::filter(!is.nan(cancov))
cor(varsall)
vars <- joined %>% dplyr::select(cancov)
###############
set.seed(1234)
data <- mvabund(data)
#1000 m buffer
mod.sl_for <- manyglm(formula = data ~ cancov, data=vars, family="binomial")
plot(mod.sl_for)
anova(mod.sl_for)
anova(mod.sl_for, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(cancov = rep(seq(0,54, by=0.55)))

pred <- predict(mod.sl_for, newdata, interval="prediction")
conf <- predict(mod.sl_for, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("cancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
vert.colors <- hcl.colors(9, palette = "tofino")
vert.colors
"#D6E0FF" "#99A5E0" "#666E9A" "#373D58" "#111111" "#2D4424" "#557A47" "#82B56F" "#C2EFB4"
martencancov <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=cancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=cancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=cancov)) +
  xlab("Canopy cover") + ylab("") + 
  ylim(0,1)+
  xlim(0, 55)+
  facet_wrap(~taxa,ncol = 6) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
martencancov
ggsave(coyotecancov, filename="CoyoteDiet_PreyProbs_NDVI.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(snow)
###############
mod.sl_shrub <- manyglm(formula = data ~ snow, data=vars, family="binomial")
plot(mod.sl_shrub)
anova(mod.sl_shrub)
anova(mod.sl_shrub, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(snow = rep(seq(138,840, by = 8)))

pred <- predict(mod.sl_shrub, newdata, interval="prediction")
conf <- predict(mod.sl_shrub, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("snow","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

martensnowplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=snow, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=snow, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=snow)) +
  xlab("Maximum snowpack (m)") + ylab("") + 
  ylim(0,1)+
  xlim(130,840)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
martensnowplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_percentshrub.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(sdcancov)
###############
data <- mvabund(data)
#1000 m buffer
mod.sl_edge <- manyglm(formula = data ~ sdcancov, data=vars, family="binomial")
plot(mod.sl_edge)
anova(mod.sl_edge)
anova(mod.sl_edge, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(sdcancov = rep(seq(0,20, by = 0.2)))

pred <- predict(mod.sl_edge, newdata, interval="prediction")
conf <- predict(mod.sl_edge, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("sdcancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
martensdcancovplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=sdcancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=sdcancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=sdcancov)) +
  xlab(bquote("Canopy cover (SD)")) + ylab("Presence in marten scat\n") + 
  ylim(0,1)+
  xlim(0,9)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
martensdcancovplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_edgeden.tiff", height=6, width=12, units="in", dpi=300, compression="lzw")

martendiet <- ggsave(ggarrange(martencancov, martensdcancovplot, martensnowplot, 
                            ncol = 1,
                            nrow = 3), 
                  filename=paste("./prelimresults/MartenDiet_PreyProbs_covresponse", Sys.Date(), ".jpeg", sep=""), 
                  height=12, width=20, units="in", dpi=600)


##### cougar #####
dietdata <- diet[!duplicated(diet[,c("taxasizegroup","labid.y")]),]
dietdata <- dietdata %>% filter(speciesid_osu == "Puma concolor")
dietdata <- spread(dietdata[,c("labid.y","taxasizegroup","numreads")], taxasizegroup, numreads, fill=0) #dims: 202 rows, 55 species, should be 128 rows
dietdata[is.na(dietdata)] <- 0 #replace NAs with 0
dietdata[,2:5][dietdata[,2:5] > 0] <- 1

joined <- left_join(scats_buffer, dietdata, 
                    by = join_by("labid" == "labid.y")) %>% dplyr::filter(speciesid_osu == "Puma concolor")
joined <- joined %>% filter(!is.na(BirdsMedium))%>% dplyr::filter(!is.nan(cancov)) 
joined <- as.data.frame(joined)
data <- joined %>% dplyr::select(BirdsMedium, MammalsLarge, MammalsMedium, MammalsSmall)
varsall <- joined %>% dplyr::select(snow, cancov, sdcancov) %>% dplyr::filter(!is.nan(cancov))
cor(varsall)
vars <- joined %>% dplyr::select(cancov)
###############
set.seed(1234)
data <- mvabund(data)
#1000 m buffer
mod.sl_for <- manyglm(formula = data ~ cancov, data=vars, family="binomial")
plot(mod.sl_for)
anova(mod.sl_for)
anova(mod.sl_for, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(cancov = rep(seq(0,54, by=0.55)))

pred <- predict(mod.sl_for, newdata, interval="prediction")
conf <- predict(mod.sl_for, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("cancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
vert.colors <- hcl.colors(9, palette = "tofino")
vert.colors
"#D6E0FF" "#99A5E0" "#666E9A" "#373D58" "#111111" "#2D4424" "#557A47" "#82B56F" "#C2EFB4"
cougcancov <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=cancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=cancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=cancov)) +
  xlab("Canopy cover") + ylab("") + 
  ylim(0,1)+
  xlim(0, 55)+
  facet_wrap(~taxa,ncol = 6) + 
  scale_color_manual(values = c("#2D4424", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
cougcancov
ggsave(coyotecancov, filename="CoyoteDiet_PreyProbs_NDVI.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(snow)
###############
mod.sl_shrub <- manyglm(formula = data ~ snow, data=vars, family="binomial")
plot(mod.sl_shrub)
anova(mod.sl_shrub)
anova(mod.sl_shrub, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(snow = rep(seq(138,840, by = 8)))

pred <- predict(mod.sl_shrub, newdata, interval="prediction")
conf <- predict(mod.sl_shrub, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("snow","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

cougsnowplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=snow, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=snow, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=snow)) +
  xlab("Maximum snowpack (m)") + ylab("") + 
  ylim(0,1)+
  xlim(130,840)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424","#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
cougsnowplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_percentshrub.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(sdcancov)
###############
data <- mvabund(data)
#1000 m buffer
mod.sl_edge <- manyglm(formula = data ~ sdcancov, data=vars, family="binomial")
plot(mod.sl_edge)
anova(mod.sl_edge)
anova(mod.sl_edge, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(sdcancov = rep(seq(0,20, by = 0.2)))

pred <- predict(mod.sl_edge, newdata, interval="prediction")
conf <- predict(mod.sl_edge, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("sdcancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
cougsdcancovplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=sdcancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=sdcancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=sdcancov)) +
  xlab(bquote("Canopy cover (SD)")) + ylab("Presence in cougar scat\n") + 
  ylim(0,1)+
  xlim(0,9)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
cougsdcancovplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_edgeden.tiff", height=6, width=12, units="in", dpi=300, compression="lzw")

cougdiet <- ggsave(ggarrange(cougcancov, cougsdcancovplot, cougsnowplot, 
                               ncol = 1,
                               nrow = 3), 
                     filename=paste("./prelimresults/CougarDiet_PreyProbs_covresponse", Sys.Date(), ".jpeg", sep=""), 
                     height=12, width=20, units="in", dpi=600)

##### gray fox #####
dietdata <- diet[!duplicated(diet[,c("taxasizegroup","labid.y")]),]
dietdata <- dietdata %>% filter(speciesid_osu == "Urocyon cinereoargenteus")
dietdata <- spread(dietdata[,c("labid.y","taxasizegroup","numreads")], taxasizegroup, numreads, fill=0) #dims: 202 rows, 55 species, should be 128 rows
dietdata[is.na(dietdata)] <- 0 #replace NAs with 0
dietdata[,2:8][dietdata[,2:8] > 0] <- 1

joined <- left_join(scats_buffer, dietdata, 
                    by = join_by("labid" == "labid.y")) %>% dplyr::filter(speciesid_osu == "Urocyon cinereoargenteus")
joined <- joined %>% filter(!is.na(BirdsSmall))%>% dplyr::filter(!is.nan(cancov)) 
joined <- as.data.frame(joined)
data <- joined %>% dplyr::select(BirdsMedium, BirdsSmall, MammalsLarge, MammalsMedium, MammalsSmall)
varsall <- joined %>% dplyr::select(snow, cancov, sdcancov) %>% dplyr::filter(!is.nan(cancov))
cor(varsall)
vars <- joined %>% dplyr::select(cancov)
###############
set.seed(1234)
data <- mvabund(data)
#1000 m buffer
mod.sl_for <- manyglm(formula = data ~ cancov, data=vars, family="binomial")
plot(mod.sl_for)
anova(mod.sl_for)
anova(mod.sl_for, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(cancov = rep(seq(0,54, by=0.55)))

pred <- predict(mod.sl_for, newdata, interval="prediction")
conf <- predict(mod.sl_for, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("cancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
vert.colors <- hcl.colors(9, palette = "tofino")
vert.colors
"#D6E0FF" "#99A5E0" "#666E9A" "#373D58" "#111111" "#2D4424" "#557A47" "#82B56F" "#C2EFB4"
foxcancov <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=cancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=cancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=cancov)) +
  xlab("Canopy cover") + ylab("") + 
  ylim(0,1)+
  xlim(0, 55)+
  facet_wrap(~taxa,ncol = 6) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
foxcancov
ggsave(coyotecancov, filename="CoyoteDiet_PreyProbs_NDVI.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(snow)
###############
mod.sl_shrub <- manyglm(formula = data ~ snow, data=vars, family="binomial")
plot(mod.sl_shrub)
anova(mod.sl_shrub)
anova(mod.sl_shrub, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(snow = rep(seq(138,840, by = 8)))

pred <- predict(mod.sl_shrub, newdata, interval="prediction")
conf <- predict(mod.sl_shrub, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("snow","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

foxsnowplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=snow, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=snow, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=snow)) +
  xlab("Maximum snowpack (m)") + ylab("") + 
  ylim(0,1)+
  xlim(130,840)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
foxsnowplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_percentshrub.tiff", height=4, width=15, units="in", dpi=300, compression="lzw")

######
# elevation 
vars <- joined %>% dplyr::select(sdcancov)
###############
data <- mvabund(data)
#1000 m buffer
mod.sl_edge <- manyglm(formula = data ~ sdcancov, data=vars, family="binomial")
plot(mod.sl_edge)
anova(mod.sl_edge)
anova(mod.sl_edge, p.uni="adjusted")

#plot predictive plots
summary(vars)
newdata <- data.frame(sdcancov = rep(seq(0,20, by = 0.2)))

pred <- predict(mod.sl_edge, newdata, interval="prediction")
conf <- predict(mod.sl_edge, newdata, se.fit=T)
conf <- conf$se.fit
names(conf) <- names(pred)

wide.pred <- cbind(newdata, pred)
#convert to long format so can plot all taxa together
long.pred <- gather(wide.pred, taxa, pred, BirdsMedium:MammalsSmall, factor_key=T) #factor_key = treat new key column as a factor
long.pred$real <- exp(long.pred$pred)/(1+exp(long.pred$pred)) #undo logit link

wide.conf <- cbind(newdata, conf)
names(wide.conf) <- names(wide.pred)
long.conf <- gather(wide.conf, taxa, conf, BirdsMedium:MammalsSmall, factor_key=T)
long.conf <- merge(long.conf, long.pred, by=c("sdcancov","taxa"))
long.conf$low <- exp(long.conf$pred - long.conf$conf*1.96)/(1+exp(long.conf$pred - long.conf$conf*1.96))
long.conf$high <- exp(long.conf$pred + long.conf$conf*1.96)/(1+exp(long.conf$pred + long.conf$conf*1.96))

# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Insecta, col=scatseason), lwd=1) + theme_bw(base_size=20)
# ggplot(data=cbind(newdata, pred)) + geom_line(aes(x=percentage_inside.r1000, y=Amphibia, col=scatseason), lwd=1) + theme_bw(base_size=20)

long.conf$taxa <- factor(long.conf$taxa, levels=c("BirdsMedium","BirdsSmall", "MammalsLarge", "MammalsMedium", "MammalsSmall"))

library(viridis)
foxsdcancovplot <- ggplot(data=long.conf) + 
  geom_ribbon(aes(x=sdcancov, ymin=low, ymax=high, col = taxa, fill = taxa), alpha=0.25) +
  geom_line(aes(x=sdcancov, y=real, col = taxa, fill = taxa), lwd=1, alpha=0.75) + 
  geom_rug(data=vars, aes(x=sdcancov)) +
  xlab("Canopy cover (SD)") + ylab("Presence in gray fox scat\n") + 
  ylim(0,1)+
  xlim(0,9)+
  facet_wrap(~taxa,ncol = 5) + 
  scale_color_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  scale_fill_manual(values = c("#2D4424", "#557A47", "#99A5E0","#666E9A","#373D58"))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(3,3,3,3)),
        legend.position = "none",
        #strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(8, "mm"))
foxsdcancovplot
ggsave(p, filename="Figures/MartenDiet_PreyProbs_edgeden.tiff", height=6, width=12, units="in", dpi=300, compression="lzw")

foxdiet <- ggsave(ggarrange(foxcancov, foxsdcancovplot, foxsnowplot, 
                               ncol = 1,
                               nrow = 3), 
                     filename=paste("./prelimresults/GreyFoxDiet_PreyProbs_covresponse", Sys.Date(), ".jpeg", sep=""), 
                     height=12, width=20, units="in", dpi=600)


### chord diagram ####
library(circlize)
library(dplyr)
setwd("~/YosemiteDiet/")
data <- read.csv("YosemiteDiet_AllspeciesSummaryFOOwPOO_240823.csv")
df <- data.frame(from = data$PreyGroup,
                 to = data$Defecator,
                 value = data$wPOO) 

df <- df %>% 
  dplyr::group_by(from, to) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::arrange(desc(from)) %>% 
  dplyr::filter(!from == "Reptiles")

chordDiagram(df)

vert.colors <- hcl.colors(9, palette = "tofino")
carn.colors <- hcl.colors(5, palette = "viridis")
"#4B0055" "#00588B" "#009B95" "#53CC67" "#FDE333"
"#D6E0FF" "#99A5E0" "#666E9A" "#373D58" "#111111" "#2D4424" "#557A47" "#82B56F" "#C2EFB4"
grid.col = c(Bobcat = "#4B0055",
             Cougar= "#00588B",
             Coyote = "#009B95",
             Grayfox = "#53CC67", 
             Marten = "#FDE333",
             Amphibians = "#C2EFB4", 
             BirdsLarge = "#82B56F", 
             BirdsMedium = "#557A47", 
             BirdsSmall = "#2D4424", 
             Fishes = "#111111", 
             MammalsLarge = "#373D58",
             MammalsMedium = "#666E9A",
             MammalsSmall = "#99A5E0")

circos.par(gap.after = c("Bobcat" = 1,
                         "Cougar" = 1,
                         "Coyote" = 1,
                         "Grayfox" = 1, 
                         "Marten" = 1,
                         "Amphibians" = 5, 
                         "Reptiles" = 5, 
                         "BirdsLarge" = 5,
                         "BirdsMedium" = 5, 
                         "BirdsSmall" = 5, 
                         "Fishes" = 5, 
                         "MammalsLarge" = 1,
                         "MammalsMedium" = 1,
                         "MammalsSmall" = 1),
           big.gap = 20)
chordDiagram(df, grid.col = grid.col, big.gap = 30, scale = TRUE, transparency = 0.15,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.01))

circos.clear()
