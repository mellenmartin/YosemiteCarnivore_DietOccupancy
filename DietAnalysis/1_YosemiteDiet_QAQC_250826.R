###################
##=============================
# Diet data - simplifies data and rids few reads
#=============================

#### Scats collected by Rogue Detection Teams (2019-2020)
#### Funded by OSU, YOSE NP, Yosemite Conservancy - Sean Matthews, Sarah Stock, David Green, Ben Sacks, Taal Levi 
#### Modified by Marie Martin
#### Modified from code written by Brent Barry     
###################

library(dplyr)
library(tidyr)
library(ggplot2)
remove(list=ls()) 
gc()
### Step 1, Load from csv

setwd("~/YosemiteDiet/") # your working directory
diet <- read.csv("YOSE_ResultsAll_12sDietMetabarcoding_Levi_mem230703.csv")
sizecats <- read.csv("YOSE_SpeciesSummary_SizeTaxonomy_230703.csv")
scatlocs <- read.csv("ScatsProjected_NAD83Z11.csv")
scatlocs <-  scatlocs %>% dplyr::select(SampleIDFi, POINT_X,POINT_Y)
scatlocs <- scatlocs %>% dplyr::rename(fieldid = "SampleIDFi")
scatmetadata <- read.csv("YOSE_ResultsAll_12sDietMetabarcoding_scatmetadata.csv")
scatmetadata <-  scatmetadata %>% dplyr::select(SampleIDFi, SpeciesIDL, SpeciesID_OSU, Date)
scatmetadata <- scatmetadata %>% dplyr::rename(fieldid = "SampleIDFi")

library(plyr)
names(diet) = tolower(names(diet))
diet <- diet %>% filter(!speciessimpleidgroomed=="Enhydra lutris",
                        !best_match_pid.mean <95)

diet2 <- diet %>% tidyr::pivot_longer(cols=c('repa', 'repb', 'repc'),
                             names_to='rep',
                             values_to='counts') 

diet2$key <- paste0(diet2$labid,"_",diet2$rep,diet2$speciessimpleidgroomed) # our merging key
View(diet2)
diet2$counts[is.na(diet2$counts)] <- 0

diet2 <- diet2 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(speciessimpleidgroomed == defecator,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)


## Calculate the percent reads for a prey item per scat
total_reads <- diet2 %>% group_by(labid) %>% summarize_at(vars(counts), sum) #total reads per scat
spec_reads <- diet2 %>% group_by(labid,speciessimpleidgroomed) %>% summarize_at(vars(counts), sum) #total reads per species per scat
key_reads <- diet2 %>% group_by(labid,rep) %>% summarize_at(vars(counts), sum) # merging key

### defecator ###
#defecator <- spec_reads %>%
#  group_by(labid) %>%
#  slice(which.max(counts))
#

diet2 <- merge(diet2, sizecats, by ="speciessimpleidgroomed")
scatlocs <- merge(scatmetadata, scatlocs, by = "fieldid")
diet2 <- left_join(diet2, scatlocs, by = "fieldid")
#combine the reads dataframes to have a single object with total reads per scat, reads per species per scat, 
# labid, and rep
reads <- merge(total_reads,key_reads, by="labid")
reads <- merge(reads,spec_reads,by="labid")
head(reads)

colnames(reads) <- c("labid","total_reads","rep", "rep_reads", "speciessimpleidgroomed","spec_reads")
reads$p_reads <- reads$spec_reads/reads$total_reads # % reads of a given species relative to the whole scat
reads$key <- paste0(reads$labid,"_",reads$rep,reads$speciessimpleidgroomed) # make our key to remerge

# merge back with the larger dataframe
db <- merge(diet2,reads,by= "key",all.y = F)

write.csv(db,"ProcessedYosemite_DietGeoData_240416.csv", row.names = FALSE)

# d1 now has the number of times 
names(db)
d1 = ddply(db,
           .(geographic_region,
             labid.y,
             fieldid,
             POINT_X,
             POINT_Y,
             SpeciesIDL,
             SpeciesID_OSU,
             defecator, 
             commonname,
             species,
             genus,
             family,
             order,
             class,
             phylum,
             speciessimpleidgroomed.y,
             idlevel,
             genusgroup,
             familygroup,
             ordergroup,
             generalgroup,
             taxasizegroup,
             best_match_pid.mean,
             total_reads,
             rep_reads,
             spec_reads,
             p_reads), summarise,
             replicates = ifelse(rep_reads == 0, 0, 1)) #length here is technically off by a bit 
#because there could several instances of the same species in replicate, but it's close enough for now
head(d1)

# filter by our needs to get our final data frame. we can use this dataframe to make plots
draft.dat <- d1[which(d1$p_reads>=0.1),] # remove spp w/ less 1 percent of reads within a replicate
#draft.dat <- draft.dat[which(draft.dat$replicates>1),] # remove species that appear in less than 2 replicates, again this isn't a perfect way to do this but is close enough for now

write.csv(draft.dat,"ProcessedYosemite_DietGeoDataSimplified_240414.csv", row.names = FALSE)

# Extract by species in 7b
# Create frequency of occurrence and plot 7c

