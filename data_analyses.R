#libraries 
library(tidyverse)
library(stats)
library(dbplyr)
library(ggplot2)
library(RColorBrewer) #To create color palettes
library(boot) #to calculate SD of diversity index
library(kableExtra) # for fancy tables
library(taxize) #for taxonomical annotation

options(digits = 3, scipen = 999)

### load raw data ### 
#microbial data 
micro1 <- read.csv(file = "data/micro.o1.csv", check.names = FALSE)
micro2 <- read.csv(file = "data/micro.o2.csv", check.names = FALSE)
#adjust colume name 
colnames(micro1)[1] <- "sample"
colnames(micro2)[1] <- "sample"

#resistome data 
res1 <- read.csv(file = "data/resistome1.csv", check.names = FALSE)
res2 <- read.csv(file = "data/resistome2.csv", check.names = FALSE)
#adjust column name
colnames(res1)[1] <- "sample"
colnames(res2)[1] <- "sample"


#### Realative abundance micro ####
# join data 
temp_micro1 <- micro1 %>%
  pivot_longer(cols = -sample, names_to = "microbe", values_to = "abundance")

temp_micro2 <- micro2 %>% 
  pivot_longer(cols = -sample, names_to = "microbe", values_to = "abundance")

#create dataframe 
data_micro <- rbind(temp_micro2, temp_micro1) %>% 
  filter(abundance > 0, microbe != "unknown") %>%
  group_by(sample) %>% 
  mutate(relative_abundance = abundance/sum(abundance)) %>%
  select(-abundance) %>% 
  pivot_wider(names_from = sample, values_from = relative_abundance)

#change NA to 0's 
data_micro[is.na(data_micro)] <- 0

#samples names
samples_names <- colnames(data_micro[,-1])

#define sample's groups 
groups_names <- c("farm_1", "farm_1", "farm_2", "farm_2", "farm_3", "farm_3")

#samples + groups 
samples_groups <- tibble(sample = samples_names, farm = groups_names)

#delete objects 
rm(temp_micro1, temp_micro2, micro1, micro2, samples_names, groups_names)



####  Taxonomical annotation #### 

### get complete taxonomic annotation ###

#Orders names 
microbes_names <- unique(data_micro$microbe)
#tax class. 
reclass_taxa <- classification(microbes_names, db = "ncbi")

micro_tax <- reclass_taxa[!is.na(reclass_taxa)]

#create data frame with taxonomic details 
micro_tax <- tibble(names = names(micro_tax), micro_tax)%>%
  unnest(cols = c(micro_tax)) %>% 
  filter(rank %in% c("superkingdom", "phylum","class","order")) %>% # filter for tax level of interest 
  select(-id) %>% 
  spread(rank, name) %>% 
  mutate(root = "ROOT") %>%
  select(names, root, superkingdom, phylum, class, order) %>% 
  unite(col = "lineage", root:order, sep = ";")

#rename column 
colnames(micro_tax)[1] <- "microbe"


#microbes missing information 
miss_micro <- setdiff(microbes_names, micro_tax$microbe)


# missing linage: added manually 
miss_lineage <- c("ROOT;Bacteria;Firmicutes;Bacilli;Turicibacterales", 
"ROOT;Bacteria;Spirochaetes;Spirochaetia;Sphaerochaetales")

#create a table with missing information 
missing_data <- tibble(microbe = miss_micro, 
                       lineage = miss_lineage)

#add to micro_tax 
micro_tax <- rbind(micro_tax, missing_data)

# add taxa to micro data 
data_micro <- data_micro %>% 
  left_join(micro_tax)

#remove objects 
rm(miss_lineage, miss_micro, missing_data, reclass_taxa)

#save file
save(data_micro, file = "rda/data_micro.rda")
save(micro_tax, file = "rda/micro_tax.rda")
save(microbes_names, file = "rda/microbes_names.rda")
save(samples_groups, file = "rda/samples_groups.rda")




#### Resistome data #### 

#put data together 
temp_res1 <- res1 %>%
  pivot_longer(cols = -sample, names_to = "gene_id", values_to = "abundance")

temp_res2 <- res2 %>% 
  pivot_longer(cols = -sample, names_to = "gene_id", values_to = "abundance")

#data resistome 
data_resistome <- rbind(temp_res2, temp_res1)

#load CARD database 
load(file = "rda/CARD_2021.rda")

# select columns 
CARD_2021 <- CARD_2021 %>% 
  select(gene_id, antibiotic_class, mechanism)

#remove duplicated args in this data base 
CARD_2021 <- CARD_2021 %>% 
  filter(!duplicated(gene_id))

#adjust entries with more than one mechanism and adjust names (remove "antibiotic") 
CARD_2021$mechanism <- CARD_2021$mechanism %>% 
  str_remove_all(";[a-z]+\\s[a-z]+\\s[a-z]+") %>% 
  str_remove_all("antibiotic\\s") %>% 
  str_remove_all("\\sto\\santibiotic")

# joing resistome data and card db 
data_resistome <- data_resistome %>% 
  left_join(CARD_2021) %>% 
  arrange(sample)

# save files 
save(data_resistome, file = "rda/data_resistome.rda")

#delete objects 
rm(temp_res1, temp_res2, res1, res2, CARD_2021)


