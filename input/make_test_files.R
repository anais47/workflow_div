## Create test files for workflow

library(stringr) 
library(dplyr)
library(tibble)
library(tidyr)

# Path to write files

path <- "C:/Users/arey/Documents/POSTDOC/DUBII_cours/module_projet_tuto/workflow_div_v1/input/"

# OTU table ----
load("C:/Users/arey/Documents/POSTDOC/ADNe_UVC_div_taxo/Output/otu_all_12S.RData")

tmp_otu <- otu_all_12S %>%
  # we sum the number of reads of the pcr replicates
  group_by(num_otu, 
           sample_norep, 
           Echantillon,
           Extraction_replicat) %>%
  summarize(count_sum_pcr=sum(count)) %>%
  # we do the mean for the sample with two series of DNA extraction
  group_by(num_otu, 
           Echantillon) %>%
  summarize(count_tmp=mean(count_sum_pcr)) %>%
  spread(num_otu, count_tmp, fill=0) %>%
  mutate_if(is.numeric, round, 0) 

write.csv2(tmp_otu, file=paste(path, "eDNA_test_otu_table.csv", sep=""), row.names = FALSE)

# METADATA table ----

metadata <- read.csv2(file="C:/Users/arey/Documents/POSTDOC/ADNe_UVC_div_taxo/Data/metadata_FISHDNA.csv") 
metadata <- metadata %>% 
  rename(sample_norep=Spygen_code)
# change Periode into September and Aout separately
metadata <- metadata %>%
  mutate(Season=case_when(str_detect(Date.ech.eau, "/09/2018") ~ "September",
                          str_detect(Date.ech.eau, "/06/2018") ~ "June",
                          str_detect(Date.ech.eau, "/08/2018") ~ "August")) %>%
  mutate(year="2018") %>%
  rename(eventDate=Date.ech.eau) %>%
  as.data.frame() %>%
  distinct(Echantillon, Site, Season, Methode, Merge_replicat, eventDate, year) %>%
  filter(Site!="Bloscon Port")


uvc_tc_env <- read.csv(file="C:/Users/arey/Documents/POSTDOC/Données_UVC/NEW_data/TC_siteXenv.csv", sep=";", dec=",")

coord_gps <- uvc_tc_env %>%
  # change names of site
  rename(Site=code_site) %>%
  mutate(Site=gsub("COCH", "COCHONS NOIRS", Site),
         Site=gsub("AST", "ASTAN", Site),
         Site=gsub("CORB", "CORBEAU", Site),
         Site=gsub("FIG", "FIGUIER", Site)) %>%
  # keep only sampled sites with eDNA
  filter(Site=="COCHONS NOIRS" |
           Site=="ASTAN" |
           Site=="CORBEAU"|
           Site=="FIGUIER") %>%
  select(Site, longitude_DD, latitude_DD) %>%
  rename(decimalLongitude=longitude_DD, decimalLatitude=latitude_DD) %>%
  unique() 


tmp_metaD <- metadata %>%
  left_join(coord_gps, by="Site")


write.csv2(tmp_metaD, file=paste(path, "eDNA_test_metadata_table.csv", sep=""), row.names = FALSE)

# TAXO table ----

tmp_scname <- tmp_otu %>%
  select(-Echantillon) %>%
  colnames() %>%
  as.data.frame() %>%
  rename(num_otu=".") %>%
  inner_join(otu_all_12S) %>%
  distinct(num_otu, scname_final_forWORMS) %>%
  rename(user_supplied_name=scname_final_forWORMS) %>%
  rename(OTU_ID=num_otu)


write.csv2(tmp_scname, file=paste(path, "eDNA_test_taxo_table.csv", sep=""), row.names = FALSE)

# TRAITS table ----

traits.tmp <-  read.csv(file="C:/Users/arey/Documents/POSTDOC/Données_UVC/TraitCollectionFishNAtlanticNEPacificContShelf.csv", sep=";", na.strings=c("","NA")) #

#taxon.comm <- taxo.table.class %>% distinct(valid_name) %>% rename(taxon=valid_name)

#sp_com <- intersect(taxon.comm$taxon, unique(traits.tmp$taxon))

traits.tmp.2 <- traits.tmp %>% #filter(taxon %in% sp_com) %>%
  filter(LME==24) %>% # en faisant comme ca on perds des espèces mais pour le cas d'étude on s'en fiche
  distinct(taxon, habitat, feeding.mode, tl, body.shape, fin.shape, spawning.type, fecundity, length.infinity, length.max) 


write.csv2(traits.tmp.2, file=paste(path, "eDNA_test_trait_table.csv", sep=""), row.names = FALSE)

# DIMENSIONS ----

dim(tmp_otu)
dim(tmp_scname)
dim(tmp_metaD)



