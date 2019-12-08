#library(taxize)
library(ape)
#iucn_key = '01524b67f4972521acd1ded2d8b3858e7fedc7da5fd75b8bb2c5456ea18b01ba'

phylacine_phy = read.nexus('/Users/tobias/GitHub/iucn_predictions/data/raw/phylogenies/mammals/phylacine/Complete_phylogeny.nex')
phy_species_list = phylacine_phy$UNTITLED$tip.label

iucn_species_list = as.character(read.csv('/Users/tobias/GitHub/iucn_predictions/data/raw/iucn/species_list_mammals.txt',header = FALSE)$V1)
iucn_species_list = gsub(" ", "_", iucn_species_list)

# since both data sources are based on IUCN taxonomy, there is no need to change taxonomy
matches_phy_iucn = intersect(phy_species_list,iucn_species_list)

# get the fraction of overlap
length(matches_phy_iucn)/length(iucn_species_list)



head(iucn_species_list)


# subset = phy_species_list[1:10]
# #synonyms(subset, db='itis')
# src = 'IUCN Red List of Threatened Species'
# # get id of target database
# subset(gnr_datasources(), title %in% src)
# databse_id = 163
# subset = c('Equus khur')
# a = gnr_resolve(names=subset,data_source_ids = c(databse_id),with_context = TRUE,fields = 'all')
# 
