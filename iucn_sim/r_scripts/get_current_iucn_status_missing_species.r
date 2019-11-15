library(rredlist)

args = commandArgs(trailingOnly = TRUE)
species_list_file = args[1]
iucn_key = args[2]
outdir = args[3]

#species_list_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/output/species_list.txt'
#iucn_key = '01524b67f4972521acd1ded2d8b3858e7fedc7da5fd75b8bb2c5456ea18b01ba'
#outdir = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/output/'

data = read.csv(species_list_file,header = FALSE)
species_list = data$V1
status_list = c()
for (species in species_list){
  #print(species)
  iucn_info = rl_search(species,key = iucn_key)
  category = iucn_info$result$category
  if (is.null(category)){
    category = NaN
  }
  status_list = c(status_list,category)
}

species_cat_df = cbind(as.character(species_list),status_list)
write.table(species_cat_df,file=paste0(outdir,'/current_status_missing_species.txt'),quote=FALSE,sep = '\t',col.names = FALSE,row.names=FALSE)

            