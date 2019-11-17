library(rredlist)

args = commandArgs(trailingOnly = TRUE)
taxon_group = args[1]
group_rank = args[2]
iucn_key = args[3]
outdir = args[4]

# load all IUCN data
data = c()
for (i in seq(0, 10, 1)){
  data = c(data,c(rl_sp(key=iucn_key,page = i)))
}

# get taxon list, class list and status list from data
taxon_list = c()
group_list = c()
status_list = c()
taxon_id = c()
for (page in data){
  if (length(page) > 1){
    target_column = which(startsWith(colnames(page),group_rank))
    taxon_list = c(taxon_list,page$scientific_name)
    group_list = c(group_list,page[,target_column])
    status_list = c(status_list,page$category)
    taxon_id = c(taxon_id,page$taxonid)
  }
}

# exclude extinct taxa
boolean = !grepl('EX|EW',status_list)
taxon_list = taxon_list[boolean]
group_list = group_list[boolean]
status_list = status_list[boolean]
taxon_id = taxon_id[boolean]

# remove all non-species level identifications
boolean = !grepl('subsp.|ssp.|subpopulation|Subpopulation',taxon_list)
taxon_list = taxon_list[boolean]
group_list = group_list[boolean]
status_list = status_list[boolean]
taxon_id = taxon_id[boolean]

# select mammals
selected_taxon_list = taxon_list[group_list==taxon_group]
selected_ids = taxon_id[group_list==taxon_group]
final_sorted_taxon_list = selected_taxon_list
#final_taxon_list = as.data.frame(cbind(selected_taxon_list,selected_ids))
#final_sorted_taxon_list = final_taxon_list[order(final_taxon_list$selected_taxon_list),]
write.table(final_sorted_taxon_list,file=paste0(outdir,'/',taxon_group,"_species_list.txt"), quote=F,row.names=F,sep='\t',col.names = FALSE)


# get historic data __________________________
species_list = selected_taxon_list
# create new dataframe with species as first column
historic_assessments = species_list
historic_assessments = as.data.frame(historic_assessments)
colnames(historic_assessments) = c('species')
# find historic assessments and fill into dataframe
counter = 1
for (species in species_list){
  print(paste0('Downloading IUCN history: species ',counter, ' of ',length(species_list)))
  #print(species)
  row_id = which(historic_assessments$species == species)
  hist_data <- rl_history(species,key=iucn_key)
  for (year in hist_data$result$year){
    id = which(hist_data$result$year == year)
    #some species have multiple assignmen ts for some years
    if (length(hist_data$result$code[id])>1){
      historic_assessments[row_id,year] <- hist_data$result$code[id][1]
    }
    else{
      historic_assessments[row_id,year] <- hist_data$result$code[id]
    }
  }
  counter = counter+1
}
write.table(historic_assessments,file=paste0(outdir,'/',taxon_group,"_iucn_history.txt"), quote=F,row.names=F,sep='\t')
#___________________________________

