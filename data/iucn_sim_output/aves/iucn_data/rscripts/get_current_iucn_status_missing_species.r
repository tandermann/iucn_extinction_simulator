
    library(rredlist)

    args = commandArgs(trailingOnly = TRUE)
    species_list_file = args[1]
    iucn_key = args[2]
    outdir = args[3]
    
    data = read.csv(species_list_file,header = FALSE)
    species_list = data$V1
    status_list = c()
    for (i in 1:length(species_list)){
      species = as.character(species_list[i])
      print(paste0('Extracting current status for ', species,' (',i,' of ',length(species_list),')'))
      iucn_info = rl_search(species,key = iucn_key)
      category = iucn_info$result$category
      if (is.null(category)){
        category = NaN
      }
      status_list = c(status_list,category)
    }
    
    species_cat_df = cbind(as.character(species_list),status_list)
    write.table(species_cat_df,file=paste0(outdir,'/current_status_missing_species.txt'),quote=FALSE,sep = '	',col.names = FALSE,row.names=FALSE)
    