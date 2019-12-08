#install.packages("https://cran.r-project.org/src/contrib/Archive/mvtnorm/mvtnorm_1.0-8.tar.gz", repos=NULL)
#install.packages('/Users/spacemule/bin/r_packages/Rphylopars_0.2.9.tar.gz', repos = NULL, type="source")
#install.packages("https://cran.r-project.org/src/contrib/Archive/mvtnorm/mvtnorm_1.0-6.tar.gz", repos=NULL)
#install_github("ericgoolsby/Rphylopars",dependencies = TRUE)
#install.packages('Rphylopars',dependencies = T)
require(Rphylopars)
require(ape)
#help(phylopars)

jetz_trees = read.tree(file='/Users/spacemule/Desktop/Tobias_Andermann/iucn_predictions/data/raw/phylogenies/birds/AllBirdsEricson1.tre')
trait_data = read.table('/Users/spacemule/Desktop/Tobias_Andermann/iucn_predictions/data/processed/trait_data/birds/trait_data_birds_in_phylogeny.txt',sep='\t',header = T)

TRAIT = trait_data
TREES = jetz_trees

#sp_list_tree = TREES$UNTITLED$tip.label
#sp_list_traits = as.character(TRAIT$species)
#intersect(sp_list_tree,sp_list_traits)


# MODEL TESTING_________________________________________________________________
i = 1
models = c('BM','mvOU','OU','lambda','kappa','EB','star')
AIC_score = c()
for (model in models){
  print(model)
  RESULTS=phylopars(trait_data = TRAIT,tree = TREES[[i]],pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model=model)
  AIC_score = c(AIC_score,AIC(RESULTS))
}

pdf('/Users/spacemule/Desktop/Tobias_Andermann/iucn_predictions/data/processed/trait_data/birds/model_tests/AIC_scores_aves.pdf')
plot(factor(models),AIC_score,main='Model test - Phylogenetic imputation for Aves',xlab="Model",ylab='AIC score')
dev.off()
#_______________________________________________________________________________



# run best model on n trees___________________________________________________

# create empty data frames for mean and std values to store all results
sp_list_tree = TREES[[1]]$tip.label
n = 100
gl_mean_df = data.frame(matrix(NA, nrow = length(sp_list_tree), ncol = n))
row.names(gl_mean_df) = names(sp_list_tree)
gl_std_df = data.frame(matrix(NA, nrow = length(sp_list_tree), ncol = n))
row.names(gl_std_df) = names(sp_list_tree)

# run the imputation for n replicates
counter = 0
for (i in sample(1:1000, n)){
  counter = counter+1
  print(counter)
  model='EB'
  RESULTS=phylopars(trait_data = TRAIT,tree = TREES[[i]],pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model=model)
  
  # extract imputed values (mean and variance)
  GL_mean = RESULTS$anc_recon[,'log_gl']
  GL_std = sqrt(RESULTS$anc_var)[,'log_gl']
  
  # extract only values of tip labels (no ancestral states)
  #sp_list_tree = TREES[[i]]$tip.label
  GL_mean_tips = GL_mean[names(GL_mean) %in% sp_list_tree]
  GL_std_tips = GL_std[names(GL_std) %in% sp_list_tree]

  gl_mean_df[,counter] = as.vector(GL_mean_tips)
  gl_std_df[,counter] = as.vector(GL_std_tips)
  
  # write everything to file (overwrite after every iteration)
  write.table(gl_mean_df,'/Users/spacemule/Desktop/Tobias_Andermann/iucn_predictions/data/processed/trait_data/birds/phylogenetic_imputation/birds_gl_mean.txt',sep="\t",row.names = FALSE,col.names=FALSE,quote=FALSE)
  write.table(gl_std_df,'/Users/spacemule/Desktop/Tobias_Andermann/iucn_predictions/data/processed/trait_data/birds/phylogenetic_imputation/birds_gl_std.txt',sep="\t",row.names = FALSE,col.names=FALSE,quote=FALSE)
  write.table(sp_list_tree,'/Users/spacemule/Desktop/Tobias_Andermann/iucn_predictions/data/processed/trait_data/birds/phylogenetic_imputation/birds_species_list.txt',sep="\t",row.names = FALSE,col.names=FALSE,quote=FALSE)
}





# draw 1000 values from uncertainty interval of each GL estimate 
#for (i in 1:1000){
#  x <- rnorm(length(as.vector(GL_mean_tips)), mean=as.vector(GL_mean_tips), sd=as.vector(GL_std_tips))
#}
