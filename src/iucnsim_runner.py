import pandas as pd

# read input species list
species_list_file = 'data/precompiled/gl_data/aves_gl.txt'
species_list = pd.read_csv(species_list_file,sep='\t',header=None)
# define reference group
reference_group = "Aves"
reference_rank = "class"
iucn_key='01524b67f4972521acd1ded2d8b3858e7fedc7da5fd75b8bb2c5456ea18b01ba'


# get IUCN history of reference group
outdir = 'private/iucn_sim_test/iucn_data'
# get iucn history of reference group
get_iucn_history(reference_group=reference_group,
                           reference_rank=reference_rank,
                           outdir=outdir)

iucn_history_file = iucn_sim.get_iucn_history(   reference_group=reference_group,
                                                 reference_rank=reference_rank,
                                                 iucn_key=iucn_key,
                                                 outdir=outdir)
# get idea of reference group stats
counted_status_transition_events = evaluate_iucn_history(iucn_history_file)

# get most recent status for each taxon in target species list
extant_taxa_current_status = get_most_recent_status_target_species(species_list=species_list,
                                                                   iucn_history_file=iucn_history_file,
                                                                   iucn_key=iucn_key,
                                                                   outdir=outdir)








