import os
import pandas as pd
import numpy as np
import csv
import scipy as sp
import scipy.stats
from collections import Counter
import pickle
import math
import sys
import matplotlib.pyplot as plt



def format_iucn_df(statuses_through_time_file,start_year,current_year):
    new_df = pd.DataFrame()
    new_df['species'] = statuses_through_time_file['species']
    for number in range(start_year,current_year+1):
        year = str(number)
        if year in statuses_through_time_file.columns:
            pass
        else:
            statuses_through_time_file[year] = np.nan
        new_df[year] = statuses_through_time_file[year]
    new_df.replace({'NE': np.nan},inplace=True)
    print('Completed: Adding missing columns and selecting target years.')
    return new_df

def exclude_extinct(iucn_dataframe, status_column_name):
    new_df = iucn_dataframe.copy()
    new_df = new_df[(new_df[status_column_name] != 'EW') & (new_df[status_column_name] != 'EX')]
    # a few odd species have come back from extinction, meaning that even though their current status is neither 'EX' or 'EW' they had 'EX' as a previous status
    new_df.replace({'EX': 'NaN'},inplace=True)
    new_df.drop(new_df.columns[[-2, -1]], axis=1, inplace=True)
    return new_df

def count_status_changes(iucn_dataframe,valid_status_dict):
    change_types = {}
    for row in iucn_dataframe.iterrows():
        species_index = row[0]
        species_name = iucn_dataframe.species[species_index]
        current_status = ''
        if species_name in valid_status_dict:
            valid_status = list(valid_status_dict[species_name])
            for i in range(len(valid_status)):        
                # compile all change types and count
                current_status = valid_status[i]
                #this takes care of when we are getting out of the list index for the following steps and when we are at the end of the list anyways
                if i+1 >= len(valid_status):
                    pass
                elif valid_status[i+1] != current_status:
                    next_status = valid_status[i+1]
                    change = "%s->%s" %(current_status,next_status)
                    change_types.setdefault(change,[])
                    change_types[change].append(1)
    for category in change_types:
        change_types[category] = sum(change_types[category])
    return change_types

def get_years_spent_in_each_category(extant_iucn_df,valid_status_dict):
    species_list = extant_iucn_df.species.values
    only_year_columns = extant_iucn_df.iloc[:,1:].copy()
    max_index = only_year_columns.shape[-1]
    years_in_status = [np.diff(np.append(np.where(only_year_columns[extant_iucn_df.species==species].notna().values[0])[0],max_index)) for species in species_list]
    status_master_list = np.array([[item for sublist in [[j]*years_in_status[index_i][index_j] for index_j,j in enumerate(valid_status_dict[i])] for item in sublist]+['NA']*(max_index-sum(years_in_status[index_i])) for index_i,i in enumerate(species_list)])
    status_master_flat_array = np.concatenate(status_master_list)
    statuses, counts = np.unique(status_master_flat_array,return_counts=True)
    years_in_each_category = dict(zip(statuses, counts))
    years_in_each_category.pop('NA')
    return years_in_each_category

def replace_iucn_status_with_int(change_types,sum_years):
    iucn_code = {'LC':0, 'NT':1, 'VU':2, 'EN':3, 'CR':4}
    new_change_types = {}
    new_sum_years = {}
    for element in change_types:
        states = element.split('->')
        new_states = []
        new_states.append(str(iucn_code[states[0]]))
        new_states.append(str(iucn_code[states[1]]))
        new_states = '->'.join(new_states)
        new_change_types.setdefault(new_states,change_types[element])
    for status in sum_years:
        new_status = iucn_code[status]
        new_sum_years.setdefault(new_status,sum_years[status])
    return new_change_types,new_sum_years


def save_obj(obj, file_name):
    with open(file_name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)

def qmatrix(rates, status_specific_p_e):   
    Q_matrix = np.zeros((5,6))
    iucn_code = {'LC':0, 'NT':1, 'VU':2, 'EN':3, 'CR':4}
    for row_id,change_type in enumerate(rates.index.values):
        states = change_type.split('->')
        if not 'DD' in states:
            original_state = iucn_code[states[0]]
            derived_state = iucn_code[states[1]]
            rate = rates[row_id]
            Q_matrix[original_state,derived_state] = rate
            Q_matrix[original_state,5] = status_specific_p_e[original_state]   
    # This fills the diagonal line of the matrix with negative values, so that the line has a sum of 0
    np.fill_diagonal(Q_matrix, -np.sum(Q_matrix,axis=1))
    np.set_printoptions(precision=8)
    return Q_matrix

def treat_dd_species(iucn_dataframe,change_type_dict,all_lc = False):
    new_df = iucn_dataframe.copy()
    #get the change frequencies from DD species to a real category and based on that draw new, valid categories for DD data
    dd_change_vector = np.zeros(5)
    status_index_dict = {'LC':0,'NT':1,'VU':2,'EN':3,'CR':4}
    for key in list(change_type_dict.keys()):
        if key.startswith('DD'):
            old_status,new_status = key.split('->')
            new_status_index = status_index_dict[new_status]
            new_value = change_type_dict[key]
            dd_change_vector[new_status_index] = new_value
    dd_change_probabilities_vector = dd_change_vector/sum(dd_change_vector)
    #extant_iucn_df.replace('DD', list(status_index_dict.keys())[list(status_index_dict.values()).index(random_choice_P(dd_change_vector)[1])])
    dd_coords = np.where(new_df == 'DD')
    new_draws = np.random.choice(['LC','NT','VU','EN','CR'], size=len(dd_coords[0]), replace=True, p=dd_change_probabilities_vector)
    if all_lc:
        new_draws = np.array(['LC']*len(dd_coords[0]))
  #  iucn_dataframe.iloc[dd_coords[0],dd_coords[1]]
    new_df = np.array(new_df)
    new_df[new_df=="DD"] = new_draws
    new_df = pd.DataFrame(new_df)
    new_df.columns = np.array(list(iucn_dataframe.columns))
    return new_df

def extract_valid_statuses(formatted_status_through_time_file):
    taxon_list = list(formatted_status_through_time_file.species.values)
    #valid_status_matrix = [list(line[line.isin(['NE','EX','EW','DD','CR','EN','VU','NT','LC'])].values) for it,line in formatted_status_through_time_file.iterrows()]
    df_array = formatted_status_through_time_file.values
    valid_status_matrix = [list(line[np.isin(line,['DD','EX','CR','EN','VU','NT','LC'])]) for line in df_array]
    # if taxon has no single valid status in IUCN history, model as NE at present
    valid_status_matrix = [['NE'] if len(i) == 0 else i for i in valid_status_matrix]
    valid_status_dict = dict(zip(taxon_list, valid_status_matrix))
    current_status_list = [line[-1] for line in valid_status_matrix]
    most_recent_status_dict = dict(zip(taxon_list, current_status_list))
    return valid_status_dict,most_recent_status_dict,current_status_list,taxon_list

def random_choice_P(vector): # randomly sample an element of 'vector' based on their values, input needs to be float
    probDeath=vector/sum(vector) # (larger values have higher prob of being sampled)
    r=np.random.choice(probDeath, p=probDeath)
    ind=np.where(probDeath==r)[0][0]
    return [vector[ind], ind]

def round_up(value):
    try:
        return math.ceil(float(value))
    except ValueError:
        return value

def simulate_extinction_and_status_change(delta_t,list_of_all_current_species_statuses,species_list,outdir,qmatrix_dict,status_change=False,dynamic_qmatrix=True):
	# write the species name and the index to a separate txt file
    # this is necessary because two years get lost when summarizing the simulated data into year bins
    #final_year = final_year+2
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    iucn_code = {'LC':0, 'NT':1, 'VU':2, 'EN':3, 'CR':4}
    species_cat_arg = [iucn_code[x] for x in list_of_all_current_species_statuses]
    species_cat = np.array(species_cat_arg)+0
    current_state=species_cat
    root = -delta_t
    n_columns = abs(root)
    status_array = np.tile(np.array([current_state]).transpose(), (1, n_columns))
    species_list_log = open("%s/species_list.txt" %outdir, "w")
    species_list_file=csv.writer(species_list_log, delimiter='\t')
    #print(q_matrix)
    extinction_time_list = []
    for taxon in range(len(species_list)):     
        species = species_list[taxon]
        if dynamic_qmatrix:
            q_matrix = qmatrix_dict[species].copy()
        else:
            q_matrix = qmatrix_dict
        #print(q_matrix)
        if not status_change:
            # if the status change is deactivated set all status shift rates in the q-matrix to 0
            q_matrix[0:5,0:5] = 0.0
            np.fill_diagonal(q_matrix, -np.sum(q_matrix,axis=1))
        species_list_file.writerow([taxon,species])
        current_state = int(status_array[taxon,0])
        D = -np.diagonal(q_matrix)
        extinction = 'extant'
    	#start at 1 so that the current status will always be protected/kept stable
        t=0
        while t < abs(root):
    		# draw randomely when the next change is going to happen, based on the negative inverse diagonal value of the q-matrix (corrected by scale)
            delta_t =np.random.exponential(1./(D[current_state]))
    		# fill all columns between t and delta_t with the current_status (gets updated after this)
            status_array[taxon,np.arange(int(t),min(int(t+delta_t), abs(root)))] = current_state
            t = min(t+delta_t, abs(root))
    		# get the q_row of the current status containing all change rates
            Qrow = q_matrix[current_state]+0
    		# eliminate the negative value in the diagonal and set to 0 since we don't want to sample this one
            Qrow[Qrow<0]= 0 #=1-sum(Qrow[Qrow>0])
            [rate, ind] = random_choice_P(Qrow)
            current_state = ind
            if ind == 5:
                # go extinct
                time_of_extinction = delta_t+t
                if time_of_extinction > abs(root):
                    break
                else:
                    status_array[taxon,np.arange(int(time_of_extinction), abs(root))] = 5
                    extinction = time_of_extinction
                    break
        extinction_time_list.append(extinction)
    a = np.array(species_list)
    b = np.array(extinction_time_list)
    extinction_array = np.stack((a, b))
    list_ex = list(map(round_up, extinction_array[1]))
    extinctions_per_year = Counter(list_ex)
    extinctions_per_year = dict(extinctions_per_year)
    return status_array, extinction_array, extinctions_per_year
    
def get_dtt_array_from_extinction_per_year_dict(extinction_dict_sim_out,current_year,final_year):
    year_bins = np.arange(current_year,final_year+1)
    ext_count= [extinction_dict_sim_out[i] if i in extinction_dict_sim_out.keys() else 0 for i in year_bins]
    diversity_through_time = sum(extinction_dict_sim_out.values()) - np.cumsum(ext_count)
    diversity_through_time
    return diversity_through_time

def run_multi_sim(n_rep,delta_t,species_list_status,dd_probs,qmatrix_dict_list,outdir,all_lc=False,status_change=True,dynamic_qmatrix=True):
    iucn_code = {'LC':0, 'NT':1, 'VU':2, 'EN':3, 'CR':4}
    current_spec_list = species_list_status.species.values
    current_status_list = species_list_status.current_status.values
    extinct_per_year_array = np.zeros([n_rep,delta_t+1])
    te_array = np.zeros((len(species_list_status),n_rep+1)).astype(object)
    te_array[:,0]=current_spec_list
    status_through_time_dict = {}
    status_through_time = np.zeros([6,delta_t+1,n_rep])
    for n in range(n_rep):
        sys.stdout.write('\rRunning simulation rep %i/%i' %(n+1,n_rep))
        if dynamic_qmatrix:
            qmatrix_dict = qmatrix_dict_list[n]
        else:
            qmatrix_dict = qmatrix_dict_list
    
        # these are the simulator functions that need to be repeated
        # new modeling of DD species every rep
        current_status_list_new_dd = current_status_list.copy()
        dd_indices = np.where([current_status_list_new_dd=='DD'])[1]    
        dd_prob_vector = dd_probs.T[n]
        if all_lc:
	        new_draws = np.array(['LC']*len(dd_indices))
        else:
            new_draws = np.random.choice(['LC','NT','VU','EN','CR'], size=len(dd_indices), replace=True, p=dd_prob_vector)
        current_status_list_new_dd[dd_indices] = new_draws
        # new modeling of NE species every rep
        status_count_dict = Counter(current_status_list)
        counts = np.array([status_count_dict[key] for key in status_count_dict.keys() if key not in ['DD','NE']])
        ne_probs = counts/sum(counts)
        status_array_count = [key for key in status_count_dict.keys() if key not in ['DD','NE']]        
        ne_indices = np.where([current_status_list_new_dd=='NE'])[1]
        if all_lc:
	        new_draws = np.array(['LC']*len(ne_indices))
        else:
            new_draws = np.random.choice(status_array_count, size=len(ne_indices), replace=True, p=ne_probs)
        current_status_list_new_dd[ne_indices] = new_draws
    
        future_status_array, extinction_array, extinction_dict = simulate_extinction_and_status_change(delta_t,current_status_list_new_dd,current_spec_list,outdir,qmatrix_dict,status_change=status_change,dynamic_qmatrix=dynamic_qmatrix)
        # diversity_through time array
        for year in extinction_dict.keys():
            if not year == 'extant':
                pos = int(year)
                extinct_per_year_array[n,pos] = extinction_dict[year]
        # status through time array
        year = 1
        for column in future_status_array.T:
            the_dict = dict(Counter(list(column)))
            for status in range(6):
                key = '%s_%s'%(str(year),str(status))
                status_through_time_dict.setdefault(key,[])
                if status in the_dict.keys():
                    status_through_time_dict[key].append(the_dict[status])
                else:
                    status_through_time_dict[key].append(0)
            year += 1
        status_list = [iucn_code[status] for status in current_status_list_new_dd]
        counts_present = dict(Counter(status_list))
        for i in counts_present.keys():
            status_through_time[i,0,n]=counts_present[i]
        # fill te_array
        for species in range(len(extinction_array[1])):
            extinction_date = ''
            if extinction_array[1][species] == 'extant':
                extinction_date = np.nan
            else:
                extinction_date = np.round(float(extinction_array[1][species]),2)
            te_array[species,n+1]=extinction_date
    # finish diversity through time array
    diversity_through_time = sum(extinction_dict.values()) - np.cumsum(extinct_per_year_array,axis=1)
    # finish status_through_time_array
    for key in status_through_time_dict:
        x = int(key.split('_')[1])
        y = int(key.split('_')[0])
        status_through_time[x,y,:] = status_through_time_dict[key]
    # write data-objects to output folder
    #with open(os.path.join(outdir,'diversity_through_time.pkl'), 'wb') as f:
    #    pickle.dump(diversity_through_time, f, pickle.HIGHEST_PROTOCOL)
    #with open(os.path.join(outdir,'te_array.pkl'), 'wb') as f:
    #    pickle.dump(te_array, f, pickle.HIGHEST_PROTOCOL)
    #with open(os.path.join(outdir,'status_through_time.pkl'), 'wb') as f:
    #    pickle.dump(status_through_time, f, pickle.HIGHEST_PROTOCOL)
    return diversity_through_time,te_array,status_through_time

def calcHPD(data, level):
    assert (0 < level < 1)
    d = list(data)
    d.sort()
    nData = len(data)
    nIn = int(round(level * nData))
    if nIn < 2 :
        raise RuntimeError("not enough data")
    i = 0
    r = d[i+nIn-1] - d[i]
    for k in range(len(d) - (nIn - 1)):
         rk = d[k+nIn-1] - d[k]
         if rk < r :
             r = rk
             i = k
    assert 0 <= i <= i+nIn-1 < len(d)
    return (d[i], d[i+nIn-1])


def write_r_scripts(output_folder):
	
    script_1_content = """
    library(rredlist)
    
    args = commandArgs(trailingOnly = TRUE)
    taxon_group = args[1]
    group_rank = args[2]
    iucn_key = args[3]
    outdir = args[4]
    exclude_extinct = FALSE
    
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
    
    # exclude extinct taxa if needed
    if (exclude_extinct){
      boolean = !grepl('EX|EW',status_list)
      taxon_list = taxon_list[boolean]
      group_list = group_list[boolean]
      status_list = status_list[boolean]
      taxon_id = taxon_id[boolean]
    }

    # remove all non-species level identifications
    boolean = !grepl('subsp.|ssp.|subpopulation|Subpopulation',taxon_list)
    taxon_list = taxon_list[boolean]
    group_list = group_list[boolean]
    status_list = status_list[boolean]
    taxon_id = taxon_id[boolean]
    
    # select target taxa
    selected_taxon_list = taxon_list[group_list==taxon_group]
    selected_ids = taxon_id[group_list==taxon_group]
    final_sorted_taxon_list = selected_taxon_list
    #final_taxon_list = as.data.frame(cbind(selected_taxon_list,selected_ids))
    #final_sorted_taxon_list = final_taxon_list[order(final_taxon_list$selected_taxon_list),]
    write.table(final_sorted_taxon_list,file=paste0(outdir,'/',taxon_group,"_species_list.txt"), quote=F,row.names=F,sep='	',col.names = FALSE)
    
    
    # get historic data __________________________
    # create new dataframe with species as first column
    historic_assessments = selected_taxon_list
    historic_assessments = as.data.frame(historic_assessments)
    colnames(historic_assessments) = c('species')
    # find historic assessments and fill into dataframe
    counter = 1
    for (i in seq(1, length(selected_taxon_list), 1)){
      species = selected_taxon_list[i]
      species_id = selected_ids[i]
      print(paste0('Downloading IUCN history: species ',counter, ' of ',length(selected_taxon_list)))
      #print(species)
      row_id = which(historic_assessments$species == species)
      hist_data <- rl_history(id=species_id,key=iucn_key)
      for (year in hist_data$result$year){
        id = which(hist_data$result$year == year)
        #some species have multiple assignments for some years
        if (length(hist_data$result$code[id])>1){
          historic_assessments[row_id,year] <- hist_data$result$code[id][1]
        }
        else{
          historic_assessments[row_id,year] <- hist_data$result$code[id]
        }
      }
      counter = counter+1
    }
    write.table(historic_assessments,file=paste0(outdir,'/',taxon_group,"_iucn_history.txt"), quote=F,row.names=F,sep='	')

    #___________________________________    
    """
    
    script_2_content = """
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
    write.table(species_cat_df,file=paste0(outdir,'/current_status_missing_species.txt'),quote=FALSE,sep = '\t',col.names = FALSE,row.names=FALSE)
    """
    
    rscript_out = os.path.join(output_folder,'rscripts')
    if not os.path.exists(rscript_out):
        os.makedirs(rscript_out)

    with open(os.path.join(rscript_out,'get_iucn_status_data_and_species_list.r'),'w') as file:
        file.write(script_1_content)
        file.close()
    with open(os.path.join(rscript_out,'get_current_iucn_status_missing_species.r'),'w') as file:
        file.write(script_2_content)
        file.close()