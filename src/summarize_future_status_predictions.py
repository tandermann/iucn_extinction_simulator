import pickle as pkl
import numpy as np

def load_obj(file_name):
    with open(file_name, 'rb') as f:
        return pkl.load(f)

def summarize_future_status_array(future_status_file):
    future_status_list = load_obj(future_status_file)
    future_status_list = np.array(future_status_list)
    status_probs = np.array([np.count_nonzero(future_status_list==i,axis=0)/future_status_list.shape[0] for i in np.arange(future_status_list.max()+1)])
    most_probable_status_per_year = np.argmax(status_probs,axis=0)+1
    return(status_probs,most_probable_status_per_year)



future_status_file = '/Users/tobiasandermann/GitHub/iucnsim/data/iucn_sim/future_simulations_future_status/future_status_array_list.pkl'
status_probs,most_probable_status_per_year = summarize_future_status_array(future_status_file)



import matplotlib.pyplot as plt

plt.plot(status_probs[0, :, :].T)

plt.plot(most_probable_status_per_year.T)
