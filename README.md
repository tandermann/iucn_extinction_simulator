<img src="https://github.com/dsilvestro/PyRate/blob/master/pyrate_lib/PyRate_logo1024.png" align="left" width="80">

# IUCN Extinction Simulator

A program for simulating future extinctions and extinction rates for a given set of species, based on IUCN threat assessments.

[![downloads](https://anaconda.org/bioconda/iucn_sim/badges/downloads.svg)](https://anaconda.org/bioconda/iucn_sim)
[![updated](https://anaconda.org/bioconda/iucn_sim/badges/latest_release_date.svg)](https://anaconda.org/bioconda/iucn_sim)
[![conda](https://anaconda.org/bioconda/iucn_sim/badges/installer/conda.svg)](https://anaconda.org/bioconda/iucn_sim)
[![conda](https://anaconda.org/bioconda/iucn_sim/badges/license.svg)](https://anaconda.org/bioconda/iucn_sim) 


## Installation

`iucn_sim` is available as a conda package, which helps installing all required Python and R dependencies, without you having to worry about taking care of this manually.
The conda package manager creates an extremely light-weight virtual environment that enables anybody to run `iucn_sim` on their computer, independently of the operating system and independently of any previous installations of Python and R.

1. Download [miniconda](https://docs.conda.io/en/latest/miniconda.html) for your operating system.

2. Once miniconda is installed, open a command line terminal (e.g. `Terminal` on macOS). Windows users will need to open the freshly installed **Anaconda Powershell Prompt** instead of the regular Command Prompt for this purpose.

3. Add the conda-forge and bioconda channels to conda, where some of the required packages are hosted:
	- `conda config --add channels conda-forge`
	- `conda config --add channels bioconda`

4. Install `iucn_sim` and all it's dependencies by typing the following into the command line, **but see the recommendation below before you do so!**
	- `conda install iucn_sim`

	> Recommendation: Installing iucn_sim with the command above is the simplest solution. However, to be safe that you don't change the standard path on your computer to a different R or Python version than what you have been previously using, it is recommendable to create your own environment for `iucn_sim`. This is not any more complicated, except that everytime before you use iucn_sim, you will need to connect to the created environment (just a single command, see below).
	> - In order to install `iucn_sim` into a new environment, simply type `conda create -n iucn_sim_env iucn_sim`.
	> - Now, everytime you want to use `iucn_sim`, type `conda activate iucn_sim_env` to connect to the environment (in some cases this command doesn't work, but you need to type `source activate iucn_sim_env` instead). Once connected to the virtual environment you can repeatedly use `iucn_sim`. When you are done and want to disconnect from the environment, type `conda deactivate`.

5. Test if installation worked by typing `iucn_sim -h` (if you created a virtual environment, you need to be connected to it for any `iucn_sim` command to work). This should show an overview of the available arguments of `iucn_sim`.

6. If step 4 caused an error something went wrong along the way. If you are a Linux or Mac user, you can instead install `iucn_sim` by [downloading this GitHub repo](https://github.com/tobiashofmann88/iucn_extinction_simulator/archive/master.zip) and building the software by typing `python setup.py install`. You will need to make sure yourself that Python3 and R are installed, including the R package `rredlist` and several Python packages, which can be installed with `pip install NAME_OF_PACHAGE` (the program will tell you which packages need to be installed).

## Running `iucn_sim`

Once installed, `iucn_sim` will be installed in your standard path, so you can simply type `iucn_sim` in your command line (use **Anaconda Powershell Prompt** if you are a Windows user) to call the program. (Again: If you installed `iucn_sim` in a separate environment, first connect to the environment by typing `conda activate iucn_sim_env` to be able to use `iucn_sim`).

- `iucn_sim -h` --> Open help page showing available functions

- `iucn_sim get_iucn_data -h` --> Open help page for specific function

The -h command will show and explain all available flags for each specific `iucn_sim` function. An example command could look like this:


## Quick tutorial

To fully understand the methodology behind `iucn_sim` we recommend you to have a look at the published `iucn_sim` manuscript at [https://doi.org/10.1101/2019.12.16.878249](https://doi.org/10.1101/2019.12.16.878249) (Andermann et al., 2020).

In the following tutorial we will predict future extinctions and extinction rates for all species of the order Carnivora (**target species list**). We will use the the whole class Mammalia as **reference group**. In `iucn_sim` the target species list contains a list of species names for which you want to simulate future extinctions. The reference group on the other hand is a group of species which are being used to estimate **status transition rates** (i.e. the rates of how often species change from one IUCN status to another) based on the IUCN history of the group. This reference group should be sufficiently large (>1,000 species) to increase the accuracy of the estimated transition rates.

This tutorial uses pre-compiled IUCN data, without requiring an IUCN API key. If you plan on running `iucn_sim` on your own target species list and reference group, you will first need to apply for an IUCN key (see information below). However, `iucn_sim` has access to a range of pre-compiled reference groups, which enable processing without requiring an IUCN API key (see overview of precompiled groups [here](https://github.com/tobiashofmann88/iucn_extinction_simulator/tree/master/data/precompiled/iucn_history), but no need to download these).

### Required input data:

The only file you need for running this tutorial is the [carnivora_gl.txt](https://github.com/tobiashofmann88/iucn_extinction_simulator/blob/master/data/precompiled/gl_data/carnivora_gl.txt) file stored in the `data/precompiled/gl_data` folder in this GitHub repo. You can download it by clicking on the file-link and copying the content to your text-editor or by executing `wget https://github.com/tobiashofmann88/iucn_extinction_simulator/blob/master/data/precompiled/gl_data/carnivora_gl.txt` in your command line, in case you have `wget` installed.

This file contains a list of all Carnivora species (IUCN 2019-v2), including 100 generation length (GL) estimates for each species (scaled in years). The purpose of having 100 estimates per species is to include the uncertainty of the GL value for those species where missing GL data was modeled based on phylogenetic imputation.

### Get IUCN data:

The first step is downloading available IUCN data, which includes the IUCN history of the reference group, the current status information for all species in the target species list, and a list of possibly extinct species belonging to the reference group.

(**Remember to activate your `iucn_sim` environment first**, in case you installed it in its own environment: `conda activate iucn_sim_env` or `source activate iucn_sim_env`)

```
iucn_sim get_iucn_data \
	--reference_group mammalia \
	--target_species_list data/precompiled/gl_data/carnivora_gl.txt \
	--outdir data/iucn_sim_output/carnivora/iucn_data/
```

### Estimate status transition rates

Now we want to estimate the rates of how often any type of status change occurs in the IUCN history of the reference group. This is done by sampling these rates from the counts of each type of status change, using a Markov chain Monte Carlo algorithm (MCMC). Additionally to the **status transition rates** we also estimate the rates at which species of any given status become extinct (**EX transition rates**). For estimating these rates `iucn_sim` offers two different methods.

- **EX mode 0** (sensu Mooers et al., 2008): This method for estimating EX transition rates applies pre-defined IUCN extinction probabilities, which are defined for the criterion E of threatened species (see [IUCN Red List guidelines](http://cmsdocs.s3.amazonaws.com/RedListGuidelines.pdf)). These probabilties are extrapolated for the non-threatened statuses Least Concern (LC) and Near Threatened (NT). Finally, if GL data are provided (as in the [carnivora_gl.txt](https://github.com/tobiashofmann88/iucn_extinction_simulator/blob/master/data/precompiled/gl_data/carnivora_gl.txt) example file), these data are being considered when calculating the EX transition rates for statuses Endangered (EN) and Critically Endangered (CR), as intended per IUCN definition.

	```
	iucn_sim transition_rates \
		--species_data data/iucn_sim_output/carnivora/iucn_data/species_data.txt \
		--iucn_history data/iucn_sim_output/carnivora/iucn_data/MAMMALIA_iucn_history.txt \
		--outdir data/iucn_sim_output/carnivora/transition_rates_0 \
		--extinction_probs_mode 0 \
	```

- **EX mode 1** (sensu Monroe et al., 2019): In this method EX transition rates are being estimated from the observed transitions in the IUCN history of the reference group towards the statuses Extinct in the Wild (EW) and Extinct (EX). The estimation of these rates is done in the same manner as for the other **status transition rates**. Additionally the user can provide a list of possibly extinct taxa (PEX), which is automatically downloaded by the `iucn_sim get_iucn_data` function ([source](https://nc.iucnredlist.org/redlist/content/attachment_files/2020_1_RL_Stats_Table_9.pdf)), and is applied to correct the usually underestimated number of observed extinctions in the IUCN history.

	```
	iucn_sim transition_rates \
		--species_data data/iucn_sim_output/carnivora/iucn_data/species_data.txt \
		--iucn_history data/iucn_sim_output/carnivora/iucn_data/MAMMALIA_iucn_history.txt \
		--outdir data/iucn_sim_output/carnivora/transition_rates_1 \
		--extinction_probs_mode 1 \
		--possibly_extinct_list data/iucn_sim_output/carnivora/iucn_data/possibly_extinct_reference_taxa.txt
	```


### Simulate future extinctions and estimate species-specific extinction rates

In this final step of `iucn_sim` we simulate future status changes and extinctions for the species in our target species list (all Carnivora in this case) over a specified time frame. From the simulated extinciton dates of individual species over several simulation replicates, `iucn_sim` estimates the extinction rates of each species. These rate estimates inherently contain the probabilities of a given species to change conservation status, as well as the GL data for this species (in case of EX mode 0). A minimum of 10,000 simulation replicates is recommended for accurate extinction rate estimates from the simulated data.

You can turn off the rather time intensive estimation of species-specific extinction rates by setting `--extinction_rates 0`, in case you are only interested in the projected diversity. Otherwise turn it on by using `--extinction_rates 1` (default).

```
iucn_sim run_sim --input_data data/iucn_sim_output/carnivora/transition_rates_1/simulation_input_data.pkl --outdir data/iucn_sim_output/carnivora/future_sim_1 --n_years 100 --n_sim 10000 --extinction_rates 1
```

## Apply for IUCN API token

To use the full functionality of `iucn_sim` you will have to apply for an [IUCN API token](https://apiv3.iucnredlist.org/api/v3/token). This key is necessary to download data from IUCN, which is done internally in `iucn_sim`. It is easy o apply for an API key, just [follow this link](https://apiv3.iucnredlist.org/api/v3/token), it will then take up to a couple of days before you receive your API key. Once you have received your IUCN token, provide it when using the `get_rates` function with the `--iucn_key` flag.

## Options without IUCN API token
If for some reason you have problems obtaining an IUCN API key or don't want to wait until receiving it, there are several **options of avoiding IUCN API key usage**. However, the available taxon options are limited and we strongly recommend to apply for your own API key.

These are your options without an IUCN API key: You can e.g. run the tutorial above or you can run your own commands by choosing the name of one of the [pre-compiled reference groups](https://github.com/tobiashofmann88/iucn_extinction_simulator/tree/master/data/precompiled/iucn_history), as the `--reference_group` in `iucn_sim get_iucn_data`. In that case you can turn off the downloading of current status data for your target species list by setting `--target_species_list 0` or use the same taxa as those in the reference group by setting `--target_species_list 1`. Alternatively you can also provide a txt file for `--target_species_list` with a subset of the taxa names of the reference group (as in the tutorial above).


## References

**Andermann et al. 2020**. iucn_sim - Improved Predictions of Future Extinctions Using IUCN Status Assessments. - biorxiv, doi: [10.1101/2019.12.16.878249](https://doi.org/10.1101/2019.12.16.878249).

**Mooers, A. Ø. et al. 2008**. Converting endangered species categories to probabilities of extinction for phylogenetic conservation prioritization. - PLoS ONE 3: 1–5, doi: 10.1371/journal.pone.0003700.

**Monroe, M. J. et al. 2019**. The dynamics underlying avian extinction trajectories forecast a wave of extinctions. - Biology Letters 15: 20190633, doi: 10.1098/rsbl.2019.0633.