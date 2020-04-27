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

4. Install `iucn_sim` and all it's dependencies by typing the following into the command line, **but see the tip below before you do so!**
	- `conda install iucn_sim`

	> Tip: Installing iucn_sim with the command above is the simplest solution. However, to be safe that you don't change the standard path on your computer to a different R or Python version than what you have been previously using, it is recommendable to create your own environment for `iucn_sim`. This is not any more complicated, except that everytime before you use iucn_sim, you will need to connect to the created environment (just a single command, see below).
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

In the following tutorial we will predict future extinctions and extinction rates for all species of the order Carnivora. We will use the use

#### Get IUCN data:

`iucn_sim get_iucn_data --reference_group mammalia --target_species_list data/precompiled/gl_data/carnivora_gl.txt --outdir data/iucn_sim_output/carnivora/iucn_data/`

`iucn_sim get_iucn_data --reference_group aves --reference_rank class --target_species_list data/precompiled/gl_data/aves_gl.txt --outdir data/iucn_sim_output/aves/iucn_data --iucn_key <IUCN-key>`

#### Estimate status transition rates and extinction probabilities for all taxa

`iucn_sim transition_rates --species_data data/iucn_sim_output/aves/iucn_data/species_data.txt --iucn_history data/iucn_sim_output/aves/iucn_data/AVES_iucn_history.txt --extinction_probs_mode 0 --possibly_extinct_list data/iucn_sim_output/aves/iucn_data/possibly_extinct_reference_taxa.txt --rate_samples 100 --random_seed 1234 --outdir data/iucn_sim_output/aves/transition_rates_0`

#### Run simulations for future

`iucn_sim run_sim --input_data data/iucn_sim_output/aves/transition_rates_0/simulation_input_data.pkl --outdir data/iucn_sim_output/aves/future_sim_0 --n_years 100 --n_sim 100 --extinction_rates 0`

See below for further explanation of the required input.

## Required input data:

For running `iucn_sim` you will need to provide the following input:

- Species list, preferably containing **generation length (GL)** data for all taxa (scaled in years and tab-separated). Can contain multiple GL values per species, representing indpendent samples from the uncertainty surrounding the GL estimate. See example at `data/example_data/gl_data_all_mammals.txt`.
- Setting a reference group: Name of a taxon group (e.g. `'Mammalia'`) or alternatively path to text file containing list of species, to be used for estimating rates of change between different IUCN categories (should be sufficiently large group, we recommend > 1,000 species).


## Apply for IUCN API token

To use the full functionality of `iucn_sim` you will have to apply for an [IUCN API token](https://apiv3.iucnredlist.org/api/v3/token). This key is necessary to download data from IUCN, which is done internally in `iucn_sim`. It is easy o apply for an API key, just [follow this link](https://apiv3.iucnredlist.org/api/v3/token), it will then take a couple of days before you receive your API key. Once you have received your IUCN token, provide it when using the `get_rates` function with the `--iucn_key` flag.

However, if for some reason you have problems obtaining an IUCN token or don't want to wait until receiving it, you can run the `get_rates` function without donwloading data from IUCN. In that case you will need to provide an extra file containing the current IUCN status of all your input species (`--status_list`, but see **Note** below) and you will further need to choose a reference group (`--reference_group`) for which the IUCN history is already pre-compiled ([see available precompiled groups here](https://github.com/tobiashofmann88/iucn_extinction_simulator/tree/master/data/precompiled/iucn_history)).

**Note:** If all species from the `--input_data` file are present in one of the pre-compiled IUCN history files in the [GitHub repo](https://github.com/tobiashofmann88/iucn_extinction_simulator/tree/master/data/precompiled/iucn_history), it is **not necessary** to provide an additional file containing the current IUCN status of all input species (`--status_list`). E.g. in the code example above for running the `get_rates` function, it is technically not necessary to provide an IUCN key, since all species in the input Carnivora file are present in the precompiled IUCN history file for Mammalia. However, due to the ever-changing taxonomy it is not a given that all your species are contained in the precompiled file, even if you have a list of mammal names and use the precompiled Mammalia reference group. Long story short: Make your life easier and [apply for an IUCN token](https://apiv3.iucnredlist.org/api/v3/token)!

## Other user settings:

- `--n_reps` Number of simulation replicates, e.g. `100` will lead to 100 independent simulation replicates. If `N` different  GL values are provided per species, the simulator will by default run `N` simulation replicates, using a different GL value for each replicate.
- `--n_years` How many years to simulate into the future, e.g. `80` if the next 80 years (starting at present) should be simulated.


## Info about IUCN colors and font
https://www.withoutnations.com/portfolio/iucn-red-list/
