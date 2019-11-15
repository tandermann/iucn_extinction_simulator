# IUCN Extinction Simulator

A program for estimating extinction probabilities and dates for a given set of species, based on IUCN threat assessments.


## Installation

`iucn_sim` is available as a conda package, which helps installing all required Python and R dependencies, without you having to worry about taking care of this yourself.
The conda package manager creates an extremely light-weight virtual environment that enables anybody to run `iucn_sim` on their computer, independently of the operating system and independently of any previous installations of Python and R.

1. Download [miniconda](https://docs.conda.io/en/latest/miniconda.html) for your operating system.

2. Once miniconda is installed, open a command line terminal. Windows users will need to open the freshly installed **Anaconda Powershell Prompt** instead of the regular Command Prompt for this purpose.

3. Install `iucn_sim` and all it's dependencies by typing `conda install iucn_sim` into the command line.

4. Test if installation worked by typing `iucn_sim -h`. This should show an overview of the available arguments of `iucn_sim`.

5. If step 4 caused an error something went wrong along the way. If you are a Linux or Mac user, you can instead install `iucn_sim` by [downloading this GitHub repo](https://github.com/tobiashofmann88/iucn_extinction_simulator/archive/master.zip) and building the software by typing `python setup.py install`. You will need to make sure yourself that Python3 and R are installed, including the R package `rredlist`.

## Required input data:

For running `iucn_sim` you will need to provide the following input:

- Species list and **generation length (GL)** data (scaled in years and tab-separated). Can contain multiple GL values per species, representing indpendent samples from the uncertainty surrounding the GL estimate. See example at `data/example_data/gl_data_all_mammals.txt`.
- Name of taxon group (e.g. `'Mammalia'`) or alternatively path to text file containing list of species, to be used for estimating rates of change between different IUCN categories (should be sufficiently large group, we recommend > 1,000 species).
- How many years to simulate into the future, e.g. `80` if the next 80 years (starting at present) should be simulated.
- Number of simulation replicates, e.g. `100` will lead to 100 independent simulation replicates. If `N` different  GL values are provided per species, the simulator will by default run `N` simulation replicates, using a different GL value for each replicate.


