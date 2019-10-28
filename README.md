# IUCN Extinction Simulator

A program for estimating extinction probabilities and dates for a given set of species, based on IUCN threat assessments.

Required input data:

- Species list and **generation length (GL)** data (scaled in years and tab-separated). Can contain multiple GL values per species, representing indpendent samples from the uncertainty surrounding the GL estimate. See example at `data/example_data/gl_data_all_mammals.txt`.
- Name of taxon group (e.g. `'Mammalia'`) or alternatively path to text file containing list of species, to be used for estimating rates of change between different IUCN categories (should be sufficiently large group, we recommend > 3,000 species).
- How many years to simulate into the future, e.g. `80` if the next 80 years (starting at present) should be simulated.
- Number of simulation replicates, e.g. `100` will lead to 100 independent simulation replicates. If `N` different  GL values are provided per species, the simulator will by default run `N` simulation replicates, using a different GL value for each replicate.


