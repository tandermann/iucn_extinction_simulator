#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MCMC-estimation of status transition rates from IUCN record

Created on Mon Oct 28 14:43:44 2019
@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""


import iucn_sim.iucn_sim as iucn_sim

def add_arguments(parser):
    parser.add_argument(
        '--species_data',
        required=True,
        metavar='<path>',
        help="File containing species list and current IUCN status of species, as well as generation length (GL) data estimates if available. GL data is only used for '--extinction_probs_mode 0' ('species_data.txt' output from get_iucn_data function).",
    )
    parser.add_argument(
        '--iucn_history',
        required=True,
        metavar='<path>',
        help="File containing IUCN history of the reference group for transition rate estimation ('*_iucn_history.txt' output of get_iucn_data function)."
    )
    parser.add_argument(
        '--outdir',
        required=True,
        metavar='<path>',
        help="Provide path to outdir where results will be saved."
    )
    parser.add_argument(
        '--extinction_probs_mode',
        default=0,
        metavar='N',
        help="Set to '0' to use the critE EX mode to determine extinction probabilities for each status (e.g. Mooers et al, 2008 approach). Set to '1' to use empirical EX mode, based on the recorded extinction in the IUCN history of the reference group (e.g. Monroe et al, 2019 approach). GL data can only be used in the critE EX mode ('0')."
    )
    parser.add_argument(
        '--possibly_extinct_list',
        default=0,
        metavar='<path>',
        help="File containing list of taxa that are likely extinct, but that are listed as extant in IUCN, including the year of their assessment as possibly extinct ('possibly_extinct_reference_taxa.txt' output from get_iucn_data function). These species will then be modeled as extinct by the esimate_rates function, which will effect the estimated extinction probabilities when chosing `--extinction_probs_mode 1`",
    )
    parser.add_argument(
        '--species_specific_regression',
        action='store_true',
        help='Enables species-specific regression fitting to model LC, NT, and VU extinction probabilities. Only applicable with --extinction_probs_mode 0 (critE mode) and if GL is provided.',
        default=False
    )
    parser.add_argument(
        '--rate_samples',
        default=100,
        metavar='N',
        help="How many rates to sample from the posterior transition rate estimates. These rates will be used to populate transition rate q-matrices for downstream simulations. Later on you can still chose to run more simulation replicates than the here specified number of produced transition rate q-matrices, in which case the `run_sim` function will randomely resample from the available q-matrices (default=100, this is ususally sufficient, larger numbers can lead to very high output file size volumes)."
    )
    parser.add_argument(
        '--n_gen',
        default=100000,
        metavar='N',
        help="Number of generations for MCMC for transition rate estimation (default=100000)."
    )
    parser.add_argument(
        '--burnin',
        default=1000,
        metavar='N',
        help="Burn-in for MCMC for transition rate estimation (default=1000)."
    )
    parser.add_argument(
        '--seed',
        default=None,
        help="Set starting seed for the MCMC."
    )


def main(args):
    tr_rates = iucn_sim.transition_rates(
        species_data = args.species_data,
        iucn_history = args.iucn_history,
        outdir = args.outdir,
        extinction_probs_mode = args.extinction_probs_mode,
        possibly_extinct_list = args.possibly_extinct_list,
        species_specific_regression = args.species_specific_regression,
        rate_samples = args.rate_samples,
        n_gen = args.n_gen,
        burnin = args.burnin,
        seed = args.seed,
        load_from_file = True # the load_from_file option should always be true when run from command line, but can be turned off when running from R or python
        )
    
    