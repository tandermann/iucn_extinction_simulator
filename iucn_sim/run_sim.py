#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run future simulations based on IUCN data and status transition rates

Created on Wed Oct 30 20:59:28 2019
@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import iucn_sim.iucn_sim as iucn_sim


def add_arguments(parser):
    parser.add_argument(
        '--input_data',
        required=True,
        help="Path to 'simulation_input_data.pkl' file created by transition_rates function."
    )
    parser.add_argument(
        '--outdir',
        required=True,
        help="Provide path to outdir where results will be saved."
    )
    parser.add_argument(
        '--n_years',
        default=100,
        help="How many years to simulate into the future."
    )
    parser.add_argument(
        '--n_sim',
        default=10000,
        help="How many simulation replicates to run. At least 10,000 simulations are recommended for accurate rate estimation (default). If the number of simulation replicates exceeds the number of available transition rate estimates (produced by the 'transition_rates' function), these rates will be randomely resampled for the remaining simulations."
    )
    parser.add_argument(
        '--status_change',
        default=1,
        help="Model IUCN status changes in future simulations. 0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--conservation_increase_factor',
        default=1,
        help="The transition rates leading to improvements in IUCN conservation status are multiplied by this factor."
    )
    parser.add_argument(
        '--threat_increase_factor',
        default=1,
        help="Opposite of conservation_increase_factor, multiplies the transition rates leading to worsening in IUCN conservation status."
    )
    parser.add_argument(
        '--model_unknown_as_lc',
        default=0,
        help="Model new status for all DD and NE species as LC (best case scenario). 0=off, 1=on (default=0)."
    )
    parser.add_argument(
        '--until_n_taxa_extinct',
        default=0,
        help="Setting this value will stop the simulations when n taxa have gone extinct. This can be used to simulate the expected time until n extinctions. The value of the --n_years flag in this case will be interpreted as the maximum possible time frame, so set it large enough to cover a realistic time-frame for these extinctions to occur. Set to 0 to disable this function (default=0)."
    )
    parser.add_argument(
        '--extinction_rates',
        default=1,
        help="Estimation of extinction rates from simulation results: 0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--n_gen',
        default=100000,
        help="Number of generations for MCMC for extinction rate estimation (default=100000)."
    )
    parser.add_argument(
        '--burnin',
        default=1000,
        help="Burn-in for MCMC for extinction rate estimation (default=1000)."
    )
    parser.add_argument(
        '--plot_diversity_trajectory',
        default=1,
        help="Plots the simulated diversity trajectory: 0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--plot_status_trajectories',
        default=1,
        help="Plots the simulated IUCN status trajectory: 0=off, 1=on (default=0)."
    )
    parser.add_argument(
        '--plot_histograms',
        default=0,
        help="Plots histograms of simulated extinction times for each species: 0=off, 1=on (default=0)."
    )
    parser.add_argument(
        '--plot_posterior',
        default=0,
        help="Plots histograms of posterior rate estimates for each species: 0=off, 1=on (default=0)."
    )
    parser.add_argument(
        '--plot_status_piechart',
        default=1,
        help="Plots pie charts of status distribution: 0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--seed',
        default=None,
        help="Set starting seed for future simulations."
    )


def main(args):
    simulation_output = iucn_sim.run_sim(
        input_data = args.input_data,
        outdir = args.outdir,
        n_years = args.n_years,
        n_sim = args.n_sim,
        status_change = args.status_change,
        conservation_increase_factor = args.conservation_increase_factor,
        threat_increase_factor = args.threat_increase_factor,
        model_unknown_as_lc = args.model_unknown_as_lc,
        until_n_taxa_extinct = args.until_n_taxa_extinct,
        plot_diversity_trajectory = args.plot_diversity_trajectory,
        plot_status_trajectories = args.plot_status_trajectories,
        plot_histograms = args.plot_histograms,
        plot_status_piechart = args.plot_status_piechart,
        seed = args.seed,
        load_from_file = True
        )

    if args.extinction_rates:
        ext_rates = iucn_sim.estimate_extinction_rates(
            simulation_output._extinction_times,
            args.n_years,
            args.outdir,
            n_gen = args.n_gen,
            burnin = args.burnin,
            plot_posterior = args.plot_posterior,
            seed = simulation_output._seed,
            load_from_file=False # since in this case the input is parsed as an object
            )


