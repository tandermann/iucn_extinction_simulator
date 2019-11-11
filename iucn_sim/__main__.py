# encoding: utf-8

#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

import os
import sys
from argparse import ArgumentParser
import logging
import warnings
from . import __version__
import importlib


__author__ = "Tobias Andermann"

# List of all subcommands. A module of the given name must exist and define
# add_arguments() and main() functions.

COMMANDS = [
		'get_rates',
		'run_sim'
]


def main(arguments=None):
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	parser = ArgumentParser(description=__doc__, prog='iucn_sim')
	parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

	subparsers = parser.add_subparsers()
	for command_name in COMMANDS:
		module = importlib.import_module('.' + command_name, 'iucn_sim')
		subparser = subparsers.add_parser(command_name,
			help=module.__doc__.split('\n')[1], description=module.__doc__)
		subparser.set_defaults(func=module.main)
		module.add_arguments(subparser)
        
	args = parser.parse_args(arguments)
	if not hasattr(args, 'func'):
		parser.error('Please provide the name of a subcommand to run')
	else:
		args.func(args)


if __name__ == '__main__':
	main()
