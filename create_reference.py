#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import click
from utils.utils import create_results_dir_


@click.command()
@click.option('--raw-reference', type=click.Path(), help='Reference genome.', default="resources/reference_genome/hg19.fa")
@click.option('--reference-name', type=str, help='Project name.', default="TEST")
@click.option('--regions', type=str, help='Bed file with regions to select.', default=".")
@click.option('--output', type=str, help='Project name.', default=".")

def main(raw_reference, project_name, output):
    raw_elements_dir, final_project_dir = create_results_dir_(output)
    