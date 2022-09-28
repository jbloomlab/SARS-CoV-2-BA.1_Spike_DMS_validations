
"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import pandas as pd

# Configuration  --------------------------------------------------------------

configfile: 'config.yaml'

# Functions -------------------------------------------------------------------

def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Target rules ---------------------------------------------------------------

localrules: all

rule all:
    input:
        'results/summary/virus_titers_for_mAb_escape.md',
        'results/summary/virus_titers_functional_mutants.md',
        'results/summary/spike_neutralization.md',
        'results/summary/VSVG_neutralization.md',
        'results/summary/Lycov1404_yeast_lenti_dms_comparison.md'
        


# Rules ---------------------------------------------------------------------

rule get_virus_titers_mAb:
    """calculate virus titers for antibody escape mutants"""
    input:
        config['virus_titers_antibody']
    output:
        nb_markdown=nb_markdown('virus_titers_for_mAb_escape.ipynb')
    params:
        nb='virus_titers_for_mAb_escape.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule get_virus_titers_functional:
    """calculate virus titers for functional mutants"""
    input:
        config['virus_titers_functional']
    output:
        nb_markdown=nb_markdown('virus_titers_functional_mutants.ipynb')
    params:
        nb='virus_titers_functional_mutants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule plot_neuts_spike:
    """plot neut curves for spike pseudotyped virus"""
    input:
        depletion_neuts=config['mAb_neuts']
    output:
        nb_markdown=nb_markdown('spike_neutralization.ipynb')
    params:
        nb='spike_neutralization.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule plot_neuts_VSVG:
    """plot neut curves for VSV-G pseudotyped virus"""
    input:
        depletion_neuts=config['mAb_neuts_VSVG']
    output:
        nb_markdown=nb_markdown('VSVG_neutralization.ipynb')
    params:
        nb='VSVG_neutralization.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule compare_yeastDMS_vs_lentiDMS:
    """compare yeast and lentivirus DMS for Ly-CoV1404"""
    input:
        yest_DMS=config['yeast_dms_lycov1404_Star'],
        lenti_DMS=config['lenti_dms_lycov1404']
    output:
        nb_markdown=nb_markdown('Lycov1404_yeast_lenti_dms_comparison.ipynb')
    params:
        nb='Lycov1404_yeast_lenti_dms_comparison.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

