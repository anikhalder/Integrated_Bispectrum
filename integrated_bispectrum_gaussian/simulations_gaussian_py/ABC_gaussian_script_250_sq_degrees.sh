#!/bin/bash

python A_healpy_patches_producer_gaussian_del.py 250
python B_treecorr_patches_correlator_gaussian_del.py 250
python C_patches_analyser_gaussian_del.py 250
