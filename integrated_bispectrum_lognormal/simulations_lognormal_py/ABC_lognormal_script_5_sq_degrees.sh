#!/bin/bash

python A_healpy_patches_producer_lognormal_del.py 5
python B_treecorr_patches_correlator_lognormal_del.py 5
python C_patches_analyser_lognormal_del.py 5
