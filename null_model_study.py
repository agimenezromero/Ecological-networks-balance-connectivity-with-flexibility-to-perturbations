import network_efficiency_lib

import numpy as np

import glob
import sys

filenames = glob.glob("Data/networks/*.graphml")

τs = np.logspace(-2, 2, 100)

null_model = sys.argv[
    1
]  # Get the null model type from command line argument ("ER", "CM", or "CB")

print("\nComputing null models...\n")

network_efficiency_lib.compute_network_efficiency_null_model(
    filenames, τs, model=null_model, num_randomizations=100, overwrite=False
)
