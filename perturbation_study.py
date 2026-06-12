import network_efficiency_lib

import numpy as np
import pandas as pd

import os
import glob

import sys

k = sys.argv[1]  # Get the index of the file to process from command line argument

realizations = 100
τs = np.logspace(-2, 2, 100)

filenames = glob.glob("Data/networks/*.graphml")

filename = filenames[int(k)]

id = filename.split("/")[-1].split(".")[0]

for p in range(1, 11):  # Perturbation percentages from 1% to 10%

    output_dir = "Results/efficiency/perturbations/percentage_{}".format(p)

    output_file = "{}/{}.csv".format(output_dir, id)

    if os.path.exists(output_file):
        print(f"{output_file} already exists. Skipping...")
        continue

    print(f"\nProcessing {filename} with {p}% perturbation...\n")

    (
        diff_extra_mean,
        diff_one_less_mean,
        diff_extra_std_err,
        diff_one_less_std_err,
        mean_diff_more_links,
        mean_diff_less_links,
    ) = network_efficiency_lib.compute_efficiency_with_perturbations(
        filename, τs, realizations, percentage=p
    )

    if not os.path.exists("Results/efficiency/perturbations/percentage_{}".format(p)):
        os.makedirs(
            "Results/efficiency/perturbations/percentage_{}".format(p), exist_ok=True
        )

    # Save relative difference results
    df = pd.DataFrame(
        {
            "id": [id],
            "percentage": [p],
            "mean_diff_more_links": [mean_diff_more_links],
            "mean_diff_less_links": [mean_diff_less_links],
        }
    )

    df.to_csv(
        output_file,
        index=False,
    )
