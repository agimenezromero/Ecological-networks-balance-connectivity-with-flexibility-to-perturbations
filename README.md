# Ecological networks balance connectivity with flexibility to perturbations

Code accompanying the manuscript *Ecological networks balance connectivity with flexibility to perturbations* by Àlex Giménez-Romero, Daniel Oro, Jordi Bascompte, and Meritxell Genovart.

This repository contains the scripts and notebooks used to construct and characterize empirical ecological interaction networks, quantify the thermodynamic efficiency of network diffusion, compare empirical networks with null models, assess the effects of link perturbations, and produce the manuscript figures. The associated network data, derived metrics, and precomputed results are archived separately on [Zenodo](https://doi.org/10.5281/zenodo.20670601).

## Abstract

Food webs become progressively sparser as species richness increases, yet the number of interactions grows faster than linearly, placing communities between the limiting expectations of constant degree and constant connectance. Similar scaling has also been documented in other interaction networks, but whether they share a common architectural regime, and what general principle could generate it, remain unresolved. Here, we show that 1,533 ecological networks spanning different interaction types broadly converge on an intermediate link--species scaling regime, with an overall exponent of 1.43. We interpret this pattern through a thermodynamic trade-off in which increasing network connectivity enhances signal propagation while reducing the diversity of independent modes contributing to the network's response to perturbations. Optimizing this trade-off yields a scaling exponent of 1.45, closely matching the empirical estimate, while randomly adding or removing interactions generally reduces thermodynamic efficiency in empirical networks. These results support a propagation--response-diversity trade-off as a plausible systems-level principle underlying how ecological network complexity grows with diversity.

## Repository contents

| File | Description |
| --- | --- |
| `Mangal.ipynb` | Julia notebook for retrieving Mangal networks, exporting network files, and assembling metadata. |
| `studies_networks_dataset.ipynb` | Processes the ecological networks and computes structural network metrics, including connectance, degree statistics, modularity, and diffusion time. |
| `network_efficiency_lib.py` | Python functions for spectral and thermodynamic quantities, diffusion times, ER/CM/Curveball null models, and link perturbations. |
| `null_model_study.py` | Computes thermodynamic efficiency for randomized network ensembles. |
| `perturbation_study.py` | Computes changes in thermodynamic efficiency following random link additions and removals. |
| `network_efficiency.ipynb` | Generates the conceptual illustration of propagation, response diversity, and thermodynamic efficiency (main-text Fig. 2). |
| `Figures.ipynb` | Performs dataset selection, scaling analyses, comparisons with null-model outputs, and plotting for the main text and Supplementary Information. |
| `LICENSE` | GNU General Public License, version 3. |

## Data

The networks originate from [Web of Life](https://www.web-of-life.es/) and [Mangal](https://mangal.io/). The initial collection comprised **1,801 network files**. After applying the study's filtering criteria and removing **36 duplicate networks** identified across the source databases, the final structural analysis comprises **1,533 networks**. Analyses of thermodynamic efficiency use a subset of **1,268 networks** with more than five species.

The [Zenodo archive](https://doi.org/10.5281/zenodo.20670601) includes the underlying network files, metadata, derived structural metrics, and precomputed thermodynamic, null-model, and perturbation outputs. It retains files for networks excluded from the final analyses; the selection and cross-database deduplication are implemented in `Figures.ipynb`. Consequently, the number of GraphML or result files in the archive is not the number of networks included in every analysis.

The analysis notebooks expect the Zenodo `Data/` and `Results/` folders to be placed in the repository's working directory, with paths such as:

```text
Data/
├── networks/                         # Network files (*.graphml)
├── mangal/                           # Source data and network metadata
├── web_of_life/                      # Source data and study references
└── NatComm_data.csv                  # External network comparison data

Results/
├── field_studies_metrics_network.csv # Computed structural metrics
├── Communities_A/                    # Individual degree distributions
└── efficiency/
    ├── null_models/                  # ER, CM and CB results
    └── perturbations/                # Link addition/removal results
```

## Analyses

The common Python implementation is in `network_efficiency_lib.py`. It uses the graph Laplacian to calculate diffusion-based entropy, free-energy changes, thermodynamic efficiency, and characteristic diffusion times. It also implements three network randomization approaches:

- **Erdős–Rényi (ER):** randomly places links while preserving network size and expected interaction density; bipartite randomizations restrict links to the two node sets.
- **Configuration model (CM):** accounts for empirical degree heterogeneity. For topologically bipartite networks, the canonical bipartite configuration model (BiCM) preserves each node's degree *in expectation*; for non-bipartite networks, a configuration-model multigraph is converted into a simple graph by removing self-loops and multiple edges.
- **Curveball (CB):** randomizes topologically bipartite networks through fixed-marginal trades, preserving the degree of every node and the total number of links exactly.

The perturbation experiment compares the efficiency of empirical networks with that of networks obtained by randomly adding or removing a specified percentage of links. The current scripts process the GraphML files they find in `Data/networks/`; the final analysis subset is selected in `Figures.ipynb`.

### Running the scripts

With the Zenodo data in the expected directories, the null-model script accepts one of three model labels:

```bash
python null_model_study.py ER
python null_model_study.py CM
python null_model_study.py CB
```

It uses 100 realizations per network and writes ensemble summaries under `Results/efficiency/null_models/`. CB is applied only when the graph is topologically bipartite.

The perturbation script processes one GraphML network selected by its index in the file list:

```bash
python perturbation_study.py 0
```

It uses 100 realizations for each perturbation percentage from 1% to 10%, saving summaries to `Results/efficiency/perturbations/percentage_<p>/`. The index is currently based on the file order returned by Python's `glob.glob`, which is not guaranteed to be stable across environments; check or explicitly sort the input file list when selecting a particular network.

The figure notebooks use the precomputed results supplied on Zenodo. `network_efficiency.ipynb` generates the conceptual figure, while `Figures.ipynb` assembles the empirical scaling, efficiency, null-model, perturbation, and supplementary analyses. The latter includes the filtering and duplicate-removal steps used to construct the final analytical dataset.

## Software

The analysis code uses Python with NumPy, pandas, NetworkX, SciPy, Matplotlib, seaborn, statsmodels, powerlaw, cartopy, and uncvalue, together with a Jupyter notebook environment. `Mangal.ipynb` uses Julia and the `Mangal`, `DataFrames`, `CSV`, `Graphs`, and `PyCall` packages. The Julia notebook was authored using a Julia 1.12 kernel.

The required Python packages depend on which analyses are run: the numerical scripts use a smaller subset, whereas `Figures.ipynb` includes the plotting and statistical packages listed above.

## Citation

If you use this code or the associated processed dataset, please cite the accompanying manuscript and the archived data record:

- Giménez-Romero, À., Oro, D., Bascompte, J., and Genovart, M. *Ecological networks balance connectivity with flexibility to perturbations*. Manuscript.
- [Associated data and computational results (Zenodo)](https://doi.org/10.5281/zenodo.20670601).

## License

The code in this repository is released under the [GNU General Public License v3.0](LICENSE). Data originating from third-party repositories remain subject to their respective attribution and reuse terms.
