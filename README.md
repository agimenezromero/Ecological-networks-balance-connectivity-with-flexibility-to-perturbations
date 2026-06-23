# Ecological networks balance connectivity with flexibility to perturbations

This repository contains the code associated with the manuscript: *Ecological networks balance connectivity with flexibility to perturbations*

The study analyzes ecological interaction networks from the Web of Life and Mangal databases. The corresponding data, processed networks, computed results, and figure inputs are archived on [Zenodo](https://doi.org/10.5281/zenodo.20670601).

## Abstract
The complexity--stability debate in ecology remains unresolved in part because its empirical basis is limited. Most evidence for the predicted decline of connectance with species richness comes from food webs, leaving unclear whether this pattern extends across the full spectrum of ecological interactions. Moreover, existing results remain conceptually unresolved: connectance decreases with diversity, yet both the total number of interactions and the number of interactions per species increase. Here, we analyze 1,500 ecological interaction networks spanning diverse habitats and interaction types. We show that these patterns are broadly shared across ecological interaction networks and can be interpreted through a recent theory of information dynamics in complex networks, in which sparsity is favored by a trade-off between signal propagation and response diversity. Our results suggest that the structural component of the debate may indeed reflect a general architectural regularity of ecological communities rather than a contradiction between theory and nature.

## Repository contents

This GitHub repository contains only the analysis code. The data required to run the analyses are available from the Zenodo record associated with the manuscript.

The main files are:

```text
.
├── network_efficiency_lib.py
├── null_model_study.py
├── perturbation_study.py
├── Figures.ipynb
├── README.md
└── LICENSE
```

The code is organized as follows:

- `network_efficiency_lib.py`: shared functions used across the analyses;
- `null_model_study.py`: script to compute null-model network efficiencies;
- `perturbation_study.py`: script to compute link-addition and link-removal perturbation analyses;
- `Figures.ipynb`: notebook used to reproduce the manuscript and Supplementary Information figures.

## Data

The data are not stored in this GitHub repository.

To reproduce the analyses, download the corresponding Zenodo archive and place the relevant folders in the root directory of this repository. The expected structure is:

```text
.
├── Data/
│   ├── networks/
│   │   └── *.graphml
│   ├── mangal/
│   ├── web-of-life/
│   └── processed/
│       └── network_metrics.csv
├── Results/
│   └── efficiency/
│       ├── null_models/
│       └── perturbations/
├── network_efficiency_lib.py
├── null_model_study.py
├── perturbation_study.py
└── Figures.ipynb
```

The Zenodo archive contains:

- raw ecological-network data from Web of Life and Mangal;
- the final set of ecological networks used in the study, stored in GraphML format;
- processed structural network metrics;
- precomputed outputs required to reproduce the figures;
- the complete data release associated with the manuscript.

The final analysed dataset contains 1,569 ecological interaction networks. The efficiency analyses were performed on the subset of 1,304 networks with more than five species.

## Code

The main analysis functions are defined in:

```text
network_efficiency_lib.py
```

This file includes functions to:

- compute the thermodynamic spectrum of a network;
- compute network efficiency across propagation scales;
- estimate diffusion time;
- generate Erdős--Rényi null networks;
- generate configuration-model null networks;
- generate Curveball fixed-marginal null networks for bipartite networks;
- perform link-addition and link-removal perturbation analyses.

The null-model and perturbation scripts import functions from this shared library.

## Software requirements

The analyses were run in Python. The main dependencies are:

```text
numpy
pandas
networkx
```

A minimal conda environment can be created with:

```bash
conda create -n network-efficiency
conda activate network-efficiency
conda install -c conda-forge numpy pandas networkx
```

Alternatively, install the dependencies in any existing Python environment or with `pip`.

## Reproducing the analyses

All commands below assume that they are run from the root directory of the repository after downloading and placing the Zenodo data folders in the expected location.

### 1. Compute null-model efficiencies

Null-model analyses can be run with:

```bash
python null_model_study.py ER
python null_model_study.py CM
python null_model_study.py CB
```

where:

- `ER` is the Erdős--Rényi null model;
- `CM` is the configuration-model null model;
- `CB` is the Curveball fixed-marginal null model.

For bipartite networks, the ER and CM randomizations preserve the bipartite structure. The Curveball null model is defined only for bipartite networks and preserves row and column sums of the bipartite incidence matrix.

The analyses use 100 null realizations per empirical network unless otherwise specified in the code.

### 2. Compute perturbation analyses

Perturbation analyses can be run with:

```bash
python perturbation_study.py <network_index>
```

where `<network_index>` is the index of the network in the sorted list of GraphML files in `Data/networks/`.

The perturbation analysis randomly adds or removes a given percentage of links and compares the resulting efficiency with that of the empirical network. For bipartite networks, link additions are restricted to cross-guild interactions.

For cluster execution, this script can be submitted as an array job, with one index per network.

### 3. Reproduce figures

The manuscript and Supplementary Information figures can be reproduced with the notebook:

```text
Figures.ipynb
```

The notebook reads the processed data and computed results from the Zenodo archive and regenerates the figures used in the manuscript and Supplementary Information.

## Notes on null models

The study uses a hierarchy of null models.

The Erdős--Rényi null model preserves the number of nodes and expected connectance. For bipartite networks, links are placed only between the two node sets.

The configuration-model null model preserves the degree sequence. For bipartite networks, rewiring is restricted to links between partitions, preserving the degree of each species while maintaining the bipartite structure.

The Curveball null model is applied only to bipartite networks. It randomizes the bipartite incidence matrix through repeated pairwise trades while preserving both row and column sums. Thus, each species retains its empirical number of interaction partners within its guild.

## Reproducibility notes

Some computations, especially null-model and perturbation analyses, can be computationally expensive. The Zenodo archive includes the precomputed outputs used to generate the manuscript figures. The code in this repository can be used to recompute these outputs from the final GraphML networks.

The scripts sort the list of GraphML files before running analyses so that command-line indices correspond to the same networks across systems.

## Data availability

The original ecological interaction networks are available from:

- [Web of Life](https://www.web-of-life.es/)
- [Mangal](https://mangal.io/)

The final processed data release associated with this manuscript is available on [Zenodo](https://doi.org/10.5281/zenodo.20670601).

Zenodo DOI:
```text
    https://doi.org/10.5281/zenodo.20670601
```

## Citation

If you use this code, please cite the associated manuscript:

```text
Giménez-Romero, À. et al. Ecological networks balance connectivity with flexibility to perturbations.
```

## License

Please see the `LICENSE` file for reuse terms.
