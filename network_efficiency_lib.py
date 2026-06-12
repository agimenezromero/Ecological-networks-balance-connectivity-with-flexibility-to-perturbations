import numpy as np
import pandas as pd
import networkx as nx

import os, glob


def thermo_spectrum(G, τ):
    """
    inputs : the adjacency matrix, the propagation scale beta (also indicated as tau)
    outputs : the change in entropy (dS), free enrgy (dF) and the eta of network formation
    """

    Ls = np.sort(nx.laplacian_spectrum(G))

    N = len(Ls)

    p = np.exp(-τ * Ls)

    Ls = np.delete(Ls, np.where(p < 10**-12))
    p = np.delete(p, np.where(p < 10**-12))

    Z = np.sum(p)

    p = p / Z

    dF = np.log(N) - np.log(Z)
    dS = np.sum(-p * np.log(p)) - np.log(N)

    eta = (dF + dS) / dF

    return dS, dF, eta


def diffusion_time(G):

    # Diffusion time is the inverse of the smallest non-zero eigenvalue of the Laplacian
    Ls = np.sort(nx.laplacian_spectrum(G))

    # Remove zero eigenvalues (if any)
    Ls = Ls[Ls > 1e-12]

    if len(Ls) == 0:
        return np.inf  # No diffusion possible if no non-zero eigenvalues

    return 1 / Ls[0]


def compute_network_efficiency(filenames, τs):

    # τs = np.logspace(-2, 2, 50)

    if not os.path.exists("Results/efficiency/"):
        os.makedirs("Results/efficiency/")

    for filename in filenames:

        # Print progress
        print(
            "Computing {}/{}".format(filenames.index(filename) + 1, len(filenames)),
            end="\r",
        )

        id = filename.split("/")[-1].split(".")[0]

        G = nx.read_graphml(filename)

        entropy_arr = np.nan * np.ones(len(τs))
        free_energy_arr = np.nan * np.ones(len(τs))
        ηs_arr = np.nan * np.ones(len(τs))

        j = 0
        for τ in τs:

            try:

                dS, dF, η = thermo_spectrum(G, τ)

                if η < 0 or η > 1:
                    raise ValueError(f"η out of bounds: {η}")

                entropy_arr[j] = dS
                free_energy_arr[j] = dF
                ηs_arr[j] = η

            except Exception as e:
                # print(f"\n\t\tError processing τ = {τ}: {e}")
                break

            j += 1

        # Save results as txt file
        np.savetxt(
            f"Results/efficiency/{id}_efficiency.txt",
            np.column_stack((τs, entropy_arr, free_energy_arr, ηs_arr)),
            header="tau\tentropy\tfree_energy\tnetwork_efficiency",
        )


def get_bipartite_sets(G):
    """
    Get the two node sets of a bipartite graph.

    First tries to use a node attribute defining the bipartition.
    If no such attribute is available, falls back to NetworkX bipartite coloring.

    Returns
    -------
    U, V : list
        The two node sets.
    """

    # Try common node attributes first
    possible_attrs = ["bipartite_set", "bipartite", "guild", "set", "partite", "level"]

    for attr in possible_attrs:
        values = nx.get_node_attributes(G, attr)

        if len(values) == G.number_of_nodes():
            unique_values = list(set(values.values()))

            if len(unique_values) == 2:
                U = [
                    node for node, value in values.items() if value == unique_values[0]
                ]
                V = [
                    node for node, value in values.items() if value == unique_values[1]
                ]

                return U, V

    # Fallback: infer bipartition from graph topology
    if not nx.is_bipartite(G):
        raise ValueError("Curveball null model requires a bipartite graph.")

    color = nx.algorithms.bipartite.color(G)

    U = [node for node, value in color.items() if value == 0]
    V = [node for node, value in color.items() if value == 1]

    return U, V


def er_null_model(G, seed=None):
    """
    Erdős--Rényi null model preserving the number of nodes and expected density.
    For bipartite networks, edges are allowed only between the two partitions.
    """
    rng = np.random.default_rng(seed)

    if nx.is_bipartite(G):
        U, V = get_bipartite_sets(G)

        n_u = len(U)
        n_v = len(V)
        E = G.number_of_edges()

        if n_u == 0 or n_v == 0:
            return nx.Graph()

        p = E / (n_u * n_v)

        G_er = nx.Graph()
        G_er.add_nodes_from(G.nodes(data=True))

        for u in U:
            for v in V:
                if rng.random() < p:
                    G_er.add_edge(u, v)

        G_er.graph.update(G.graph)
        return G_er

    else:
        N = G.number_of_nodes()
        E = G.number_of_edges()

        if N < 2:
            return G.copy()

        p = 2 * E / (N * (N - 1))

        G_er = nx.erdos_renyi_graph(N, p, seed=seed)

        # Relabel integer nodes back to original node labels
        mapping = dict(zip(G_er.nodes(), list(G.nodes())))
        G_er = nx.relabel_nodes(G_er, mapping)

        # Preserve node attributes
        for node, attrs in G.nodes(data=True):
            G_er.nodes[node].update(attrs)

        G_er.graph.update(G.graph)
        G_er.remove_edges_from(nx.selfloop_edges(G_er))

        return G_er


def bipartite_degree_preserving_rewire(G, n_swaps=None, seed=None, max_tries=None):
    """
    Degree-preserving bipartite rewiring using double-edge swaps.

    Starting from the empirical bipartite graph, repeatedly replaces
    (u1, v1), (u2, v2) with (u1, v2), (u2, v1), preserving all node degrees
    and keeping the graph simple and bipartite.
    """
    rng = np.random.default_rng(seed)

    U, V = get_bipartite_sets(G)
    U = set(U)
    V = set(V)

    G_new = G.copy()

    # Keep only correctly oriented bipartite edges
    edges = []
    for a, b in G_new.edges():
        if a in U and b in V:
            edges.append((a, b))
        elif a in V and b in U:
            edges.append((b, a))

    # Cannot perform a double-edge swap with fewer than two edges
    if len(edges) < 2:
        print(
            f"Skipping bipartite CM rewiring for small network: {G.graph.get('name', 'unknown')}"
        )
        return G_new

    if n_swaps is None:
        n_swaps = 10 * len(edges)

    if n_swaps <= 0:
        return G_new

    if max_tries is None:
        max_tries = 100 * n_swaps

    swaps = 0
    tries = 0

    while swaps < n_swaps and tries < max_tries:
        tries += 1

        idx1, idx2 = rng.choice(len(edges), size=2, replace=False)

        a, b = edges[idx1]
        c, d = edges[idx2]

        # Need four distinct nodes
        if len({a, b, c, d}) < 4:
            continue

        new_e1 = (a, d)
        new_e2 = (c, b)

        # Avoid duplicate edges
        if G_new.has_edge(*new_e1) or G_new.has_edge(*new_e2):
            continue

        G_new.remove_edge(a, b)
        G_new.remove_edge(c, d)
        G_new.add_edge(*new_e1)
        G_new.add_edge(*new_e2)

        # Update edge list
        edges[idx1] = new_e1
        edges[idx2] = new_e2

        swaps += 1

    return G_new


def cm_null_model(G, seed=None):
    """
    Degree-preserving null model.

    For bipartite networks, rewiring is restricted to edges between partitions,
    preserving the degree of every node exactly. For non-bipartite networks,
    uses a simple configuration-model approximation.
    """
    if nx.is_bipartite(G):
        return bipartite_degree_preserving_rewire(
            G, n_swaps=10 * G.number_of_edges(), seed=seed
        )

    degree_sequence = [d for _, d in G.degree()]

    G_cm = nx.configuration_model(degree_sequence, seed=seed)
    G_cm = nx.Graph(G_cm)
    G_cm.remove_edges_from(nx.selfloop_edges(G_cm))

    # Relabel integer nodes back to original node labels
    mapping = dict(zip(G_cm.nodes(), list(G.nodes())))
    G_cm = nx.relabel_nodes(G_cm, mapping)

    for node, attrs in G.nodes(data=True):
        if node in G_cm:
            G_cm.nodes[node].update(attrs)

    G_cm.graph.update(G.graph)

    return G_cm


def curveball_trade(row_sets, rng):
    """
    Perform one Curveball trade between two randomly selected rows.

    row_sets[i] is the set of column indices containing 1s in row i.
    The operation preserves all row and column sums.
    """

    if len(row_sets) < 2:
        return row_sets

    i, j = rng.choice(len(row_sets), size=2, replace=False)

    row_i = row_sets[i]
    row_j = row_sets[j]

    common = row_i & row_j
    only_i = row_i - row_j
    only_j = row_j - row_i

    pool = list(only_i | only_j)

    if len(pool) == 0:
        return row_sets

    rng.shuffle(pool)

    new_only_i = set(pool[: len(only_i)])
    new_only_j = set(pool[len(only_i) :])

    row_sets[i] = common | new_only_i
    row_sets[j] = common | new_only_j

    return row_sets


def curveball_null_model(G, n_trades=None, seed=None):
    """
    Generate a Curveball fixed-marginal null model for a bipartite network.

    This null model preserves:
    - bipartite structure
    - row sums
    - column sums
    - the degree of every node
    - the total number of interactions

    It is only defined for bipartite graphs.
    """

    if not nx.is_bipartite(G):
        raise ValueError("Curveball null model requires a bipartite graph.")

    rng = np.random.default_rng(seed)

    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))

    U, V = get_bipartite_sets(G)

    U = list(U)
    V = list(V)

    v_to_col = {v: j for j, v in enumerate(V)}

    # Build row representation of the bipartite incidence matrix
    row_sets = []

    for u in U:
        cols = set()

        for v in G.neighbors(u):
            if v in v_to_col:
                cols.add(v_to_col[v])

        row_sets.append(cols)

    if n_trades is None:
        # Conservative default. Increase if you want stronger mixing.
        n_trades = 10 * len(U)

    for _ in range(n_trades):
        row_sets = curveball_trade(row_sets, rng)

    # Reconstruct randomized graph
    G_cb = nx.Graph()

    # Preserve nodes and their attributes
    G_cb.add_nodes_from(G.nodes(data=True))

    for i, u in enumerate(U):
        for j in row_sets[i]:
            G_cb.add_edge(u, V[j])

    G_cb.graph.update(G.graph)

    return G_cb


def compute_network_efficiency_null_model(
    filenames, τs, model, num_randomizations=100, overwrite=False
):

    num_to_compute = num_randomizations

    if not os.path.exists("Results/efficiency/null_models/{}".format(model)):
        os.makedirs("Results/efficiency/null_models/{}".format(model))

    print(f"\nComputing null models for {model}...\n")

    for filename in filenames:

        G_empirical = nx.read_graphml(filename)

        if model == "CB" and not nx.is_bipartite(G_empirical):
            print(
                f"Skipping {filename}: Curveball null model requires a bipartite network."
            )
            continue

        # If file already exists, skip
        id = filename.split("/")[-1].split(".")[0]
        output_file = f"Results/efficiency/null_models/{model}/{id}_efficiency.txt"
        if os.path.exists(output_file) and not overwrite:
            print(f"File {output_file} already exists, skipping...")
            continue

        entropy_arr = np.zeros((num_to_compute, len(τs)))
        free_energy_arr = np.zeros((num_to_compute, len(τs)))
        ηs_arr = np.zeros((num_to_compute, len(τs)))

        # Print progress
        print(
            "\nComputing {}/{}".format(filenames.index(filename) + 1, len(filenames)),
        )

        for k in range(num_to_compute):

            print(f"\tRandomization {k+1}/{num_to_compute}", end="\r")

            id = filename.split("/")[-1].split(".")[0]

            G_empirical = nx.read_graphml(filename)

            if G_empirical.number_of_edges() == 0:
                print(f"\n\t\tGraph has no edges, skipping...")
                break
            elif G_empirical.number_of_nodes() < 2:
                print(f"\n\t\tGraph has less than 2 nodes, skipping...")
                break

            if model == "ER":
                G = er_null_model(G_empirical, seed=k)

            elif model == "CM":
                G = cm_null_model(G_empirical, seed=k)

            elif model == "CB":
                if not nx.is_bipartite(G_empirical):
                    print(
                        "\n\t\tGraph is not bipartite, skipping Curveball null model..."
                    )
                    break

                G = curveball_null_model(
                    G_empirical,
                    n_trades=None,
                    seed=k,
                )

            else:
                raise ValueError(f"Unknown null model: {model}")

            j = 0
            for τ in τs:

                try:

                    dS, dF, η = thermo_spectrum(G, τ)

                    if η < 0 or η > 1:
                        raise ValueError(f"η out of bounds: {η}")

                    entropy_arr[k, j] = dS
                    free_energy_arr[k, j] = dF
                    ηs_arr[k, j] = η

                except Exception as e:
                    # print(f"\n\t\tError processing τ = {τ}: {e}")
                    break

                j += 1

        mean_entropy = np.nanmean(entropy_arr, axis=0)
        mean_free_energy = np.nanmean(free_energy_arr, axis=0)
        mean_ηs = np.nanmean(ηs_arr, axis=0)

        std_entropy = np.nanstd(entropy_arr, axis=0)
        std_free_energy = np.nanstd(free_energy_arr, axis=0)
        std_ηs = np.nanstd(ηs_arr, axis=0)

        # Save results as txt file
        np.savetxt(
            f"Results/efficiency/null_models/{model}/{id}_efficiency.txt",
            np.column_stack(
                (
                    τs,
                    mean_entropy,
                    mean_free_energy,
                    mean_ηs,
                    std_entropy,
                    std_free_energy,
                    std_ηs,
                )
            ),
            header="tau\tentropy\tfree_energy\tnetwork_efficiency\tentropy_std\tfree_energy_std\tnetwork_efficiency_std",
        )


def zscore(empirical, null_ensemble):
    mu = np.nanmean(null_ensemble, axis=0)
    sigma = np.nanstd(null_ensemble, axis=0)
    return (empirical - mu) / sigma


def add_link(G, percentage=1):
    """
    Add a percentage of links to the graph G.
    """
    G_new = G.copy()

    N = G.number_of_nodes()
    E = G.number_of_edges()

    num_links_to_add = max(1, int((percentage / 100) * E))

    if nx.is_bipartite(G_new):
        U, V = get_bipartite_sets(G_new)
        possible_edges = {(u, v) for u in U for v in V if not G_new.has_edge(u, v)}
    else:
        possible_edges = set(nx.non_edges(G_new))

    if len(possible_edges) == 0:
        return G_new  # No edges can be added

    edges_to_add = np.random.choice(
        len(possible_edges),
        size=min(num_links_to_add, len(possible_edges)),
        replace=False,
    )

    for idx in edges_to_add:
        edge = list(possible_edges)[idx]
        G_new.add_edge(*edge)

    return G_new


def remove_link(G, percentage=1):
    """
    Remove a percentage of links from the graph G.
    """
    G_new = G.copy()

    E = G.number_of_edges()

    num_links_to_remove = max(1, int((percentage / 100) * E))

    existing_edges = list(G_new.edges())

    if len(existing_edges) == 0:
        return G_new  # No edges can be removed

    edges_to_remove = np.random.choice(
        len(existing_edges),
        size=min(num_links_to_remove, len(existing_edges)),
        replace=False,
    )

    for idx in edges_to_remove:
        edge = existing_edges[idx]
        G_new.remove_edge(*edge)

    return G_new


def compute_efficiency_with_perturbations(filename, τs, realizations, percentage=1):

    G = nx.read_graphml(filename)

    entropy_arr = np.nan * np.ones(len(τs))
    free_energy_arr = np.nan * np.ones(len(τs))
    ηs_arr = np.nan * np.ones(len(τs))

    j = 0
    for τ in τs:

        try:

            dS, dF, η = thermo_spectrum(G, τ)

            if η < 0 or η > 1:
                raise ValueError(f"η out of bounds: {η}")

            entropy_arr[j] = dS
            free_energy_arr[j] = dF
            ηs_arr[j] = η

        except Exception as e:
            # print(f"\n\t\tError processing τ = {τ}: {e}")
            break

        j += 1

    # Modified networks with 1% more and 1% less links
    ηs_arr_extra = np.zeros((realizations, len(τs)))
    ηs_arr_one_less = np.zeros((realizations, len(τs)))

    for r in range(realizations):

        print(f"\tRealization {r+1}/{realizations}", end="\r")

        G_with_extra_link = add_link(G, percentage=percentage)
        G_with_one_less_link = remove_link(G, percentage=percentage)

        j = 0
        for τ in τs:
            try:
                dS_extra, dF_extra, η_extra = thermo_spectrum(G_with_extra_link, τ)
                dS_one_less, dF_one_less, η_one_less = thermo_spectrum(
                    G_with_one_less_link, τ
                )

                ηs_arr_extra[r, j] += η_extra
                ηs_arr_one_less[r, j] += η_one_less

            except Exception as e:
                pass

            j += 1

    diff_extra = ηs_arr_extra - ηs_arr
    diff_one_less = ηs_arr_one_less - ηs_arr

    diff_extra_mean = np.nanmean(diff_extra, axis=0)
    diff_one_less_mean = np.nanmean(diff_one_less, axis=0)

    diff_extra_std_err = np.nanstd(diff_extra, axis=0) / np.sqrt(realizations)
    diff_one_less_std_err = np.nanstd(diff_one_less, axis=0) / np.sqrt(realizations)

    # Average relative differences across τs to get a single summary metric
    mean_diff_more_links = np.nanmean(diff_extra_mean)
    mean_diff_less_links = np.nanmean(diff_one_less_mean)

    return (
        diff_extra_mean,
        diff_one_less_mean,
        diff_extra_std_err,
        diff_one_less_std_err,
        mean_diff_more_links,
        mean_diff_less_links,
    )
