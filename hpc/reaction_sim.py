from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

import CBRdb


def ax_plot(fig: plt.Figure, ax: plt.Axes, xlab: str, ylab: str, xs: int = 14, ys: int = 14) -> None:
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=ys - 2, direction='in', length=6, width=2)
    ax.tick_params(axis='both', which='minor', labelsize=ys - 2, direction='in', length=4, width=2)
    ax.tick_params(axis='both', which='both', top=True, right=True)
    ax.set_xlabel(xlab, fontsize=xs)
    ax.set_ylabel(ylab, fontsize=ys)
    fig.tight_layout()
    return None


def plot_kde(
        data,
        bandwidth=None,
        grid_size=1000,
        y_scale='log',
        xlab="Value",
        ylab="Frequency",
        fig=None,
        ax=None,
        fig_size=(8, 5),
        fontsize=16,
):
    """
    Plots a Kernel Density Estimate (KDE) for the given data using Matplotlib.

    Parameters:
    data (array-like): The data for which the KDE is to be plotted.
    bandwidth (float or None, optional): The bandwidth of the KDE. If None, it is automatically determined. Defaults to None.
    grid_size (int, optional): The number of points in the grid for evaluating the KDE. Defaults to 1000.
    y_scale (str or None, optional): The scale for the y-axis (e.g., 'log'). If None, no scaling is applied. Defaults to 'log'.
    xlab (str, optional): Label for the x-axis. Defaults to "Value".
    ylab (str, optional): Label for the y-axis. Defaults to "Frequency".
    fig (matplotlib.figure.Figure or None, optional): The Matplotlib figure object. If None, a new figure is created. Defaults to None.
    ax (matplotlib.axes.Axes or None, optional): The Matplotlib Axes object. If None, a new Axes is created. Defaults to None.
    fig_size (tuple, optional): The size of the figure in inches (width, height). Defaults to (8, 5).
    fontsize (int, optional): Font size for axis labels. Defaults to 16.

    Returns:
    tuple: A tuple containing the Matplotlib figure and axis objects (fig, ax).
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=fig_size)
    data = np.asarray(data)
    n = len(data)

    # Fit KDE
    kde = gaussian_kde(data, bw_method=bandwidth)

    # Evaluate KDE on a grid
    x_min, x_max = data.min(), data.max()
    xs = np.linspace(x_min, x_max, grid_size)
    ys = kde(xs)

    # Convert density to expected counts
    dx = xs[1] - xs[0]
    counts = ys * n * dx

    ax.set_xlim(x_min, x_max)
    if y_scale is not None:
        # Find the nearest order of magnitude to the maximum count
        order_of_magnitude = 10 ** np.floor(np.log10(max(counts)))
        # Set the y-axis limit to the next order of magnitude
        ax.set_ylim(1, order_of_magnitude * 10)
        ax.set_yscale('log')

    # Overlay KDE scaled to counts
    ax.plot(xs, counts, lw=2)
    ax_plot(fig, ax, xlab=xlab, ylab=ylab, xs=fontsize, ys=fontsize)

    return fig, ax


def find_similar_reactions(query_reaction, reaction_list):
    func_sim = partial(CBRdb.get_reaction_similarity, query_reaction)
    similarities = CBRdb.mp_calc(func_sim, reaction_list)
    # Sort by similarity score in descending order
    similarities.sort(key=lambda x: x, reverse=True)
    return similarities


def find_minmax_similar_reactions(query_reaction, reaction_list):
    func_sim = partial(CBRdb.get_reaction_similarity, query_reaction)
    similarities = CBRdb.mp_calc(func_sim, reaction_list)
    return min(similarities), max(similarities)


if __name__ == "__main__":
    print(flush=True)
    data_c = pd.read_csv('../CBRdb_C.csv', low_memory=False)
    data_r = pd.read_csv('../CBRdb_R.csv', low_memory=False)

    print("Generating SMARTS for reactions...")
    # if the smarts column does not exist, create it
    if 'smarts' not in data_r.columns:
        func_smarts = partial(CBRdb.to_smarts_rxn_line, data_c=data_c, add_stoich=False)
        eq_list = data_r['reaction'].tolist()
        data_r['smarts'] = CBRdb.mp_calc(func_smarts, eq_list)

    # Only select the reactions where the id starts with 'R'
    sel = data_r['id'].str.startswith('R')
    # Split the data into KEGG and ATLAS reactions
    data_r_atlas = data_r[~sel]
    data_r_kegg = data_r[sel]

    # print the number of reactions
    print(f"Number of reactions (KEGG) : {len(data_r_kegg)}")
    print(f"Number of reactions (ATLAS): {len(data_r_atlas)}")

    smarts_list_kegg = data_r_kegg['smarts'].tolist()
    smarts_list_atlas = data_r_atlas['smarts'].tolist()

    n = 100
    sim_min = np.zeros(n)
    sim_max = np.zeros(n)
    for i in range(n):
        print(f"Processing reaction {i + 1}/{n}")
        sim_min[i], sim_max[i] = find_minmax_similar_reactions(smarts_list_atlas[i], smarts_list_kegg)

    plt.hist(sim_max, bins=50, alpha=0.5)
    plt.xlabel('Similarity score')
    plt.ylabel('Frequency')
    plt.show()

    # smarts_list = data_r_atlas['smarts'].tolist()
    # scores_atlas = find_similar_reactions(smarts_0, smarts_list)
    #
    # print(f"Max similarity score (KEGG) : {max(scores_kegg)}")
    # print(f"Max similarity score (ATLAS): {max(scores_atlas)}")
    # print(f"Min similarity score (KEGG) : {min(scores_kegg)}")
    # print(f"Min similarity score (ATLAS): {min(scores_atlas)}")
    #
    # plt.hist(scores_kegg, bins=50, alpha=0.5, label='KEGG')
    # plt.legend()
    # plt.xlabel('Similarity score')
    # plt.ylabel('Frequency')
    # plt.title('Histogram of reaction similarity scores')
    # plt.show()
    #
    # plt.hist(scores_atlas, bins=50, alpha=0.5, label='ATLAS')
    # plt.legend()
    # plt.xlabel('Similarity score')
    # plt.ylabel('Frequency')
    # plt.title('Histogram of reaction similarity scores')
    # plt.show()
    #
    # plot_kde(scores_kegg, y_scale=None)
    # plt.show()
    #
    # plot_kde(scores_atlas, y_scale=None)
    # plt.show()
