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


def plot_histogram(data,
                   bins=30,
                   xlab='Best Similarity score',
                   ylab='Frequency',
                   figsize=(8, 5),
                   fontsize=16,
                   ):
    """
    Plots a histogram for the given data.

    Parameters:
    data (list or array-like): The data to be plotted in the histogram.
    bins (int, optional): Number of bins for the histogram. Default is 30.
    xlab (str, optional): Label for the x-axis. Default is 'Values'.
    ylab (str, optional): Label for the y-axis. Default is 'Frequency'.
    figsize (tuple, optional): Size of the figure in inches (width, height). Default is (8, 6).
    fontsize (int, optional): Font size for axis labels. Default is 16.

    Returns:
    tuple: A tuple containing the Matplotlib figure and axis objects (fig, ax).
    """
    fig, ax = plt.subplots(figsize=figsize)
    plt.hist(data, bins=bins, color='blue', edgecolor='black', alpha=0.8)
    ax_plot(fig, ax, xlab=xlab, ylab=ylab, xs=fontsize, ys=fontsize)
    return fig, ax


if __name__ == "__main__":
    print(flush=True)
    data_c = pd.read_csv('../CBRdb_C.csv', low_memory=False)
    data_r = pd.read_csv('../CBRdb_R.csv', low_memory=False)

    # Trim them down to only the necessary columns
    data_r = pd.DataFrame(data_r[['id', 'reaction']])

    # if the smarts column does not exist, create it
    if 'smarts' not in data_r.columns:
        print("Generating SMARTS for reactions...", flush=True)
        func_smarts = partial(CBRdb.to_smarts_rxn_line, data_c=data_c, add_stoich=False)
        data_r['smarts'] = CBRdb.mp_calc(func_smarts, data_r['reaction'].tolist())

    print("Generating fingerprints for reactions...", flush=True)
    data_r['fp_struct'] = CBRdb.mp_calc(CBRdb.get_rxn_fingerprint, data_r['smarts'].tolist())

    # Only select the reactions where the id starts with 'R'
    sel = data_r['id'].str.startswith('R')
    # Split the data into KEGG and ATLAS reactions
    data_r_atlas = data_r[~sel]
    data_r_kegg = data_r[sel]

    # print the number of reactions
    print(f"Number of reactions (KEGG) : {len(data_r_kegg)}")
    print(f"Number of reactions (ATLAS): {len(data_r_atlas)}")

    fp_kegg = data_r_kegg['fp_struct'].tolist()
    fp_atlas = data_r_atlas['fp_struct'].tolist()

    func_sim = partial(CBRdb.find_max_similar_rxn, rxn_db=fp_kegg)
    sim_max, sim_max_idx = zip(*CBRdb.mp_calc(func_sim, fp_atlas))

    # convert the idx into the corresponding reaction id
    sim_max_id = [data_r_kegg.iloc[idx]['id'] for idx in sim_max_idx]

    # add the results to the atlas dataframe
    data_r_atlas['sim_max'] = sim_max
    data_r_atlas['sim_max_id'] = sim_max_id

    # Print the lowest 10 similarity scores
    print("Lowest 10 similarity scores:", flush=True)
    print(data_r_atlas.nsmallest(10, 'sim_max')[['id', 'sim_max', 'sim_max_id']], flush=True)

    plot_histogram(sim_max)
    plt.savefig('ATLAS_R_sim_hist.png', dpi=300)
    plt.savefig('ATLAS_R_sim_hist.pdf', dpi=300)
    plt.show()

    plot_histogram(sim_max)
    plt.yscale('log')
    plt.savefig('ATLAS_R_sim_log_hist.png', dpi=300)
    plt.savefig('ATLAS_R_sim_log_hist.pdf', dpi=300)
    plt.show()

    data_r_atlas.to_csv('ATLAS_R_sim.csv.gz', index=False, compression='gzip')
    print('Done.', flush=True)
    print(flush=True)
