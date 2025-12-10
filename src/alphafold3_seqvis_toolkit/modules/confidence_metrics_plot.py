# this module provides a simple overview of all confidence measures in an AlphaFold 3 prediction.

import json
import os 
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from typing import Dict, Optional, Tuple
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap


# function to load JSON data from a file
def load_json_data(json_file_path):
    try:
        with open(json_file_path, "r") as f:
            data = json.load(f)
        return data
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return None



# 1, For Global confidence measures and ipTM matrix
def plot_global_confidence(confid_json_file_path, output_path):
    """
    Description
    -----------
    plot global confidence metrics, such as ipTM matrix and chain-wise ipTM/pTM scores.

    Args
    ----
    confid_json_file_path : str
        path to the JSON file containing global confidence metrics.
    output_path : str
        path to save the output plots and data files.

    Returns
    -------
    output a text file of scalar values,
    a tsv file of list values,
    and several pdf files of plots for list and 2D ndarray values.
        
    Notes
    -----
    - 1, Currently, the function of customizing output file names is not yet supported. 
    All generated figures will be saved in a single dedicated folder, 
    and the folder will be named by the system for your convenience.
    """

    global_confidence = load_json_data(confid_json_file_path)

    # extract global confidence measures
    chain_iptm = np.asarray(global_confidence['chain_iptm'], dtype=float) # list
    chain_pair_iptm = np.asarray(global_confidence['chain_pair_iptm'], dtype=float) # 2D ndarray 
    chain_pair_pae_min = np.asarray(global_confidence['chain_pair_pae_min'], dtype=float) # 2D ndarray
    chain_ptm = np.asarray(global_confidence['chain_ptm'], dtype=float) # list
    fraction_disordered = global_confidence['fraction_disordered'] # scalar 
    has_clash = global_confidence['has_clash'] # boolean
    iptm = global_confidence['iptm'] # scalar
    num_recycles = global_confidence['num_recycles'] # scalar
    ptm = global_confidence['ptm'] # scalar
    ranking_score = global_confidence['ranking_score'] # scalar

    # (1) For Scalar value, we directly print them and write them in a text file
    # Maybe we can plot them in structure viewer
    # print(f"Fraction Disordered: {fraction_disordered}")
    # print(f"Has Clash: {has_clash}")
    # print(f"ipTM: {iptm}")
    # print(f"Number of Recycles: {num_recycles}")
    # print(f"pTM: {ptm}")
    # print(f"Ranking Score: {ranking_score}")
    
    # write them to a tsv file
    job_name = re.match(r'fold_(.*)_summary_confidences_\d+\.json', os.path.basename(confid_json_file_path)).group(1) 
    with open(f"{output_path}/{job_name}_global_confidence_SCALAR_measures.tsv", "w") as f:
        f.write(f"Fraction Disordered\t{fraction_disordered}\n")
        f.write(f"Has Clash\t{has_clash}\n")
        f.write(f"ipTM\t{iptm}\n")
        f.write(f"Number of Recycles\t{num_recycles}\n")
        f.write(f"pTM\t{ptm}\n")
        f.write(f"Ranking Score\t{ranking_score}\n")

    # nunmerical chain index to chain id conversion
    def number_to_chain_idx(num):
        """Convert a chain index to a chain ID (A, B, C, ..., Z, AA, AB, ...)."""
        letter = ""
        while num >= 0: 
            letter = chr(num % 26 + ord('A')) + letter
            num = num // 26 - 1
        return letter


    # (2) For list value, we output them as a dataframe, and plot them as a line plot or bar plot (chain id as x axis)
    # for chain_iptm or chain_ptm, we can plot them as line plot and bar plot
    x = np.arange(len(chain_iptm))
    # first line plot
    plt.figure(figsize=(10,5))
    plt.plot(chain_iptm, marker='o', color='b', alpha=0.7, label='Chain ipTM')
    plt.plot(chain_ptm, marker='o', color='r', alpha=0.7, label='Chain pTM')
    plt.title("Chain ipTM/pTM Scores")
    plt.xlabel("Chain Index")
    plt.ylabel("ipTM/pTM Score")
    plt.xticks(x, [number_to_chain_idx(i) for i in x])
    plt.legend()
    plt.grid()
    # we save the figure into a pdf file
    plt.savefig(f"{output_path}/{job_name}_global_confidence_chain_ptm_iptm_lineplot.pdf", bbox_inches='tight')
    # also save as png, dpi 300
    plt.savefig(f"{output_path}/{job_name}_global_confidence_chain_ptm_iptm_lineplot.png", bbox_inches='tight', dpi=300)
    plt.close()

    # and bar plot
    width = 0.35
    plt.figure(figsize=(10,5))
    plt.bar(x - width/2, chain_iptm, width, label='Chain ipTM', color='b', alpha=0.7)
    plt.bar(x + width/2, chain_ptm, width, label='Chain pTM', color='r', alpha=0.7)
    plt.title("Chain ipTM/pTM Scores")
    plt.xlabel("Chain Index")
    plt.ylabel("ipTM/pTM Score")
    plt.xticks(x, [number_to_chain_idx(i) for i in x])
    plt.legend()
    plt.grid()
    # we save the figure into a pdf file
    plt.savefig(f"{output_path}/{job_name}_global_confidence_chain_ptm_iptm_barplot.pdf", bbox_inches='tight')
    # also save as png, dpi 300
    plt.savefig(f"{output_path}/{job_name}_global_confidence_chain_ptm_iptm_barplot.png", bbox_inches='tight', dpi=300)
    plt.close()


    # write them to a tsv file
    with open(f"{output_path}/{job_name}_global_confidence_chain_ptm_iptm.tsv", "w") as f:
        f.write("Chain_Index\tChain_ipTM_Score\tChain_pTM_Score\n")
        for i in range(len(chain_iptm)):
            f.write(f"{number_to_chain_idx(i)}\t{chain_iptm[i]}\t{chain_ptm[i]}\n")


    # (3) For 2D ndarray value, we output them as heatmap plot (chain id as x and y axis)
    import seaborn as sns
    # for chain_pair_iptm
    plt.figure(figsize=(8,6))
    sns.heatmap(chain_pair_iptm, annot=True, fmt=".2f", cmap="coolwarm",
                xticklabels=[number_to_chain_idx(i) for i in range(chain_pair_iptm.shape[0])],
                yticklabels=[number_to_chain_idx(i) for i in range(chain_pair_iptm.shape[1])])
    plt.title("Chain Pair ipTM Scores Heatmap")
    plt.xlabel("Chain Index")
    plt.ylabel("Chain Index")
    # we save the figure into a pdf file
    plt.savefig(f"{output_path}/{job_name}_global_confidence_chain_pair_iptm_heatmap.pdf", bbox_inches='tight')
    # also save as png, dpi 300
    plt.savefig(f"{output_path}/{job_name}_global_confidence_chain_pair_iptm_heatmap.png", bbox_inches='tight', dpi=300)
    plt.close()

    # for chain_pair_pae_min
    plt.figure(figsize=(8,6))
    sns.heatmap(chain_pair_pae_min, annot=True, fmt=".2f", cmap="coolwarm",
                xticklabels=[number_to_chain_idx(i) for i in range(chain_pair_pae_min.shape[0])],
                yticklabels=[number_to_chain_idx(i) for i in range(chain_pair_pae_min.shape[1])])
    plt.title("Chain Pair PAE Min Heatmap")
    plt.xlabel("Chain Index")
    plt.ylabel("Chain Index")
    # we save the figure into a pdf file
    plt.savefig(f"{output_path}/{job_name}_global_confidence_chain_pair_pae_min_heatmap.pdf", bbox_inches='tight')
    # also save as png, dpi 300
    plt.savefig(f"{output_path}/{job_name}_global_confidence_chain_pair_pae_min_heatmap.png", bbox_inches='tight', dpi=300)
    plt.close()

    

# 2, For local confidence measures
def plot_local_confidence(full_json_file_path, output_path, chains: Optional[object]=None, tick_step: int = 100):
    """
    Description
    -----------
    plot local confidence metrics, such as PAE matrix and per-atom pLDDT scores.
    
    Args
    ----
    full_json_file_path : str
        path to the JSON file containing local confidence metrics.
    output_path : str
        path to save the output plots and data files.
    chains : str or list or tuple, optional, currently designed for PAE matrix plot only
        specify the chain id(s) to plot. If None, plot all chains. Default is None.
        - str: single chain id, e.g., 'A'
        - list or tuple: multiple chain ids, e.g., ['A', 'B']
    tick_step : int, optional
        step size for residue ticks on axes. Default is 100.

    Returns
    -------
    output a figure of PAE matrix for the specified chains.
    
    """

    # load local confidence data
    local_confidence = load_json_data(full_json_file_path)
    if local_confidence is None:
        raise ValueError("Failed to load local confidence data from JSON file.")
    
    # extract local confidence measures
    pae_matrix = np.asarray(local_confidence['pae'], dtype=float)  # 2D ndarray
    contact_probs = np.asarray(local_confidence['contact_probs'], dtype=float)  # 2D ndarray
    atom_chain_ids = np.asarray(local_confidence['atom_chain_ids'], dtype=str)  # list of str   
    atom_plddts = np.asarray(local_confidence['atom_plddts'], dtype=float)  # list of float
    token_chain_ids = np.asarray(local_confidence['token_chain_ids'], dtype=str)  # list of str
    token_res_ids = np.asarray(local_confidence['token_res_ids'], dtype=int)  # list of int


    # convert chains param into a hashable list
    if chains is None:
        selected_chains = None
    elif isinstance(chains, str):
        # single chain id as string, like 'A'
        selected_chains = [chains]
    else:
        # list or tuple of chain ids, like ['A','B'] or ('A','B')
        selected_chains = list(chains)

    # select token indices for the specified chains
    if selected_chains is not None:
        missing = [chain_id for chain_id in selected_chains if chain_id not in np.unique(token_chain_ids)]
        if missing:
            raise ValueError(f"Specified chains not found in data: {missing}")
        
        # create a boolean mask for selected chains, res in selected chains will be True
        # token/res level
        mask = np.isin(token_chain_ids, selected_chains)
        idx = np.where(mask)[0] # tuple to index array

        if idx.size == 0:
            raise ValueError("No residues found for the specified chains.")
        
        # filter pae_matrix and contact_probs based on selected chains, and other token-based arrays
        # Note: we only filter token-based arrays here, atom-based arrays will be filtered later
        pae_matrix_sub = pae_matrix[np.ix_(idx, idx)]
        contact_probs_sub = contact_probs[np.ix_(idx, idx)]
        token_chain_ids_sub = token_chain_ids[idx]
        token_res_ids_sub = token_res_ids[idx]

        # atom level
        atom_mask = np.isin(atom_chain_ids, selected_chains)
        atom_chain_ids_sub = atom_chain_ids[atom_mask]
        atom_plddts_sub = atom_plddts[atom_mask]

        if atom_plddts_sub.size == 0:
            raise ValueError("No atoms found for the specified chains.")

    # if no chains specified, use all data (default)
    else:
        # token/res level
        pae_matrix_sub = pae_matrix
        contact_probs_sub = contact_probs
        token_chain_ids_sub = token_chain_ids
        token_res_ids_sub = token_res_ids
        idx = np.arange(pae_matrix.shape[0])

        # atom level
        atom_chain_ids_sub = atom_chain_ids
        atom_plddts_sub = atom_plddts
        
    # compute tick positions and labels for the selected chains
    xticks_loc = []
    xticks_labels = []
    for index, res in enumerate(token_res_ids_sub):
        # Note that index is the index in the sub-matrix, while res is the residue id
        # for example, index 728 in submatrix may correspond to residue 1 in chain B
        
        # ⚠️ Note that res is 1-based residue id, so we convert it to 0-based for tick logic
        res_0b = res - 1
        if res == 1 or (tick_step and res_0b % tick_step == 0):
            xticks_loc.append(index)
            xticks_labels.append(int(res_0b))

    job_name = re.match(r'fold_(.*)_full_data_\d+\.json', os.path.basename(full_json_file_path)).group(1) 

    # Cause PAE and contact_probs matrices are similar in plotting style and chain_id segmentation need, so we define a general plotting function here.
    def _plot_matrix_with_chain_bars(matrix, title, xlabel, ylabel, cbar_label,cmap_name, outfilename):
        fig, ax = plt.subplots(figsize=(10,10))
        im = ax.imshow(matrix, cmap=cmap_name, origin="upper", aspect="auto")
        ax.set_title(title,fontsize=10)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_xticks(xticks_loc)
        ax.set_xticklabels(xticks_labels, fontsize=10)
        # ax.set_yticks(xticks_loc)
        # ax.set_yticklabels(xticks_labels, fontsize=10)
        ax.set_yticks([])

        # add colorbar
        cbar = fig.colorbar(im, ax=ax, fraction=0.045)
        cbar.set_label(cbar_label, fontsize=10)

        # if multiple chains present (either full file with >1 chain, or selected_chains with >1 chain)
        # draw small colored bars on top and left showing chain segmentation
        unique_chains = np.unique(token_chain_ids_sub)
        draw_bars = len(unique_chains) > 1
        # if chain left > 1, we will draw bars and segmentation lines
        if draw_bars:
            divider = make_axes_locatable(ax)
            ax_top = divider.append_axes("top", size="5%", pad=0.03)
            ax_left = divider.append_axes("left", size="5%", pad=0.03)

            # build contiguous blocks (chain_id, start, end) relative to token_chain_ids_sub order
            chain_blocks = []
            start = 0
            for i in range(1, len(token_chain_ids_sub) + 1):
                if i == len(token_chain_ids_sub) or token_chain_ids_sub[i] != token_chain_ids_sub[start]:
                    chain_blocks.append((token_chain_ids_sub[start], start, i - 1))
                    start = i
            
            # create a color map for chains, map chain ids -> integers for colors
            chain_to_int = {chain_id: i for i, chain_id in enumerate(unique_chains)}
            chain_row = np.asarray([chain_to_int[chain_id] for chain_id in token_chain_ids_sub]).reshape(1, -1)
    
            # use pastel colors for the top/left chain bars and slightly transparent
            # 20 coloars should be enough for most cases!  
            if len(unique_chains) <= 20:
                cmap = plt.get_cmap("tab20", len(unique_chains))
            else:
            # Use gist_rainbow for >20 chains to ensure distinctness
            # set 0.1 to 0.9 to avoid too light/dark colors            
                colors = plt.get_cmap("gist_rainbow")(np.linspace(0.1, 0.9, len(unique_chains)))
                cmap = ListedColormap(colors)

            ax_top.imshow(chain_row, cmap=cmap, aspect="auto", alpha=0.9)
            ax_top.set_xticks([])
            ax_top.set_yticks([])
            ax_left.imshow(chain_row.T, cmap=cmap, aspect="auto", alpha=0.9)
            ax_left.set_xticks([])
            ax_left.set_yticks([])

            # draw dashed separations on main axes at block boundaries (skip boudarys at 0 and N)
            for _, s, e in chain_blocks:
                if s != 0:
                    sep = s - 0.5
                    ax.axvline(sep, color="k", linewidth=1, linestyle="--")
                    ax.axhline(sep, color="k", linewidth=1, linestyle="--")
                    ax_top.axvline(sep, color="w", linewidth=1)
                    ax_left.axhline(sep, color="w", linewidth=1)

            # annotate chain ids centered on their contiguous blocks (bold)
            for cid, s, e in chain_blocks:
                center = (s + e) / 2.0
                ax_top.text(center, 0, str(cid), ha="center", va="center",
                            fontsize=14, weight="bold", color="#222222")
                ax_left.text(0, center, str(cid), ha="center", va="center",
                             fontsize=14, weight="bold", color="#222222", rotation=90)
                

            # save figure
            sel_name = "_".join(selected_chains) if selected_chains else "all"
            plt.savefig(f"{output_path}/{job_name}_{outfilename}_{sel_name}.pdf", bbox_inches='tight')
            # also save as png, dpi 300
            plt.savefig(f"{output_path}/{job_name}_{outfilename}_{sel_name}.png", bbox_inches='tight', dpi=300)
            plt.close(fig)


    # Now, we can directly use the above function to plot PAE matrix and contact probability matrix
    
    # 1, Plot PAE matrix for specified chains
    _plot_matrix_with_chain_bars(
        matrix=pae_matrix_sub,
        title="Predicted Aligned Error (PAE) Matrix",
        xlabel="Scored Residue",
        ylabel="Aligned Residue",
        cbar_label="Expected Position Error (Å)",
        cmap_name="Greens_r",
        outfilename="local_confidence_PAE_matrix_selected_chains"
    )

    # 2, Plot contact probability matrix for specified chains
    _plot_matrix_with_chain_bars(
        matrix=contact_probs_sub,
        title="Contact Probability Matrix",
        xlabel="Residue i",
        ylabel="Residue j",
        cbar_label="Contact Probability (≤8Å)",
        cmap_name="viridis",  # after test, we found that viridis is better than RbBu_r
        outfilename="local_confidence_contact_probability_matrix_selected_chains"
    )


   # for atom pLDDT scores, we now plot them
    def _plot_atom_plddt_chain_bar(atom_plddt, atom_chain_ids_arr, outfilename, title="Atom pLDDT Track"):
        """
        1D pLDDT line plot with colored fill-under based on pLDDT thresholds,
        plus a top chain color bar and dashed separators at chain boundaries.
        """
        n = len(atom_plddt)
        fig, ax = plt.subplots(figsize=(20,10))
        x = np.arange(n)
        ax.plot(x, atom_plddt, color='black', linewidth=0.8)
        ax.fill_between(x, atom_plddt, 0, color='skyblue', alpha=0.2)
        ax.set_xlabel("Atom Index (token order)", fontsize=15)
        ax.set_ylabel("pLDDT Score", fontsize=15)
        ax.set_title(title, fontsize=15)


        # pLDDT threshold colors 
        color_vh = (0.051, 0.341, 0.827)   # >90
        color_h  = (0.416, 0.796, 0.945)   # 90 > x > 70
        color_l  = (0.996, 0.851, 0.212)   # 70 > x > 50
        color_vl = (0.992, 0.490, 0.302)   # <=50

        # masks for thresholds (per-atom)
        p = np.asarray(atom_plddt)
        mask_vh = p > 90
        mask_h  = (p <= 90) & (p > 70)
        mask_l  = (p <= 70) & (p > 50)
        mask_vl = p <= 50

        # fill under curve with different colors for each category
        ax.fill_between(x, p, 0, where=mask_vh, facecolor=color_vh, interpolate=True, alpha=0.7)
        ax.fill_between(x, p, 0, where=mask_h,  facecolor=color_h,  interpolate=True, alpha=0.7)
        ax.fill_between(x, p, 0, where=mask_l,  facecolor=color_l,  interpolate=True, alpha=0.7)
        ax.fill_between(x, p, 0, where=mask_vl, facecolor=color_vl, interpolate=True, alpha=0.7)

        # small legend like the reference image
        patches = [
            mpatches.Patch(color=color_vh, label="Very high (pLDDT > 90)"),
            mpatches.Patch(color=color_h,  label="High (90 ≥ pLDDT > 70)"),
            mpatches.Patch(color=color_l,  label="Low (70 ≥ pLDDT > 50)"),
            mpatches.Patch(color=color_vl, label="Very low (pLDDT ≤ 50)"),
        ]
        # place legend to the right of the axes, centered vertically and just outside the plot
        ax.legend(
            handles=patches,
            bbox_to_anchor=(1.00, 0.5),  # x,y: just outside right edge, vertically centered
            loc="center left",
            borderaxespad=0.0,
            fontsize=10,
            frameon=False,
            title="Model Confidence",
        )

        # compute minimal xticks (choose ~20 ticks)
        if n > 0:
            step = max(1, n // 20)
            ax.set_xticks(np.arange(0, n, step))
            ax.set_xlim(0, n - 1)

        # build contiguous chain blocks for atoms
        chain_blocks_atoms = []
        start = 0
        for i in range(1, n + 1):
            if i == n or atom_chain_ids_arr[i] != atom_chain_ids_arr[start]:
                chain_blocks_atoms.append((atom_chain_ids_arr[start], start, i - 1))
                start = i

        unique_chains_atoms = np.unique(atom_chain_ids_arr)
        if len(unique_chains_atoms) > 1:
            # map chain -> integer color index
            chain_to_int_atoms = {cid: i for i, cid in enumerate(unique_chains_atoms)}
            chain_row = np.asarray([chain_to_int_atoms[c] for c in atom_chain_ids_arr]).reshape(1, -1)

        # Automatically generate colors based on number of chains
        if len(unique_chains_atoms) <= 20:
            cmap = plt.get_cmap("tab20", len(unique_chains_atoms))
        else:
            # use gist_rainbow for >20 chains to ensure distinctness
            # set 0.1 to 0.9 to avoid too light/dark colors
            colors = plt.get_cmap("gist_rainbow")(np.linspace(0.1, 0.9, len(unique_chains_atoms)))
            cmap = ListedColormap(colors)

        # top color bar
        divider = make_axes_locatable(ax)
        ax_top = divider.append_axes("top", size="8%", pad=0.03)
        ax_top.imshow(chain_row, cmap=cmap, aspect="auto", alpha=0.95)
        ax_top.set_xticks([])
        ax_top.set_yticks([])

        # draw dashed separations on main axis and top bar at block boundaries
        for _, s, e in chain_blocks_atoms:
            if s != 0:
                sep = s - 0.5
                ax.axvline(sep, color="k", linewidth=1, linestyle="--")
                ax_top.axvline(sep, color="w", linewidth=1)

        # annotate chain ids centered on blocks
        for cid, s, e in chain_blocks_atoms:
            center = (s + e) / 2.0
            ax_top.text(center, 0, str(cid), ha="center", va="center",
                        fontsize=15, weight="bold", color="#222222")
            
        # save figure
        sel_name = "_".join(selected_chains) if selected_chains else "all"
        plt.savefig(f"{output_path}/{job_name}_{outfilename}_{sel_name}.pdf", bbox_inches='tight')
        # also save as png, dpi 300
        plt.savefig(f"{output_path}/{job_name}_{outfilename}_{sel_name}.png", bbox_inches='tight', dpi=300)
        plt.close(fig)
    
    # 3, now we can use the above function to plot atom pLDDT scores
    _plot_atom_plddt_chain_bar(
        atom_plddt=atom_plddts_sub,
        atom_chain_ids_arr=atom_chain_ids_sub,
        outfilename="local_confidence_atom_plddt_selected_chains",
        title="Atom pLDDT Scores for Selected Chains"
    )

    # 4, Statistics of pLDDT scores for the selected chains AND all
    # Note this is not chain specific output, unrelevant to selected chains
    def _plddt_statistics(atom_plddt, outfilename):
        """
        Compute and save pLDDT statistics to a text file.
        Statistics include mean, median, std, and confidence category fractions.
        """

        # first we define a function to compute statistics
        def _stats_from_vals(vals):
            if vals.size == 0:
                return (np.nan, np.nan, np.nan, 0.0, 0.0, 0.0, 0.0)
            mean = float(np.nanmean(vals))
            median = float(np.nanmedian(vals))
            std = float(np.nanstd(vals))
            n = vals.size
            vh = float(np.sum(vals > 90)) / n
            h  = float(np.sum((vals <= 90) & (vals > 70))) / n
            l  = float(np.sum((vals <= 70) & (vals > 50))) / n
            vl = float(np.sum(vals <= 50)) / n
            return (mean, median, std, vh, h, l, vl)
        
        rows = []
        # overall (all atoms)
        rows.append(("All",) + _stats_from_vals(atom_plddt))

        unique_chains = np.unique(atom_chain_ids)
        for cid in unique_chains:
            mask = (atom_chain_ids == cid)
            vals = atom_plddt[mask]
            rows.append((str(cid),) + _stats_from_vals(vals))
        
        # write to tsv file
        with open(f"{output_path}/{job_name}_{outfilename}_plddt_statistics.tsv", "w") as f:
            f.write("Chain_ID\tMean_pLDDT\tMedian_pLDDT\tStd_pLDDT\tFraction_Very_High(>90)\tFraction_High(90-70)\tFraction_Low(70-50)\tFraction_Very_Low(<=50)\n")
            for row in rows:
                f.write("\t".join([f"{x:.2f}" if isinstance(x, float) else str(x) for x in row]) + "\n")
    # compute and save pLDDT statistics
    _plddt_statistics(
        atom_plddt=atom_plddts,
        outfilename="local_confidence_overall_atom"
    )





