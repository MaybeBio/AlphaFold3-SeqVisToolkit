import os
import re
import json
import numpy as np
from typing import Optional, List, Union, Dict, Any
from Bio.PDB import MMCIFParser, Polypeptide
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

# from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable 
# fail at multi-track plotting, that's why we import gridspec below, but can be imported again to draw Scale bar

import matplotlib.gridspec as gridspec
from alphafold3_seqvis_toolkit.utils.track_utils import parse_bed_to_track_data

def contact_map_vis_with_track(
    mmcif_file: str,
    out_path: str,
    chains: Optional[Union[str, List[str]]] = None,
    cmap: str = "RdBu", # Set RdBu so Red is close (contact), Blue is far, or coolwarm_r
    track_bed_file: Optional[str] = None,
    color_config: Optional[str] = "tab10",
    tick_step: int = 100
):
    """
    Description
    -----------
    Generate a contact map (distance matrix) for a single structure from an mmCIF file.
    Uses representative atoms: CA for Protein, C1' for Nucleic Acids, First Atom for others.
    
    Args
    -----
        mmcif_file (str): Path to the mmcif file.
        out_path (str): Path to save the output plot (e.g., "/data2").
        chains (str | List[str], optional): Chain ID(s) to include. 
                                            e.g., "A" or ["A", "B"]. 
                                            If None, all chains are included.
        vmax (float): Maximum distance (Angstrom) for colorbar scaling. Default is 95th percentile of distances.
        cmap (str): Colormap to use. Default is "RdBu_r" (Red=Close, Blue=Far).
        track_bed_file (str, optional): Path to the BED-like file containing 1D track data for overlay.see utils/track_utils.py for details.
        color_config (Union[str, Dict[str, Any]]): Color configuration for the 1D tracks/Path to color config file (JSON) or colormap name, generally will be parsed into Union[str, Dict[str, Any]] format, and kind of like below:
            - String: A colormap name (e.g. "tab10") or a single color (e.g. "orange").
            - Dict: {TrackName: ColorConfig}.
              e.g. {"IDR": "red", "Domain": {"DomainA": "blue", "DomainB": "green"}}
              or {"IDR": "red", "Domain": "tab10"}
        Anyway, we will input something whose format like above to parameter color in parse_bed_to_track_data() function in track_utils.py.
        tick_step (int): Step size for ticks on the axes. Default is 100.
    
    Notes
    ------
    - 1, We currently support one AND list additional track_data overlay.
    - 2, The 1D track data should be prepared externally in some kind certain fomrat, see utils/track_utils.py for details.
    - 3, When not providing 1D track data, we will not plot any overlays.
    - 4, Track data is parsed using parse_bed_to_track_data function from track_utils.py, which needs parameter color_config and track_bed_file.
    The logic is simple: track_bed_file + color_config -> parse_bed_to_track_data() -> track_data for plotting.
    
    parse_bed_to_track_data() will return something like below:

    track_data (Dict[str, Any], optional): Additional data (1D feature track data) to overlay on the contact map.
            Format:
            [{
                "track_name": "IDR",
                "track_type": "categiorical" or "numerical",
                "color": "red" or {"A": "red", "B": "blue"} or "tab10",
                "track_data": {"A": [0.1, ...], "B": [0.5, ...]}, # Lists must match residue count
            }]
    - 5, The track bed file must be 0-based indexed!
    - 6, The track data is bed format, and color config data is either colormap name or json file path.
    """

    # first, we prepare output path
    try:
        # if it is a standard AF3 mmcif file
        job_name = re.match(r'fold_(.*)_model_\d+\.cif', os.path.basename(mmcif_file)).group(1) 
    except AttributeError:
        # fallback to general naming in case it is not standard AF3 mmcif file
        job_name = os.path.basename(mmcif_file).split(".")[0]

    # for computing contact map, we need to load the structure and extract representative atoms
    def _load_representative_atoms(mmcif_path, target_chains=None):
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('struct', mmcif_path)
        # AF3 usually has model 0
        model = structure[0]

        coords = []
        chain_labels = [] 
        res_types = [] # Track type for info   
        res_ids = [] # Store residue IDs for chain-based indexing if needed           
        
        # Normalize target_chains
        if isinstance(target_chains, str):
            target_chains = {target_chains}
        elif isinstance(target_chains, (list, tuple)):
            target_chains = set(target_chains)
        
        found_chains = set()

        for chain in model:
            if target_chains is not None and chain.id not in target_chains:
                continue
            
            for res in chain:
                # Skip water 
                if res.id[0] == "W": continue   
                # Note that HETATM residues are included here, cause Zn2+ is HETATM, so we don't skip them
                # but if we want to skip other HETATMs, we can add more conditions here.
                # like: if res.id[0].startswith("H_"): continue
 
                rep_atom = None      
                rtype = "Unknown"

                # 1. Protein -> CA
                if Polypeptide.is_aa(res, standard=True):
                    if 'CA' in res:
                        rep_atom = res['CA']
                        rtype = "Protein"
                
                # 2. Nucleic Acid -> C1' (Distinguish DNA/RNA)
                elif "C1'" in res:
                    rep_atom = res["C1'"]
                    rname = res.resname.strip()
                    if rname in ['DA', 'DT', 'DG', 'DC']:
                        rtype = "DNA"
                    elif rname in ['A', 'U', 'G', 'C']:
                        rtype = "RNA"
                    else:
                        rtype = "Nucleic" # Fallback for modified bases
                
                # 3. Others (Ligands/Ions) -> First Atom
                else:
                    atoms = list(res.get_atoms())
                    if len(atoms) > 0:
                        rep_atom = atoms[0]
                        rtype = "Ligand"

                if rep_atom:
                    coords.append(rep_atom.get_coord())
                    chain_labels.append(chain.id)
                    res_types.append(rtype)
                    res_ids.append(res.id[1]) # residue sequence number, 1-based
                    found_chains.add(chain.id)

        if not coords:
            msg = f"No valid atoms found in {mmcif_path}."
            if target_chains: msg += f" (Searched for chains: {target_chains})"
            raise ValueError(msg)

        return np.array(coords, dtype=np.float32), chain_labels, res_ids, sorted(list(found_chains))


    # We need 1D track data support here 
    def _get_aligned_track_data(track_cfg, chain_labels_list):
        """
        Description
        ------------
        Align 1D track data to the residues in the structure.
        Aligns track data (sparse dict) to the flat list of residues in the structure.
        
        Args
        ----
            track_cfg (dict): Configuration for the track data.
            chain_labels_list (list): List of chain IDs corresponding to each residue (rep atom).
        
        Notes
        -----
        - 1, The track_cfg is 0-based index !
        """
        
        # check first
        if not track_cfg or "track_data" not in track_cfg:
            return None
        
        full_data = []
        # our data
        raw_data = track_cfg["track_data"]

        # we traverse all residues in the structure
        current_chain = None
        res_counter = 0

        # Align data, we calculate residue index per chain
        for cid in chain_labels_list:
            if cid != current_chain:
                current_chain = cid
                res_counter = 0 # Note that track_cfg is 0-based index!
            else:
                res_counter += 1

            # default value is nan
            val = np.nan

            # check the raw_data
            if cid in raw_data:
                chain_data = raw_data[cid]
                # get bed value, default is nan
                # ⚠️ res_counter is 0-based index, chain_data is from track bed file, also 0-based index
                val = chain_data.get(res_counter, np.nan)

            # now we get the value for this residue
            full_data.append(val)
        # full_data now contains aligned values for all residues/rep atoms in the structure, whether the chain is present in track_data or not
        # that means, full_data is bed value aligned to chains in config_track, and nan for missing data
        # value column in track data can be numerical or categorical (str), so we just keep it 
        return np.asarray(full_data, dtype=object)



    # 1. Load Data
    # ⚠️ Note that res_ids from mmCIF are 1-based residue numbers ！
    coords, chain_labels, res_ids, loaded_chains = _load_representative_atoms(mmcif_file, chains)
    N = coords.shape[0] # number of tokens/representative atoms/residues
    print(f"Loaded {N} tokens from chains: {loaded_chains}")

    # 2. Compute Distance Matrix
    # shape (N, 1, 3) - (1, N, 3) -> (N, N, 3)
    diff = coords[:, None, :] - coords[None, :, :]
    dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))

    # 3. Prepare Tracks
    # first, we parse the color config if it's a json file path
    if color_config and color_config.endswith(".json") and os.path.isfile(color_config):
        with open(color_config, "r") as cf:
            color_config = json.load(cf)
    # else, it can be a colormap name already
    track_data = parse_bed_to_track_data(
        bed_file = track_bed_file,
        color = color_config)
    track_list = []
    if isinstance(track_data, dict): 
        track_list = [track_data]
    # see in utils/track_utils.py, parse_bed_to_track_data() will return a list normally
    elif isinstance(track_data, list): 
        track_list = track_data

    n_tracks = len(track_list)

    # 4. Setup Figure with Gridspec
    # Layout: [Top Tracks...] + [Chain Bar] + [Heatmap]
    # We need n_tracks + 2 rows, and n_tracks + 3 columns (Left Tracks + Chain + Heatmap + Cbar)

    # Ratios
    # ⚠️ categorical tracks are usually thinner, numerical tracks can be thicker, we need to distinguish them
    # track_ratio = 0.8  # Track height relative to others
    ratio_cat = 0.4 # Thinner for categorical (bar/strip)
    ratio_num = 1.0 # Taller for numerical (line plot)
    chain_ratio = 0.5  # Chain bar is thin
    main_ratio = 10.0  # Heatmap is large
    cbar_ratio = 0.4   # Colorbar width

        # Determine ratios for each track based on type
    track_ratios = []
    for track in track_list:
        t_type = track.get("track_type", "categorical")
        # Match the logic in plotting loop: "categorical" gets bar style, others get line style
        if t_type == "categorical":
            track_ratios.append(ratio_cat)
        else:
            track_ratios.append(ratio_num)

    # The GridSpec rows/cols for tracks are filled from outside (index 0) to inside (index n_tracks-1)
    # But our plotting loop places track i (from track_list) at index `n_tracks - 1 - i`.
    # So track 0 (inner most) corresponds to the last row of the track section.
    # Track N-1 (outer most) corresponds to the first row (index 0).
    # Therefore, we need to reverse the ratios list for GridSpec definition.
    track_ratios_gs = track_ratios[::-1]

    # Rows: Top Tracks (N) -> Chain Bar (1) -> Heatmap (1)
    height_ratios = track_ratios_gs + [chain_ratio, main_ratio]
    
    # Cols: Left Tracks (N) -> Chain Bar (1) -> Heatmap (1) -> Colorbar (1)
    width_ratios = track_ratios_gs + [chain_ratio, main_ratio, cbar_ratio]
    
    # Dynamic size according to number of tracks
    fig = plt.figure(figsize=(15 + n_tracks, 15 + n_tracks)) 
    gs = gridspec.GridSpec(
        nrows=n_tracks + 2, # Rows: Top Tracks (N) -> Chain Bar (1) -> Heatmap (1)
        ncols=n_tracks + 3, # Cols: Left Tracks (N) -> Chain Bar (1) -> Heatmap (1) -> Colorbar (1)
        figure=fig,
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        wspace=0.02, hspace=0.02 # Tight gap
    )

    # Indices for the Main Heatmap
    # cause we have n_tracks top and left, then 1 chain bar, so main heatmap is at (n_tracks+1, n_tracks+1)
    main_row_idx = n_tracks + 1
    main_col_idx = n_tracks + 1
    

    # --- Calculate Ticks (⚠️) ---
    # Logic similar to confidence_metrics_plot.py
    # we parametrize tick_step for flexibility
    # tick_step = 100 # Show tick every 100 residues
    xticks_loc = []
    xticks_labels = []
    
    # Note again, res_ids from _load_representative_atoms() are 1-based residue numbers !
    for i, (r_id, c_id) in enumerate(zip(res_ids, chain_labels)):
        # Check if it is the start of a new chain
        is_chain_start = (i == 0) or (chain_labels[i] != chain_labels[i-1])
        
        # now we calculate 0-based index
        r_id_0b = r_id - 1

        # Add tick if it's chain start or a multiple of tick_step
        # for example, residue 0-99, and res 100 is another 100-length segment, so we need to put tick at res 0 and res 100
        if is_chain_start or (r_id_0b % tick_step == 0):
            xticks_loc.append(i)
            xticks_labels.append(str(r_id_0b))



    # --- Plot Main Heatmap ---
    ax_main = fig.add_subplot(gs[main_row_idx, main_col_idx])
    im = ax_main.imshow(
        dist_matrix,
        cmap=cmap,
        origin="upper",
        vmin=0,
        interpolation="nearest",
        aspect='auto' # Important for GridSpec
    )
    ax_main.set_yticks([]) # Hide Y ticks
    ax_main.set_xlabel("Residue Index", fontsize=12)
    ax_main.set_xticks(xticks_loc)
    ax_main.set_xticklabels(xticks_labels, fontsize=8) # Smaller font to avoid overlap
    
    # --- Plot Colorbar ---
    ax_cbar = fig.add_subplot(gs[main_row_idx, main_col_idx + 1])
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Distance (Å)", fontsize=12)

    # --- Prepare Chain Blocks ---
    chain_blocks = []
    start = 0
    for i in range(1, N + 1):
        if i == N or chain_labels[i] != chain_labels[start]:
            chain_blocks.append((chain_labels[start], start, i - 1))
            start = i
    
    unique_chains = sorted(list(set(chain_labels)))
    chain_to_int = {cid: i for i, cid in enumerate(unique_chains)}
    chain_row = np.array([chain_to_int[c] for c in chain_labels]).reshape(1, -1)
    
    if len(unique_chains) <= 20: 
        cmap_chains = plt.get_cmap("tab20", len(unique_chains))
    else: 
        # set 0.1 to 0.9 to avoid too light/dark colors
        cmap_chains = ListedColormap(plt.get_cmap("gist_rainbow")(np.linspace(0.1, 0.9, len(unique_chains))))

    # --- Plot Chain Bars ---
    # Top Chain Bar
    ax_chain_top = fig.add_subplot(gs[main_row_idx - 1, main_col_idx], sharex=ax_main)
    ax_chain_top.imshow(chain_row, cmap=cmap_chains, aspect="auto")
    ax_chain_top.set_yticks([])
    # ⚠️⚠️ax_chain_top.tick_params(labelbottom=False) # Hide x labels
    # FIX: Use tick_params to hide ticks visually instead of removing them, to preserve ax_main ticks
    ax_chain_top.tick_params(axis='x', which='both', bottom=False, labelbottom=False) 
    
    # Left Chain Bar
    ax_chain_left = fig.add_subplot(gs[main_row_idx, main_col_idx - 1], sharey=ax_main)
    ax_chain_left.imshow(chain_row.T, cmap=cmap_chains, aspect="auto")
    ax_chain_left.set_xticks([])
    # ⚠️⚠️ax_chain_left.tick_params(labelleft=False)
    # FIX: Use tick_params to hide ticks visually
    ax_chain_left.tick_params(axis='y', which='both', left=False, labelleft=False)

    # --- Plot 1D Tracks ---
    # Loop from inner (closest to chain bar) to outer
    for i, track_cfg in enumerate(track_list):
        # Calculate Grid Indices
        # Top tracks stack upwards: index = (main_row - 1) - 1 - i
        row_idx = (main_row_idx - 1) - 1 - i
        col_idx = main_col_idx
        
        # Left tracks stack leftwards: index = (main_col - 1) - 1 - i
        row_idx_l = main_row_idx
        col_idx_l = (main_col_idx - 1) - 1 - i
        
        # Create Axes
        # top track axes, share x with main heatmap
        ax_t_top = fig.add_subplot(gs[row_idx, col_idx], sharex=ax_main)
        # left track axes, share y with main heatmap
        ax_t_left = fig.add_subplot(gs[row_idx_l, col_idx_l], sharey=ax_main)
        
        # Get Aligned Data
        # align track data to the residue list
        track_vals = _get_aligned_track_data(track_cfg, chain_labels)
        if track_vals is None: 
            continue
            
        # Get Track Config
        # track color/type/name, defult value is useless here, so we just set something
        t_color = track_cfg.get("color", "tab10c")
        t_type = track_cfg.get("track_type", "categorical") # categorical or numerical
        t_name = track_cfg.get("track_name", "")
        
        x_indices = np.arange(N)
        
        # --- Plotting Logic ---
        if t_type == "categorical":
            # Use imshow for categorical data (cleaner than bar)
            # 1. Map categories to integers
            # Handle NaNs by assigning them to -1 or a specific index
            str_vals = [str(v) for v in track_vals]
            # Get color map
            color_map = t_color if isinstance(t_color, dict) else {}
            
            # Create a list of colors for the colormap
            # We need to ensure the integer mapping matches the colormap order
            unique_cats = sorted(list(set([v for v in str_vals if v != "nan"])))
            cat_to_int = {cat: i for i, cat in enumerate(unique_cats)}
            
            # Create integer array for imshow
            int_row = np.full((1, N), -1) # -1 for background/nan
            for idx, v in enumerate(str_vals):
                if v in cat_to_int:
                    int_row[0, idx] = cat_to_int[v]
            
            # Create Colormap
            # Colors must match the integer indices
            colors_list = [color_map.get(cat, "#808080") for cat in unique_cats]
            if not colors_list: colors_list = ["#FFFFFF"] # Fallback
            custom_cmap = ListedColormap(colors_list)
            
            # Mask the background (-1) to make it transparent
            masked_row = np.ma.masked_where(int_row == -1, int_row)

            # Plot Top
            # We need to mask -1 values to be transparent or white
            if len(unique_cats) > 0:
                # Determine vmin/vmax to ensure correct color mapping
                vmin, vmax = 0, len(unique_cats) - 1
                if vmin == vmax: # Single category case
                    vmin -= 0.5
                    vmax += 0.5

                ax_t_top.imshow(masked_row, cmap=custom_cmap, aspect="auto", interpolation="nearest", vmin=vmin, vmax=vmax)
                
                # Plot Left (Transpose)
                ax_t_left.imshow(masked_row.T, cmap=custom_cmap, aspect="auto", interpolation="nearest", vmin=vmin, vmax=vmax)
            
        else:
            # Numerical: Plot Line without Fill (there is no need to fill, cause we do not know the baseline)
            try:
                # Convert to float, invalid parsing will be nan
                vals = track_vals.astype(float)
            except:
                vals = np.full(N, np.nan)
            
            # ---- Top Plot (line) -------
            ax_t_top.plot(x_indices, vals, color=t_color, linewidth=1)
            
            # Get valid mask
            valid_mask = ~np.isnan(vals)

            # ⚠️ Deprecated !
            # Fill area under curve to min value
            '''
            if np.any(valid_mask):
                min_val = np.nanmin(vals)
                ax_t_top.fill_between(x_indices, vals, min_val, color=t_color, alpha=0.4)
            '''
            
            # Left Plot (Rotated)
            ax_t_left.plot(vals, x_indices, color=t_color, linewidth=1)
            
            # ⚠️ Deprecated also !
            # Fill area to min value
            '''
            if np.any(valid_mask):
                ax_t_left.fill_betweenx(x_indices, vals, min_val, color=t_color, alpha=0.4)
            '''
            
            # Limits
            if np.any(valid_mask):
                vmin, vmax = np.nanmin(vals), np.nanmax(vals)
                # ⚠️ Add some margin, 10% of range
                # if the output plot is not good, we can adjust this margin ratio, or just skip margin for original data range
                margin = (vmax - vmin) * 0.1
                # So that line is not at the edge, the actual data range is larger than vmin/vmax
                # the curve is vmin/vmax, but the track axes limit is extended, so line is not at the edge
                ax_t_top.set_ylim(vmin - margin, vmax + margin)
                ax_t_left.set_xlim(vmin - margin, vmax + margin) # Note xlim for left plot

        # --- Styling ---
        # Hide ticks
        # ⚠️⚠️FIX: Do not use set_xticks([]) on shared axes (ax_t_top shares X, ax_t_left shares Y)
        # This prevents wiping out the ticks on the main heatmap
        
        # ax_t_top.set_xticks([]); ax_t_top.set_yticks([])
        # ax_t_left.set_xticks([]); ax_t_left.set_yticks([])
        
        # Add Label
        # ax_t_top.set_ylabel(t_name, fontsize=9, rotation=0, ha="right", va="center")
        # Top Track: Shares X. Hide X ticks/labels visually. Y is independent.
        ax_t_top.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        # ax_t_top.set_yticks([])
        
        # Left Track: Shares Y. Hide Y ticks/labels visually. X is independent.
        # ax_t_left.set_xticks([])
        ax_t_left.tick_params(axis='y', which='both', left=False, labelleft=False)
        
        # ⚠️ Conditional tick display based on track type
        if t_type == "categorical":
            # For categorical, hide the value axis ticks (Y for top, X for left)
            ax_t_top.set_yticks([])
            ax_t_left.set_xticks([])
        else:
            # For numerical, show ticks to indicate scale
            # Top Track: Value is on Y axis
            ax_t_top.tick_params(axis='y', labelsize=6)
            ax_t_top.locator_params(axis='y', nbins=3) # Limit to ~3 ticks (min, mid, max) to avoid clutter
            
            # Left Track: Value is on X axis
            # To match the "Left" Y-axis of the Top Track (after 90 deg rotation),
            # we should place the X-axis ticks on the TOP.
            ax_t_left.xaxis.tick_top() 
            ax_t_left.tick_params(axis='x', labelsize=6, labelrotation=-90) # Rotate labels for narrow tracks
            ax_t_left.locator_params(axis='x', nbins=3)


        # Add Label
        ax_t_top.set_ylabel(t_name, fontsize=9, rotation=0, ha="right", va="center")



        # Add Chain Boundaries
        for cid, s, e in chain_blocks:
            if s != 0:
                sep = s - 0.5
                ax_t_top.axvline(sep, color='black', linestyle='--', linewidth=0.5, alpha=0.3)
                ax_t_left.axhline(sep, color='black', linestyle='--', linewidth=0.5, alpha=0.3)

    # --- Final Touches ---
    # Add Chain Boundaries to Main Heatmap
    for cid, s, e in chain_blocks:
        if s != 0:
            sep = s - 0.5
            ax_main.axvline(sep, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
            ax_main.axhline(sep, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
        
        # Add Labels to Chain Bars
        center = (s + e) / 2.0
        ax_chain_top.text(center, 0, str(cid), ha="center", va="center", fontsize=10, weight="bold")
        ax_chain_left.text(0, center, str(cid), ha="center", va="center", fontsize=10, weight="bold", rotation=90)

    # Title (Set on the top-most track axes or figure)
    fig.suptitle(f"Contact Map: {job_name}", fontsize=16, y=0.9) # default y=0.98

    plt.savefig(f"{out_path}/{job_name}_contact_map.pdf", bbox_inches='tight')
    plt.savefig(f"{out_path}/{job_name}_contact_map.png", bbox_inches='tight', dpi=300)
    plt.close(fig)
