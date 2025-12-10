import os
import re
import numpy as np
from typing import Optional, List, Union, Dict, Any
from Bio.PDB import MMCIFParser, Polypeptide
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

def contact_map_vis_without_track(
    mmcif_file: str,
    chains: Optional[Union[str, List[str]]] = None,
    out_path: Optional[str] = None,
    cmap: str = "RdBu", # Reversed RdBu_r so Red is close (contact), Blue is far
    tick_step: int = 100
):
    """
    Description
    -----------
    Generate a contact map (distance matrix) for structures from an mmCIF file without 1D TRACK overlays.
    Uses representative atoms: CA for Protein, C1' for Nucleic Acids, First Atom for others.
    
    Args
    -----
        mmcif_file (str): Path to the mmcif file.
        chains (str | List[str], optional): Chain ID(s) to include. 
                                            e.g., "A" or ["A", "B"]. 
                                            If None, all chains are included.
        out_path (str, optional): Path to save the output plot (e.g., "/data2").
        cmap (str): Colormap to use. Default is "RdBu" (Red=Close, Blue=Far).
        tick_step (int): Step size for ticks on the axes. Default is 100.
    
    Notes
    ------
    - 1, THIS version of contact map visualization does NOT support 1D TRACK overlays. IF you need that, please see contact_map_visualization_with_track.py
    """

    # first, we prepare output path
    try:
        # if it is a standard AF3 mmcif file
        job_name = re.match(r'fold_(.*)_model_\d+\.cif', os.path.basename(mmcif_file)).group(1) 
    except AttributeError:
        # fallback to general naming in case it is not standard AF3 mmcif file
        job_name = os.path.basename(mmcif_file).split(".")[0]

    def _load_representative_atoms(mmcif_path, target_chains=None):
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('struct', mmcif_path)
        # AF3 usually has model 0
        model = structure[0]

        coords = []
        chain_labels = [] 
        res_types = [] # Track type for info    
        res_ids = [] # Store residue IDs for tick labeling          
        
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
                    res_ids.append(res.id[1]) # Capture 1-based residue ID
                    found_chains.add(chain.id)

        if not coords:
            msg = f"No valid atoms found in {mmcif_path}."
            if target_chains: msg += f" (Searched for chains: {target_chains})"
            raise ValueError(msg)

        return np.array(coords, dtype=np.float32), chain_labels, res_ids, sorted(list(found_chains))

    # 1. Load Data
    # ⚠️ Note that res_ids are 1-based residue numbers from the mmCIF file
    coords, chain_labels, res_ids, loaded_chains = _load_representative_atoms(mmcif_file, chains)
    N = coords.shape[0] # number of tokens/representative atoms/residues
    print(f"Loaded {N} tokens from chains: {loaded_chains}")

    # 2. Compute Distance Matrix
    # shape (N, 1, 3) - (1, N, 3) -> (N, N, 3)
    diff = coords[:, None, :] - coords[None, :, :]
    dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))

    # 3. Plotting
    fig, ax = plt.subplots(figsize=(15, 12))

    # Calculate Ticks Logic (Same as in with_track module)
    # we parametrize tick_step for flexibility
    # tick_step = 100 # Show tick every 100 residues
    xticks_loc = []
    xticks_labels = []

    for i, (r_id, c_id) in enumerate(zip(res_ids, chain_labels)):
        # Check if it is the start of a new chain
        is_chain_start = (i == 0) or (chain_labels[i] != chain_labels[i-1])
        
        # Convert to 0-based index (PDB is 1-based)
        r_id_0b = r_id - 1

        # Add tick if it's chain start or a multiple of tick_step
        if is_chain_start or (r_id_0b % tick_step == 0):
            xticks_loc.append(i)
            xticks_labels.append(str(r_id_0b))
    
    im = ax.imshow(
        dist_matrix,
        cmap=cmap,
        origin="upper",
        vmin=0,
        interpolation="nearest" 
    )
    
    # Create divider for existing axes instance
    divider = make_axes_locatable(ax)
    
    # Append axes to the top and left for chain bars
    ax_top = divider.append_axes("top", size="5%", pad=0.03)
    ax_left = divider.append_axes("left", size="5%", pad=0.03) 
    
    # Append axes for colorbar to the right
    cax = divider.append_axes("right", size="5%", pad=0.1)
    
    # Apply Ticks to Main Heatmap
    ax.set_xticks(xticks_loc)
    ax.set_xticklabels(xticks_labels, fontsize=8, rotation=0)
    
    # Add Colorbar
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Distance (Å)", fontsize=12)

    # Prepare chain blocks and colors
    chain_blocks = []
    start = 0
    for i in range(1, N + 1):
        if i == N or chain_labels[i] != chain_labels[start]:
            chain_blocks.append((chain_labels[start], start, i - 1))
            start = i
            
    unique_chains = sorted(list(set(chain_labels)))
    chain_to_int = {cid: i for i, cid in enumerate(unique_chains)}
    chain_row = np.array([chain_to_int[c] for c in chain_labels]).reshape(1, -1)
    
    # Use tab20 for distinct chain colors
    if len(unique_chains) <= 20:
        cmap_chains = plt.get_cmap("tab20", len(unique_chains))
    else:
        # set 0.1 to 0.9 to avoid too light/dark colors
        colors = plt.get_cmap("gist_rainbow")(np.linspace(0.1, 0.9, len(unique_chains)))
        cmap_chains = ListedColormap(colors)
    
    # Plot chain bars
    ax_top.imshow(chain_row, cmap=cmap_chains, aspect="auto", alpha=0.9)
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    
    ax_left.imshow(chain_row.T, cmap=cmap_chains, aspect="auto", alpha=0.9)
    ax_left.set_xticks([])
    ax_left.set_yticks([])
    
    # Remove y-axis ticks on the main plot to prevent overlap with the left chain bar
    ax.set_yticks([])

    # Add Chain Boundaries and Labels
    for cid, s, e in chain_blocks:  
        # Draw boundary lines
        if s != 0:
            sep = s - 0.5
            ax.axvline(sep, color='black', linestyle='--', linewidth=1, alpha=0.7)
            ax.axhline(sep, color='black', linestyle='--', linewidth=1, alpha=0.7)
            # White lines on the chain bars
            ax_top.axvline(sep, color="w", linewidth=1)
            ax_left.axhline(sep, color="w", linewidth=1)
            
        # Add labels 
        center = (s + e) / 2.0 
        ax_top.text(center, 0, str(cid), ha="center", va="center", 
                    fontsize=14, weight="bold", color="#222222") 
        ax_left.text(0, center, str(cid), ha="center", va="center", 
                     fontsize=14, weight="bold", color="#222222", rotation=90)

    # Set titles and labels
    ax_top.set_title(f"Contact Map: {job_name}", fontsize=14, pad=10) 
    ax.set_xlabel("Residue Index (Protein: CA, DNA/RNA: C1', Ligand: Atom1)", fontsize=12) 
    ax.set_ylabel("Residue Index", fontsize=12) 

    # 4. Save Output 
    plt.savefig(f"{out_path}/{job_name}_contact_map.pdf", bbox_inches='tight') 
    # also save as png, dpi 300 
    plt.savefig(f"{out_path}/{job_name}_contact_map.png", bbox_inches='tight', dpi=300) 
    plt.close(fig) 