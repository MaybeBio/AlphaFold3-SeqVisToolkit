# Note: 

# Input: 2 mmcif files generated from AlphaFold3 predictions
# - 1, the object compared in 2 mmcif files should be the same protein sequence. For example, both files are predictions of the same protein sequence but under different conditions or with different ligands.

# Output: A 2x2 contact difference map
# - 1, leftupper: Contact map of structure 1, rightupper: Contact map of structure 2, leftlower: Contact difference map (structure 1 - structure 2), rightlower: Contact difference map (structure 2 - structure 1).

# reference: https://biopython.org/docs/1.75/api/Bio.PDB.html


from typing import Optional, Tuple, List, Union
from Bio.PDB import MMCIFParser, Polypeptide
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap 
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import os
import re

def contact_map_diff_multimer(
        mmcif_file_a: str,
        mmcif_file_b: str,
        chain_a: str, # Can be "A" or "A,B"
        chain_b: str, # Can be "C" or "C,D"
        region_1: Union[Tuple[int, int], List[int], str] = None,
        region_2: Optional[Union[Tuple[int, int], List[int], str]] = None,  
        region_pairs: Optional[Union[List[Tuple[Tuple[int, int], Tuple[int, int]]], List[str]]] = None,
        vmax : Optional[float] = None,
        vmax_percentile: float = 95.0,
        vdiff: Optional[float] = None,
        vdiff_percentile: float = 95.0,
        include_nonstandard_residue: bool = False,
        return_maxtrix = False,
        out_path: Optional[str] = None,
        cmap_dist = "RdBu", 
        cmap_diff = "seismic",
        tick_step: int = 100 
):
    """
    Description
    --
        Generate a 2x2 contact difference map from two mmcif files.

    Args
        mmcif_file_a (str): Path to the first mmcif file.
        mmcif_file_b (str): Path to the second mmcif file.
        chain_a (str): Chain ID to include from the first mmcif file. We currently do not support None/all-chains option here. Please provide a specific chain ID.
        chain_b (str): Chain ID to include from the second mmcif file. 
        region_1 (Tuple[int, int] | List[int, int]): 0-based Region of interest in the first structure (start, end) or list of residue indices.
            for example, (10, 50) or [10, 11, 12, ..., 50], or [10, 50]
        region_2 (Tuple[int, int] | List[int, int], optional): 0-based Region of interest in the second structure (start, end) or list of residue indices. Defaults to None, which means using region_1.
        regions_pair (List[Tuple[Tuple, Tuple]] | List[str], optional): List of region pairs for batch processing. Each item is a tuple of two regions (region_1, region_2) or a string "start1:end1,start2:end2". Defaults to None.
        vmax (float): Maximum distance value for color bar scaling. If None, it will be determined automatically based on the data. Defaults to None.
        vmax_percentile (float): Percentile value to determine vmax if vmax is None. Defaults to 95.0.
        vdiff (float): Maximum absolute difference value for color bar scaling. If None, it will be determined automatically based on the data. Defaults to None.
        vdiff_percentile (float): Percentile value to determine vdiff if vdiff is

        include_nonstandard_residue (bool): Whether to include non-standard residues like MSE. Defaults to False.
        return_maxtrix (bool): Whether to return the contact matrices of the selected regions and related information. Defaults to False.
        cmap_dist (str): Colormap for distance maps. Defaults to "RdBu".(red for close, blue for far)
        cmap_diff (str): Colormap for difference maps. Defaults to "seismic".(red for positive, blue for negative)
        out_path (str, optional): Directory to save the output figure. 
        tick_step (int): Step size for ticks on the axes. Defaults to 100.

    Returns
        None: Displays a 2x2 contact difference map.

    Notes:
    - 1, region_1 and region_2 are both 0-based residue indices, like CTCF (0,726)
    - 2, current logic only supports like: if only region_1 provided, then we draw region_1 x region_1 submaxtrix;
    if both region_1 and region_2 are provided, then we draw region_1 x region_2 submatrix.
    - 3, if regions_pair is provided, region_1 and region_2 will be ignored, and we will batch process all region pairs in regions_pair.
    Format of regions_pair can be:
        - a, List of tuple of tuples: [ ((start1, end1), (start2, end2)), ... ]
        - b, List of strings: [ "start1:end1,start2:end2", ... ], ["start1-end1,start2-end2"]
    - 4, We currently do not support None/all-chains option here. Please provide a specific chain ID.
    
    """
    
    # first, we prepare output path
    try:
        # if it is a standard AF3 mmcif file
        job_a = re.match(r'fold_(.*)_model_\d+\.cif', os.path.basename(mmcif_file_a)).group(1) 
        job_b = re.match(r'fold_(.*)_model_\d+\.cif', os.path.basename(mmcif_file_b)).group(1)
    except AttributeError:
        # fallback to general naming in case it is not standard AF3 mmcif file
        job_a = os.path.basename(mmcif_file_a).split(".")[0]
        job_b = os.path.basename(mmcif_file_b).split(".")[0]

    chain_a_str = chain_a.replace(",", "_")
    chain_b_str = chain_b.replace(",", "_")
    job_name = f"{job_a}_{chain_a_str}_vs_{job_b}_{chain_b_str}"

    # then, we parse chain inputs into lists
    chains_a_list = [c.strip() for c in chain_a.split(",")]
    chains_b_list = [c.strip() for c in chain_b.split(",")]
    if len(chains_a_list) != len(chains_b_list):
        raise ValueError(f"Number of chains must match! A: {len(chains_a_list)} vs B: {len(chains_b_list)}")

    # Finally, we need get the correct region index first
    # _parse_region for single region
    def _parse_region(r, boundaries=None):
        """
        supporting format like (start, end) or [start, end] or "start:end" or "start-end"
        Also supports "Chain:start-end" or "Chain:start:end" if boundaries is provided.
        """
        offset = 0
        # Handle Chain:Region format
        if boundaries and isinstance(r, str) and ":" in r:
            parts = r.split(":", 1)
            # Check if the first part is a valid chain ID in the boundaries
            # boundaries structure: [(chain_id, start_idx, end_idx), ...]
            chain_id = parts[0].strip()
            found_chain = next((b for b in boundaries if b[0] == chain_id), None)
            
            if found_chain:
                offset = found_chain[1] # The start index of this chain
                r = parts[1] # The rest is the region string

        if isinstance(r,(tuple,list)) and len(r) == 2:
            s, e = int(r[0]), int(r[1])
        elif isinstance(r, str) and (":" in r or "-" in r):
            sep = ":" if ":" in r else "-"
            s, e = r.split(sep)
            s, e = int(s.strip()), int(e.strip())
        else:
            raise ValueError(f"Illegal region format: {r}")
        if e < s:
            raise ValueError(f"Illegal region format: {r}, end < start.")
        return (s + offset, e + offset)

    # _parse_region_pairs for batch region pairs
    def _parse_region_pairs(pair_str, boundaries=None) -> Tuple[Tuple[int, int], Tuple[int, int]]:
        """
        supporting format like "start1:end1,start2:end2" or "start1-end1,start2-end2"
        """
        pair_str = pair_str.strip()
        sep = "," if "," in pair_str else None
        if not sep:
            raise ValueError(f"Illegal region pair format (missing ','): {pair_str}")
        left, right = pair_str.split(sep, 1)
        r1 = _parse_region(left.strip(), boundaries)
        r2 = _parse_region(right.strip(), boundaries)
        return (r1, r2)
    
    # then, we load the result of structure prediction
    # Cause we generally only compare proteins, so here, we only focus on CA atoms ! ! !
    # That means we choose CA atoms as the representative atom for each residue
    # ⚠️ Of course, we can extend this further, choose C1' for nucleotides, Atom 1 for ligands, etc. And we can compare different types of molecules before and after some treatments.
    # And note that each residue should have only one CA atom, so we can directly use residue index to access CA atom
    # For loading representative atoms from mmcif file for different molecule types, please refer to the module: contact_map_vis_no/Track.py
    '''
    def _load_ca(mmcif_file, target_chains_list, include_nonstandard_residue=False, model_index=0):
        """
        Description
        -----------
            Load CA atom coordinates for a list of chains in the specific order provided.

        Args
        ----
        mmcif_file: str, path to the mmcif file
        target_chains_list: list of str, chain IDs to include in the specific order provided. If None, all chains are included.
        include_nonstandard_residue: bool, whether to include non-standard residues like MSE
        model_index: int, index of the model to load (default: 0)

        Notes
        -----
        - 1, For loading representative atoms from mmcif file for different molecule types, please refer to the module: contact_map_vis_no/Track.py
        """
        
        # for af3, we generally need only the first model, that is model_0
        parser = MMCIFParser()
        # According to SMCRA hierarchy, we need to go through Structure -> Model -> Chain -> Residue -> Atom
        structure = parser.get_structure('struct', mmcif_file)
        models = list(structure)
        if not models:
            raise ValueError(f"No models found in the mmcif file: {mmcif_file}")
        # Generally, there is only one model in the mmcif file generated by AlphaFold3
        model = models[model_index]

        all_ca_coords = []
        all_ca_info = []
        chain_boundaries = [] # store (chain_id, start_idx, end_idx)
        current_idx = 0

        # Create a lookup dict for chains to avoid repeated looping over the model
        model_chains = {chain.id: chain for chain in model}

        # Iterate strictly according to the user-provided order
        for target_chain_id in target_chains_list:
            if target_chain_id not in model_chains:
                raise ValueError(f"Chain {target_chain_id} not found in {mmcif_file}")
            
            chain = model_chains[target_chain_id]
            chain_start_idx = current_idx

            for res in chain:
                hetfield, resseq, icode = res.id # or .get_id()
                if hetfield != " ":
                    # that means it's a hetero atom, we skip it
                    continue
                if not Polypeptide.is_aa(res, standard=not include_nonstandard_residue):
                    # that means it's a non-standard residues/amino acids, we also skip it
                    continue
                if 'CA' not in res:
                    # that means there is no CA atom in this residue, we skip it
                    continue
                ca = res['CA']
                # get the coordinate of CA atom
                all_ca_coords.append(ca.get_coord())
                all_ca_info.append(
                    {
                        "chain": chain.id, # like "A"
                        "resname": res.resname, # like "ALA"
                        "resseq": resseq, # like position 1(1-based)
                        "icode": (icode if isinstance(icode, str) and icode.strip() else "") # like 'A' for Thr 80 A，Ser 80 B
                    }
                )
                current_idx += 1

            # Record boundary if atoms were added for this chain
            if current_idx > chain_start_idx:
                chain_boundaries.append((target_chain_id, chain_start_idx, current_idx - 1))
            else:
                print(f"Warning: Chain {target_chain_id} found but no valid CA atoms loaded.")

    
        if not all_ca_coords:
            raise ValueError(f"No CA atoms found for chains {target_chains_list} in {mmcif_file}")
            
        return np.asarray(all_ca_coords, dtype=np.float32), all_ca_info, chain_boundaries
    '''

    # Generalized atom loading function
    def _load_representative_atoms(mmcif_file, target_chains_list, include_nonstandard_residue=False, model_index=0):
        """
        Description
        -----------
            Load representative atoms (CA for protein, C1' for nucleic, Atom1 for ligand) 
        for a list of chains in the specific order provided.

        Args
        ----
        mmcif_file: str, path to the mmcif file
        target_chains_list: list of str, chain IDs to include in the specific order provided. If None, all chains are included.
        include_nonstandard_residue: bool, whether to include non-standard residues like MSE
        model_index: int, index of the model to load (default: 0)

        Notes
        -----
        - 1, For loading representative atoms from mmcif file for different molecule types, please refer to the module: contact_map_vis_no/Track.py
        """
        
        # for af3, we generally need only the first model, that is model_0
        parser = MMCIFParser()
        # According to SMCRA hierarchy, we need to go through Structure -> Model -> Chain -> Residue -> Atom
        structure = parser.get_structure('struct', mmcif_file)
        models = list(structure)
        if not models:
            raise ValueError(f"No models found in the mmcif file: {mmcif_file}")
        # Generally, there is only one model in the mmcif file generated by AlphaFold3
        model = models[model_index]

        all_coords = []
        all_info = []
        chain_boundaries = [] # store (chain_id, start_idx, end_idx)

        current_idx = 0
        # Create a lookup dict for chains to avoid repeated looping over the model
        model_chains = {chain.id: chain for chain in model}

        # Iterate strictly according to the user-provided order
        for target_chain_id in target_chains_list:
            if target_chain_id not in model_chains:
                raise ValueError(f"Chain {target_chain_id} not found in {mmcif_file}")
            
            chain = model_chains[target_chain_id]
            chain_start_idx = current_idx

            for res in chain:
                # Skip water
                if res.id[0] == "W": 
                    continue

                rep_atom = None
                rtype = "Unknown"

                # 1. Protein -> CA
                if Polypeptide.is_aa(res, standard=not include_nonstandard_residue):
                    if 'CA' in res:
                        rep_atom = res['CA']
                        rtype = "Protein"
                
                # 2. Nucleic Acid -> C1'
                elif "C1'" in res:
                    rep_atom = res["C1'"]
                    rname = res.resname.strip()
                    if rname in ['DA', 'DT', 'DG', 'DC']: rtype = "DNA"
                    elif rname in ['A', 'U', 'G', 'C']: rtype = "RNA"
                    else: rtype = "Nucleic"
                
                # 3. Others (Ligands/Ions) -> First Atom
                else:
                    # Only include if we want to be permissive, or check specific types
                    # For comparison, usually we want everything in the chain
                    atoms = list(res.get_atoms())
                    if len(atoms) > 0:
                        rep_atom = atoms[0]
                        rtype = "Ligand"

                if rep_atom:
                    all_coords.append(rep_atom.get_coord())
                    all_info.append({
                        "chain": chain.id,
                        "resname": res.resname,
                        "resseq": res.id[1],
                        "icode": res.id[2],
                        "type": rtype
                    })
                    current_idx += 1

            # Record boundary if atoms were added for this chain
            if current_idx > chain_start_idx:
                chain_boundaries.append((target_chain_id, chain_start_idx, current_idx - 1))
            else:
                print(f"Warning: Chain {target_chain_id} found but no valid representative atoms loaded.")
        
        if not all_coords: 
            raise ValueError(f"No valid atoms found for chains {target_chains_list}.")
        
        return np.asarray(all_coords, dtype=np.float32), all_info, chain_boundaries



    def _pairwise_dist(ca_coords):
        """
        For a protein structure with N residues(equals N CA atoms), compute the pairwise distance matrix of shape (N, N)
        ca_coords: np.ndarray of shape (N, 3), coordinates of CA atoms
        """

        # genarally we may use double for loop to compute pairwise distance matrix
        # but here we use broadcasting to accelerate the computation
        diff = ca_coords[:, None, :] - ca_coords[None, :, :]
        return np.sqrt(np.sum(diff * diff, axis=-1))
    
    # load representative atom coordinates and info from mmcif files
    ca_coords_a, ca_info_a, boundaries_a = _load_representative_atoms(mmcif_file_a, chains_a_list, include_nonstandard_residue=include_nonstandard_residue)
    ca_coords_b, ca_info_b, boundaries_b = _load_representative_atoms(mmcif_file_b, chains_b_list, include_nonstandard_residue=include_nonstandard_residue)
    
    # Validate lengths per chain pair
    # Since we are comparing A+B vs C+D, we expect len(A)=len(C) and len(B)=len(D)
    if len(boundaries_a) != len(boundaries_b):
         raise ValueError("Number of non-empty chains loaded differs between files.")
    
    for (chain_a_id, start_a, end_a), (chain_b_id, start_b, end_b) in zip(boundaries_a, boundaries_b):
        len_a = end_a - start_a + 1
        len_b = end_b - start_b + 1
        if len_a != len_b:
            raise ValueError(f"Chain length mismatch between structures for chain pair {chain_a_id} and {chain_b_id}: {len_a} vs {len_b}. Please check whether the two mmcif files are generated from the same protein sequence.")
    
    
    
    # calculate the shape of each protein structure
    # cause we are evaluateing the same protein sequence, so Na should be equal to Nb
    Na, Nb = ca_coords_a.shape[0], ca_coords_b.shape[0]
    if Na != Nb:
        raise ValueError(f"The number of residues in the protein part of two structures are different: {Na} vs {Nb}. Please check whether the two mmcif files are generated from the same protein sequence,\n \
                        cause we are only focusing on the protein regions, so they should match.")
    N = Na  # since Na == Nb

    # parse the region of interest
    pairs: List[Tuple[Tuple[int,int], Tuple[int,int]]] = []
    is_default_region = False # Flag for default full length

    # if region_pairs is provided, we will prefer it for batch processing, otherwise we use region_1 and region_2
    if region_pairs:
        # format like [((1,2),(3,4)), ...] or ["1:2,3:4", ...]
        # if item is str, we need parse it into tuple of tuples (the first format) first
        if isinstance(region_pairs[0], str):
            pairs = [_parse_region_pairs(p, boundaries_a) for p in region_pairs]
        else:
            pairs = [(_parse_region(a, boundaries_a), _parse_region(b, boundaries_b)) for (a,b) in region_pairs]
            # pairs = region_pairs
    # else, we use region_1 and region_2
    # If no regions are specified, default to the whole structure and set the flag
    elif region_1 is None:
        pairs = [((0, N-1), (0, N-1))]
        is_default_region = True
    # if only region_1 is provided, we compare region_1 x region_1
    else:
        r1 = _parse_region(region_1, boundaries_a)
        r2 = _parse_region(region_2, boundaries_b) if region_2 is not None else r1
        pairs = [(r1, r2)]

    # check whether the region is valid
    # Note: we are using 0-based residue index here
    # Note again: we are comparing regions between two structures with the same protein sequence, so N1 = N2 
    for (ra, rb) in pairs:
        for r in (ra, rb):
            if r[0] < 0 or r[1] > N-1:
                raise ValueError(f"Region {r} is out of bounds for protein size {N} (0-based).")

    # calculate pairwise distance matrices for both structures
    dist_a = _pairwise_dist(ca_coords_a)
    dist_b = _pairwise_dist(ca_coords_b)
    diff_ab = dist_a - dist_b
    diff_ba = - diff_ab

    # calculate the color bar limits
    # Here, we want to compare the distance maps at the same scale, so we should unify vmax from the real data of both max values
    # Generally, we may use the max value of both dist maximums as the unified vmax
    # But considering that there may be some extreme values in the distance maps and they maybe noise, we use the 95th percentile value instead
    # For minimum value, we do not set 5th percentile, because the minimum distance should be 0 anyway, it is rational
    # So we set the vmax_percentile parmeter to contorl the color bar limits, default is 95.0, and 100 means using the real max value
    if vmax is None:
        vmax_use = float(max(
            np.nanpercentile(dist_a, vmax_percentile),
            np.nanpercentile(dist_b, vmax_percentile)
        ))
    else:
        vmax_use = float(vmax)

    # the same for vdiff (0 centered)
    if vdiff is None:
        abs_all = np.abs(np.concatenate([diff_ab.ravel(),diff_ba.ravel()]))
        vdiff_use = float(np.nanpercentile(abs_all, vdiff_percentile))
    else:
        vdiff_use = float(vdiff)

    # ----PLOTTING SECTION----
    # set 2x2 subplots, better layout with constrained_layout=True
    fig, axes = plt.subplots(2, 2, figsize=(20, 15))

    # Titles for each subplot
    title_a = f"A: {job_a}"
    if chain_a: 
        title_a += f" (Chain {chain_a})"

    title_b = f"B: {job_b}"
    if chain_b: 
        title_b += f" (Chain {chain_b})"
    
    titles = [
        title_a,
        title_b,
        "A - B (A->B: Red = closer, blue = farther)",
        "B - A (B->A: Red = closer, blue = farther)"
    ]

    # Draw Images!
    # plot distance map of structure 1, protein marked as A/1
    # A/1
    imA = axes[0,0].imshow(
        dist_a, 
        cmap=cmap_dist,
        vmin=0,
        vmax=vmax_use,
        interpolation="nearest",
        origin="upper"
    )
    # Update titles to show Chain ID
    # axes[0,0].set_title(f"A:{os.path.splitext(os.path.basename(mmcif_file_a))[0].split('_')[1]}")
    # title_a = f"A: {job_a}"
    # if chain_a: title_a += f" (Chain {chain_a})"
    # axes[0,0].set_title(title_a)

    # plot distance map of structure 2, protein marked as B/2
    # B/2
    imB = axes[0,1].imshow(
        dist_b,
        cmap=cmap_dist,
        vmin=0,
        vmax=vmax_use,
        interpolation="nearest",
        origin="upper"
    )
    # axes[0,1].set_title(f"B:{os.path.splitext(os.path.basename(mmcif_file_b))[0].split('_')[1]}")
    # title_b = f"B: {job_b}"
    # if chain_b: title_b += f" (Chain {chain_b})"
    # axes[0,1].set_title(title_b)

    # plot distance difference map of structure 1 - structure 2, marked as A-B/1-2
    # A-B/1-2

    # Here we use TwoSlopeNorm instead of LinearNorm to ensure that 0 is at the center of the color bar
    # Customize color mapping rules for zero-bounded difference data
    norm = TwoSlopeNorm(vmin=-vdiff_use, vcenter=0, vmax=vdiff_use)
    imAB = axes[1,0].imshow(
        diff_ab,
        cmap=cmap_diff,
        norm=norm,
        interpolation="nearest",
        origin="upper"
    )
    # axes[1,0].set_title("A - B (A->B: Red = closer, blue = farther)")

    # plot distance difference map of structure 2 - structure 1, marked as B-A/2-1
    # B-A/2-1
    imBA = axes[1,1].imshow(
        diff_ba,
        cmap=cmap_diff,
        norm=norm,
        interpolation="nearest",
        origin="upper"
    )
    # axes[1,1].set_title("B - A (B->A: Red = closer, blue = farther)")


    # The following annotated function is used to draw the starting and ending points of a region on the axis
    '''
    # annotate the region of interest:
    # '#FFA500'(orange) for region_1, '#800080'(purple) for region_2, avoid using red/blue, which are used in the distance/difference maps
    def _mark_regions(ax, r1, r2=None):
        r1s, r1e = r1
        for pos in (r1s, r1e):
            ax.axhline(pos, color='#FFA500', linestyle='--', linewidth=1.2, zorder=5)
            ax.axvline(pos, color='#FFA500', linestyle='--', linewidth=1.2, zorder=5)
        if r2 is not None:
            r2s, r2e = r2
            for pos in (r2s, r2e):
                ax.axhline(pos, color='#800080', linestyle='--', linewidth=1.2, zorder=5)
                ax.axvline(pos, color='#800080', linestyle='--', linewidth=1.2, zorder=5)
    '''

    # Mark region pairs function
    # annotate the region of interest using Rectangle patch:
    def _mark_pairs(ax, region_pairs, color="#00FF00"):
        for (ra, rb) in region_pairs:
            a_s, a_e = ra
            b_s, b_e = rb
            ha = a_e - a_s + 1
            hb = b_e - b_s + 1
            # draw ra x rb rectangle (rows=ra, cols=rb)
            ax.add_patch(
                Rectangle(
                    (b_s, a_s), hb, ha,
                    fill=False,
                    edgecolor=color,
                    linestyle="--",
                    linewidth=1.5,
                    alpha=1.0,
                    zorder=5
                )
            )
            # draw rb x ra rectangle (rows=rb, cols=ra) - symmetric block
            ax.add_patch(
                Rectangle(
                    (a_s, b_s), ha, hb,
                    fill=False,
                    edgecolor=color,
                    linestyle="--",
                    linewidth=1.5,
                    alpha=1.0,
                    zorder=5
                )
            )

    # Chain Bar Drawing Function
    '''
    def _add_chain_bars(ax, boundaries):
        divider = make_axes_locatable(ax)
        
        # Colors
        n_chains = len(boundaries)
        if n_chains <= 20:
            cmap = plt.get_cmap("tab20", n_chains)
        else:
            cmap = ListedColormap(plt.get_cmap("gist_rainbow")(np.linspace(0.1, 0.9, n_chains)))


        # Create index array
        total_len = boundaries[-1][2] + 1
        chain_row = np.zeros((1, total_len))
        for i, (_, s, e) in enumerate(boundaries):
            chain_row[0, s:e+1] = i
            
        # Top Bar
        ax_top = divider.append_axes("top", size="5%", pad=0.05)
        ax_top.imshow(chain_row, cmap=cmap, aspect="auto", alpha=0.8)
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        
        # Left Bar
        ax_left = divider.append_axes("left", size="5%", pad=0.05)
        ax_left.imshow(chain_row.T, cmap=cmap, aspect="auto", alpha=0.8)
        ax_left.set_xticks([])
        ax_left.set_yticks([])
        
        # Labels and Separators
        for i, (cid, s, e) in enumerate(boundaries):
            center = (s + e) / 2
            # Top Label
            ax_top.text(center, 0, cid, ha="center", va="center", fontsize=9, fontweight='bold')
            # Left Label
            ax_left.text(0, center, cid, ha="center", va="center", fontsize=9, fontweight='bold', rotation=90)
            
            # Separators (skip for first chain start)
            if i > 0:
                sep = s - 0.5
                # On bars
                ax_top.axvline(sep, color='w', linewidth=1)
                ax_left.axhline(sep, color='w', linewidth=1)
                # On main plot
                ax.axvline(sep, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
                ax.axhline(sep, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
    '''


    # Helper: Add Chain Bars & Title
    def _add_chain_bars_and_title(ax, boundaries, title_text):
        divider = make_axes_locatable(ax)
        n_chains = len(boundaries)
        if n_chains <= 20:
            cmap = plt.get_cmap("tab20", n_chains)
        else:
            cmap = ListedColormap(plt.get_cmap("gist_rainbow")(np.linspace(0.1, 0.9, n_chains)))

        total_len = boundaries[-1][2] + 1
        chain_row = np.zeros((1, total_len))
        for i, (_, s, e) in enumerate(boundaries):
            chain_row[0, s:e+1] = i
            
        # Top Bar
        ax_top = divider.append_axes("top", size="5%", pad=0.05)
        ax_top.imshow(chain_row, cmap=cmap, aspect="auto", alpha=0.8)
        ax_top.set_xticks([]); ax_top.set_yticks([])
        
        # ⚠️ Set title on the top bar axis to avoid overlap
        ax_top.set_title(title_text, fontsize=12, pad=8, fontweight='bold')

        # Left Bar
        ax_left = divider.append_axes("left", size="5%", pad=0.05)
        ax_left.imshow(chain_row.T, cmap=cmap, aspect="auto", alpha=0.8)
        ax_left.set_xticks([]); ax_left.set_yticks([])
        
        for i, (cid, s, e) in enumerate(boundaries):
            center = (s + e) / 2
            ax_top.text(center, 0, cid, ha="center", va="center", fontsize=9, fontweight='bold')
            ax_left.text(0, center, cid, ha="center", va="center", fontsize=9, fontweight='bold', rotation=90)
            if i > 0:
                sep = s - 0.5
                ax_top.axvline(sep, color='w', linewidth=1)
                ax_left.axhline(sep, color='w', linewidth=1)
                ax.axvline(sep, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
                ax.axhline(sep, color='k', linestyle='--', linewidth=0.5, alpha=0.5)


    # Helper: Set Custom Ticks (Per Chain 0-based)
    def _set_custom_ticks(ax, boundaries, step=100):
        locs = []
        labels = []
        for _, s, e in boundaries:
            # Always label start (0)
            locs.append(s)
            labels.append("0")
            # Label increments
            curr = s + step
            while curr <= e:
                locs.append(curr)
                labels.append(str(curr - s))
                curr += step
        
        ax.set_xticks(locs)
        ax.set_xticklabels(labels, fontsize=8)
        # Hide Y ticks as requested
        ax.set_yticks([]) 

    
    # Plotting Loop
    for i, ax in enumerate(axes.ravel()):
        # Only mark pairs if NOT default full length
        if not is_default_region:
            _mark_pairs(ax, pairs)
        
        # Determine boundaries for this subplot
        # Row 0 Col 0 (A) -> boundaries_a
        # Row 0 Col 1 (B) -> boundaries_b
        # Row 1 (Diff) -> boundaries_a (Assuming alignment is perfect, A's boundaries apply)
        
        current_boundaries = boundaries_b if ax == axes[0, 1] else boundaries_a
        
        # Add bars and pass title to be set on top axis
        _add_chain_bars_and_title(ax, current_boundaries, titles[i])
        
        # Set custom ticks
        _set_custom_ticks(ax, current_boundaries, tick_step)

        # Set limits (important to do after bars to ensure alignment)
        N = current_boundaries[-1][2] + 1
        ax.set_xlim(0, N-1)
        ax.set_ylim(N-1, 0)
        
        # set labels , note that residue index is 0-based here
        ax.set_xlabel("Residue index")
        # remove y label to avoid overlap
        # ax.set_ylabel("Residue index")

    # color bars are shared between subplots in the same row, cause we are comparing the same protein in structure1 and structure2
    cbar_top = fig.colorbar(imA, ax=[axes[0,0], axes[0,1]], fraction=0.046, pad=0.02)
    cbar_top.set_label("Cα–Cα distance (Å)")
    cbar_bot = fig.colorbar(imAB, ax=[axes[1,0], axes[1,1]], fraction=0.046, pad=0.02)
    cbar_bot.set_label("Δ distance (Å)")

    # save the figure in pdf format
    # or we can save it in png format if needed, set dpi=300 or higher
    if out_path:
        # save figure
        sel_name = "_".join(region_pairs) if region_pairs else f"{region_1}" if region_2 is None else f"{region_1}_vs_{region_2}"
        plt.savefig(f"{out_path}/{job_name}_{sel_name}.pdf", bbox_inches='tight')
        # also save as png, dpi 300
        plt.savefig(f"{out_path}/{job_name}_{sel_name}.png", bbox_inches='tight', dpi=300)
        plt.close(fig)

    
    # return the contact matrices if needed
    # first, we summary the current result 
    result = {
        "A":{
            "file": mmcif_file_a, # file path
            "N_res": Na, # number of residues of whole protein A
            "dist_matrix": dist_a, # array/matrix of pairwise distance of protein A
            "res_ca_info": ca_info_a # list of dicts of CA atom info for whole protein A
        },
        "B":{
            "file": mmcif_file_b,
            "N_res": Nb,
            "dist_matrix": dist_b,
            "res_ca_info": ca_info_b
        },
        "Dist_diff":{
            "A-B": diff_ab, # array/matrix of distance difference A - B
            "B-A": diff_ba # array/matrix of distance difference B - A
        },
        "pairs": pairs # list of region pairs compared
    }

    # then, for each region pair, we extract the sub-matrix and related info
    details = []
    for (ra, rb) in pairs:
        a_s, a_e = ra
        b_s, b_e = rb
        sub_a = dist_a[a_s:a_e+1, b_s:b_e+1]
        sub_b = dist_b[a_s:a_e+1, b_s:b_e+1]
        sub_ab = diff_ab[a_s:a_e+1, b_s:b_e+1]
        sub_ba = diff_ba[b_s:b_e+1, a_s:a_e+1]
        if sub_a.size == 0 or sub_b.size == 0:
            raise ValueError(f"The selected regions {ra} or {rb} are invalid, resulting in empty sub-matrix.")
        details.append({
            "pair": (ra, rb), # region pair
            "A_sub": sub_a, # sub-matrix of region ra vs region rb in structure A
            "B_sub": sub_b, # sub-matrix of region ra vs region rb in structure B
            "A-B_sub": sub_ab, # sub-matrix of region ra vs region rb in distance difference A - B
            "B-A_sub": sub_ba, # sub-matrix of region ra vs region rb in distance difference B - A
            "shape": sub_a.shape, # shape of the sub-matrix
            "mean_A_sub": float(np.nanmean(sub_a)), # mean value of the sub-matrix in structure A
            "mean_B_sub": float(np.nanmean(sub_b)), # mean value of the sub-matrix in structure B
            "mean_A-B_sub": float(np.nanmean(sub_ab)), # mean value of the sub-matrix in distance difference A - B
            "mean_B-A_sub": float(np.nanmean(sub_ba)) # mean value of the sub-matrix in distance difference B - A
        })
    result["Pairs_detail"] = details

    # finally, we return the result if needed
    if return_maxtrix:
        return result



