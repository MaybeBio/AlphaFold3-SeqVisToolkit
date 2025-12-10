# the function defined here is modified to support 1D track data visualization alongside 2D matrix
from typing import Literal, Optional, List, Union, Dict, Any
import os
import pandas as pd
from matplotlib.colors import ListedColormap, to_hex
import matplotlib.pyplot as plt

# *** Track based grouping is better than Chain based grouping ***
# -1. easy for plotting (track is a layer)
# -2. config is setting for track instead of chain

# 1, transfer the bed data into what we need for 1D track plotting
def parse_bed_to_track_data(
        bed_file: str,
        color: Optional[Union[str, Dict[str, Any]]] = "tab10"
) -> List[Dict[str, Any]]:
    """
    Description
    ----------
    Parse a BED file to extract 1D track data for visualization.
    Helper function to convert a BED-like file into the track_data dictionary format.

    Args
    -----

    bed_file (str): Path to the BED-like file (csv file more precisly) containing 1D track data.y
    color (Union[str, Dict[str, str]]): Color for the track plot. 
        - String: A colormap name (e.g. "tab10") or a single color (e.g. "orange").
        - Dict: {TrackName: ColorConfig}.
          e.g. {"IDR": "red", "Domain": {"DomainA": "blue", "DomainB": "green"}}
          or {"IDR": "red", "Domain": "tab10"}

    Returns
    ------
    List[Dict[str, Any]]: A list of track data dictionaries for visualization. 
    Each track dictionary corresponds to a unique track in the BED file.

    Notes
    -----
    - 1, The track bed file must be 0-based indexed!
    
    """

    # First, we load the bed-like file(tsv)
    track_df = pd.read_csv(bed_file, sep="\t", header=0)   

    # here, we assume the bed file has 5 columns: chain_id, type, start, end, value
    # col0=Chain, col1=Type, col2=Start, col3=End, col4=Value
    tracks_to_return = []

    # --- Helper: Process a single track group ---
    def _process_group(track_df, track_name, track_color, track_val_type):
        """
        Description
        -----------
        Process a single dataframe group into a track dictionary.

        Args
        ----
            track_df (pd.DataFrame): DataFrame containing the track data for a specific track.
            track_name (str): Name of the track.
            track_color (Union[str, Dict[str, str]]): Color configuration for the track.
            track_val_type (str): Type of the value in the "value" column, either "numerical" or "categorical".

        Notes
        -----
        - 1, Here, we assumen the bed file has 5 columns: chain_id, type, start, end, value
        - 2, For value column, we assume only two types: numerical (float or int) and categorical (str).
        Here we treat all str values as categorical values.
        If there is string that is not categorical, we will deal with it further.

        Todos
        -----
        - 1, Add support for more complex value types if needed in the future.
        If there is string that is not categorical, we can add a new parameter "is_categorical" to specify it.
        So string will be grouped into 2 types: categorical and non-categorical based on parameter "is_categorical".

        """
        track_dict = {}
        all_categories = set()

        # Group by Chain ID
        for chain_id, group in track_df.groupby(track_df.columns[0]):
            # here we use the second strategy
            chain_vals = {}
            
            # note that group is still a dataframe
            # and for dataframe we iterate rows using iterrows()
            for _, row in group.iterrows():
                # note that _ is the index, row is the actual row data, index is useless here
                start = int(row.iloc[2])
                end = int(row.iloc[3])
                raw_val = row.iloc[4]

                # numerical value, like disorder score
                if track_val_type == "numerical":
                    # actually float is enough for int values
                    value = float(raw_val)
                else:
                    # str type, categorical value, like hydrophobic, polar, charged
                    value = str(raw_val)
                    # NOTE AGAIN HERE! we treat all str values as categorical values
                    all_categories.add(value)

                # then we fill in the track_dict
                for i in range(start, end + 1):
                    # we assume here the bed range is closure [start, end]
                    chain_vals[i] = value
            
            # finally, we add this chain's data to track_dict
            # Now we have {res: bed value} -> chain_id -> track 
            track_dict[chain_id] = chain_vals 
        
        # For color settings, we will deal with it further
        final_color = track_color
        if track_val_type == "str":
            # For categorical data, we need to generate color map {category: color}
            # first, we will see how many categories we have
            unique_cats = sorted(list(all_categories))

            # Case1: user provide a dict for this track, assume it's {Category: Color}
            if isinstance(track_color, dict):
                final_color = track_color
            
            # Case2: user provide a colormap name or None
            # we need to generate color map ourselves
            else:
                # if colormap, we accept it; if None, we use default "tab10"(fallback)
                cmap_name = track_color if isinstance(track_color, str) else "tab10"
                try:
                    cmap = plt.get_cmap(cmap_name, len(unique_cats))
                    # then we generate {category:hex color} dict
                    final_color = {
                        cat: to_hex(cmap(i)) for i, cat in enumerate(unique_cats)
                    }
                # if the colormap name is invalid, we just use default colors, like "tab10"
                except ValueError:
                    print(f"Warning: Colormap '{cmap_name}' not found. Using 'tab10'.")
                    cmap = plt.get_cmap("tab10", len(unique_cats))
                    final_color = {cat: to_hex(cmap(i)) for i, cat in enumerate(unique_cats)}
        
        # ! for numerical data?


        # Finally, we construct the track_data dictionary
        return {
            "track_name": track_name,
            "track_type": track_val_type,
            "color": final_color,
            "track_data": track_dict
        }

    # ------- Main Logic: Group by Track Name (Column 1) -------

    # group by the type column (track name, col1)
    # Note again: we assume that the bed file has 5 columns: chain_id, type, start, end, value
    
    # number of unique tracks
    groups = list(track_df.groupby(track_df.columns[1]))

    # prepare default palette for numerical tracks (to distinguish them if no config provided)
    default_palette = plt.get_cmap("tab10", len(groups))

    for i, (track_name, track_group) in enumerate(groups):

        track_name = str(track_name)

        # 1. Auto-Detect Type (numerical or categorical) if needed
        try:
            pd.to_numeric(track_group.iloc[:, 4], errors="raise")
            val_type = "numerical"
        except ValueError:
            val_type = "categorical"

        # 2. Determine Color Config for THIS track
        # Priotity: Dict Config > Colormap Generation > Single Color
        this_track_color = None

        # if color is a dict
        if isinstance(color, dict):
            # check if there is a config for this Track Name
            this_track_color = color.get(track_name)

        # Usually, we provide colormap dict for well-annotated categorical tracks
        # for example, for charge type, we can provide {"Positive": "red", "Negative": "blue", "Neutral": "gray"}
        # cause charge is well-annotated categorical track, we know what categories are there, so we can provide color dict directly and define each category color
        # but there are some categorical tracks that are not well-annotated, for example, domain types, we may not know how many domain types are there in advance
        # ⚠️ WHAT I MEAN HERE IS THAT: if this_track_color is None, there is definitely no chance that this track is a well-annotated categorical track
        # so there are two possibilities:
        # 1, this is a numerical track, we can just generate distinct colors for distinct tracks
        # 2, this is a categorical track, but not well-annotated, we can only provide a colormap name or None, then we generate color map based on the categories
        if this_track_color is None:
            # Fallback to string/palette
            if isinstance(color, str):
                if val_type == "categorical":
                    this_track_color = color  # colormap name for categorical
                else:
                    # For numerical, we want distinct colors for distinct tracks
                    if color in plt.colormaps():
                        this_track_color = to_hex(default_palette(i))
                    else:
                        this_track_color = color  # User forced single color "orange"
            else:
                this_track_color = to_hex(default_palette(i))
        
        # 3. Process this track group
        track_obj = _process_group(
            track_group,
            track_name,
            this_track_color,
            val_type
        )
        tracks_to_return.append(track_obj)
    
    return tracks_to_return

        