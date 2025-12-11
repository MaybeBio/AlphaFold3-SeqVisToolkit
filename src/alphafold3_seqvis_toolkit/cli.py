import typer # For building CLI applications, command-line interfaces
from typing import List, Optional, Dict, Any
import os
from alphafold3_seqvis_toolkit.modules.confidence_metrics_plot import plot_local_confidence, plot_global_confidence
from alphafold3_seqvis_toolkit.modules.contact_map_comparison_monomer import contact_map_diff_monomer
from alphafold3_seqvis_toolkit.modules.contact_map_comparison_multimer import contact_map_diff_multimer
from alphafold3_seqvis_toolkit.modules.contact_map_visualization_with_track import contact_map_vis_with_track
from alphafold3_seqvis_toolkit.modules.contact_map_visualization_without_track import contact_map_vis_without_track

app = typer.Typer(help="AlphaFold3 SeqVis Toolkit", no_args_is_help=True)

@app.command(
    "confidence",
    no_args_is_help=True
)
def confidence_cmd(
    global_json: Optional[str] = typer.Option(None, "--global-json", help="fold_{YOUR-JOB-NAME}_summary_confidences_{i}.json (Required for global/all mode)", rich_help_panel="Input"),
    full_json: Optional[str] = typer.Option(None, "--full-json", help="fold_{YOUR-JOB-NAME}_full_data_{i}.json (Required for local/all mode)", rich_help_panel="Input"),
    output_path: str = typer.Option(".", "--output-path", "-o", help="Directory for outputs, default is current directory", rich_help_panel="Output"),
    mode: str = typer.Option("all", "--mode", "-m", help="Analysis mode: 'all' (default), 'global', or 'local'.", rich_help_panel="Mode"),
    chains: Optional[List[str]] = typer.Option(None, "--chains", "-c", help="Repeatable: chain IDs for local subset. Only used in 'local' or 'all' mode with --full-json.", rich_help_panel="Local Plot Options"),
    tick_step: int = typer.Option(100, "--tick-step", help="Residue tick step for local plots", rich_help_panel="Local Plot Options"),
):
    """
    Plot global (ipTM/pTM etc.) and/or local (PAE/contact/atom pLDDT) confidence metrics.

    \b
    Notes:
    - 1, Use --mode to control the scope:
        - 'all' (default): Plots everything possible based on provided JSON files. If only one file is provided, it plots that one.
        - 'global': Only plots global metrics (requires --global-json).
        - 'local': Only plots local metrics (requires --full-json).
    - 2, [--chains] option is only effective for LOCAL metrics (PAE, contact probs, atom pLDDT).
    - 3, There may be NULL values in the above metrics produced by AlphaFold3, so we will convert them to NaN (these values may appear as NA when output, and this handling also applies to plotting). Therefore, if you are confused about the output, it is recommended to first check your original data.
    - 4, All residue indices in this module are 0-based logic driven.    
    
    \b
    Examples:
    - 1, Plot EVERYTHING (Global + Local) for a job:
    af3-vis confidence \
--global-json summary_confidences.json \
--full-json full_data.json \
-o out_path
    - 2, Only plot Global metrics:
    af3-vis confidence --mode global --global-json summary_confidences.json -o out_path
    - 3, Only plot Local metrics for specific chains:
    af3-vis confidence --mode local --full-json full_data.json -c A -c B -o out_path
    """

    # Basic validation
    if mode == "all" and not global_json and not full_json:
        raise typer.Exit("Error: In 'all' mode, need at least one of --global-json or --full-json.")

    import os
    os.makedirs(output_path, exist_ok=True)

    # Logic for Global
    if mode in ["all", "global"]:
        if global_json:
            plot_global_confidence(confid_json_file_path=global_json, output_path=output_path)
        elif mode == "global":
            raise typer.Exit("Error: --global-json is required when mode is 'global'.")

    # Logic for Local
    if mode in ["all", "local"]:
        if full_json:
            plot_local_confidence(
                full_json_file_path=full_json,
                output_path=output_path,
                chains=chains,
                tick_step=tick_step,
            )
        elif mode == "local":
            raise typer.Exit("Error: --full-json is required when mode is 'local'.")


@app.command(
    "contact-map-diff",
    no_args_is_help=True
)
def contact_map_diff_cmd(
    mmcif_a: str = typer.Option(..., "--mmcif-a", help="Path to mmCIF file A", rich_help_panel="Inputs"),
    mmcif_b: str = typer.Option(..., "--mmcif-b", help="Path to mmCIF file B", rich_help_panel="Inputs"),
    chain_a: Optional[str] = typer.Option(..., "--chain-a", help="Chain ID(s) for mmCIF file A. Monomer mode: single chain (e.g. 'A'). Multimer mode: one or more chains (e.g. 'A' or 'A,B').", rich_help_panel="Inputs"),
    chain_b: Optional[str] = typer.Option(..., "--chain-b", help="Chain ID(s) for mmCIF file B. Must align with chain-a (⚠️  same number and sequence!).", rich_help_panel="Inputs"),
    region_1: Optional[str] = typer.Option(None, "--region-1", help="(Legacy) Select a region to focus on/compare (highlighted with a green box). In multimer mode, supports 'Chain:Start:End'.", rich_help_panel="Legacy"),
    region_2: Optional[str] = typer.Option(None, "--region-2", help="(Legacy) Select a second region to compare against region-1.", rich_help_panel="Legacy"),
    region_pair: List[str] = typer.Option(None, "--region-pair", help="Repeatable region pair(s) selected to focus on/compare (highlighted with green boxes): 'Start:End,Start:End' format like 'a:b,c:d' or 'a-b,c-d' for monomer mode, Chains-specified 'Chain:Start:End' format like 'A:a:b,B:c:d' or 'A:a-b,B:c-d' for multimer mode. Use multiple --region-pair to add more.", rich_help_panel="Regions"),
    vmax: Optional[float] = typer.Option(None, help="vmax for distance heatmap", rich_help_panel="Color scaling"),
    vmax_percentile: float = typer.Option(95.0, help="Percentile used if vmax is not set", rich_help_panel="Color scaling"),
    vdiff: Optional[float] = typer.Option(None, help="Max abs value for diff heatmap (0-centered)", rich_help_panel="Color scaling"),
    vdiff_percentile: float = typer.Option(95.0, help="Percentile used if vdiff is not set", rich_help_panel="Color scaling"),
    include_nonstandard_residue: bool = typer.Option(False, "--include-nonstandard/--no-include-nonstandard", help="Include non-standard amino acid residues"),
    out_path: Optional[str] = typer.Option(".", "--out-path", help="Directory to save figure files (png/pdf), default is current directory", rich_help_panel="Output"),
    mode: str = typer.Option("multimer", "--mode", "-m", help="Comparison mode: 'multimer' (default) or 'monomer'.", rich_help_panel="Mode"),
    tick_step: int = typer.Option(100, "--tick-step", help="Step size for ticks on the axes (multimer mode only)", rich_help_panel="Plot Options"),
):
    """
    Compare contact maps between two AlphaFold3/General mmCIF structures (for the same molecule with an identical sequence), and plot the distance/diff matrices. Supports both monomer and multimer modes.

    \b
    Notes:
    - 1, Designed for comparing the same protein sequence under different conditions.
    - 2, Regions are 0-based and inclusive: 'start:end' means [start, end].
    - 3, If region_2 is omitted, region_1 is applied symmetrically.
    - 4, We strongly recommend using --region-pair for multiple region comparisons, as it is more flexible and clearer. 
    And you can also use --region-pair to replace the legacy --region-1/--region-2 options. 
    - 5, Use --mode to switch between 'monomer' and 'multimer' (default) comparison logic.
    - 6, All residue indices in this module are 0-based logic driven.
    - 7, chain_a and chain_b should strictly align in order with same sequence! E.g., in multimer mode, if chain_a is "A,B,C", chain_b is "D,E,F", then it should be A aligns with D, B aligns with E, and C aligns with F.

    \b
    Examples:
    \b
    - 1, Basic diff on one region (Multimer default):
    af3-vis contact-map-diff \
--mmcif-a A.cif --mmcif-b B.cif \
--region-1 0:200 --chain-a A --chain-b A --out-path . 
    \b
    - 2, Diff between two regions (Monomer mode):
    af3-vis contact-map-diff \
--mmcif-a A.cif --mmcif-b B.cif \
--region-1 0-200 --region-2 300-500 \
--chain-a A --chain-b A --mode monomer \
--out-path .
    \b
    - 3, Multiple region pairs (Multimer, Chain A vs Chain A): 
    af3-vis contact-map-diff \
--mmcif-a A.cif --mmcif-b B.cif \
--region-pair 265:576,0:15 \
--region-pair 265:576,42:54 \
--region-pair 265:576,83:93 \
--region-pair 265:576,170:189 \
--region-pair 265:576,214:241 \
--region-pair 265:576,578:589 \
--region-pair 265:576,606:639 \
--region-pair 265:576,691:720 \
--chain-a A --chain-b A \
--out-path .
    \b
    - 4, Multimer Complex (Chain A,B vs Chain A,B) with specific chain regions: (🌟 RECOMMENDED IN ANY CASE!)
    af3-vis contact-map-diff \
--mmcif-a complex_v1.cif --mmcif-b complex_v2.cif \
--chain-a A,B --chain-b A,B \
--region-pair A:10:50,B:10:50 \
--region-pair A:100:150,A:100:150 \
--out-path .
    """
    
    if mode == "monomer":
        contact_map_diff_monomer(
            mmcif_file_a=mmcif_a,
            mmcif_file_b=mmcif_b,
            chain_a=chain_a,
            chain_b=chain_b,
            region_1=region_1,
            region_2=region_2,
            region_pairs=region_pair if region_pair else None,
            vmax=vmax,
            vmax_percentile=vmax_percentile,
            vdiff=vdiff,
            vdiff_percentile=vdiff_percentile,
            include_nonstandard_residue=include_nonstandard_residue,
            out_path=out_path,
        )
    else:  # multimer mode
        contact_map_diff_multimer(
            mmcif_file_a=mmcif_a,
            mmcif_file_b=mmcif_b,
            chain_a=chain_a,
            chain_b=chain_b,
            region_1=region_1,
            region_2=region_2,
            region_pairs=region_pair if region_pair else None,
            vmax=vmax,
            vmax_percentile=vmax_percentile,
            vdiff=vdiff,
            vdiff_percentile=vdiff_percentile,
            include_nonstandard_residue=include_nonstandard_residue,
            out_path=out_path,
            tick_step=tick_step,
        )

@app.command(
    "contact-map-vis",
    no_args_is_help=True
)
def contact_map_vis_cmd(
    mmcif_file: str = typer.Option(..., "--mmcif-file", help="Path to mmCIF file", rich_help_panel="Input"),
    chains: Optional[List[str]] = typer.Option(None, "--chains", "-c", help="Repeatable: chain IDs to include, default is all chains", rich_help_panel="Input"),
    out_path: str = typer.Option(".", "--out-path", "-o", help="Directory for outputs, default is current directory", rich_help_panel="Output"),
    mode: str = typer.Option("no-track", "--mode", "-m", help="Visualization mode: 'no-track' (default) or 'track'.", rich_help_panel="Mode"),
    track_bed_file: Optional[str] = typer.Option(None, "--track-bed-file", help="Path to BED file for custom tracks, e.g., domains、IDRs (Required if mode is 'track')", rich_help_panel="Custom Tracks"),
    color_config: Optional[str] = typer.Option("tab10", "--color-config", help="Path to color config file (JSON) or colormap name (Only used if mode is 'track')", rich_help_panel="Custom Tracks"),
    tick_step: int = typer.Option(100, "--tick-step", help="Step size for ticks on the axes", rich_help_panel="Plot Options"),
):
    """
    Visualize contact map from an AlphaFold3 mmCIF structure or a general mmCIF structure. Supports 'no-track' (simple) and 'track' (custom annotation) modes.

    \b
    Notes:
    - 1, Use --mode to switch between 'no-track' (default) and 'track' visualization.
    - 2, 'no-track' mode: Generates a standard contact map.
    - 3, 'track' mode: Generates a contact map with custom annotation tracks (e.g., domains, IDRs).
         - Requires --track-bed-file.
         - Tracks can be numerical (line plot) or categorical (bar/strip plot).
         - The track bed file must be 0-based indexed!
    - 4, By default, all chains in the mmCIF file are included. Use --chains to specify particular chains if needed.
    - 5, A color configuration file can be provided to customize the colors of categorical tracks (only for 'track' mode).
    - 6, Modify the tick_step parameter to adjust the spacing of residue ticks on the axes as needed.

    \b
    Examples:
    - 1, Basic contact map visualization (Default no-track):
    af3-vis contact-map-vis --mmcif-file model.cif -o out_path
    - 2, Contact map with custom annotation tracks:
    af3-vis contact-map-vis \
--mmcif-file model.cif \
--mode track \
--track-bed-file custom_tracks.bed \
--color-config color_config.json \
-o out_path
    """

    if mode == "track":
        if not track_bed_file:
            raise typer.Exit("Error: --track-bed-file is required when --mode is 'track'")
        
        contact_map_vis_with_track(
            mmcif_file=mmcif_file,
            chains=chains,
            out_path=out_path,
            track_bed_file=track_bed_file,
            color_config=color_config,
            tick_step=tick_step,
        )
    else:  # no-track mode
        contact_map_vis_without_track(
            mmcif_file=mmcif_file,
            chains=chains,
            out_path=out_path,
            tick_step=tick_step,
        )

if __name__ == "__main__":
    app()