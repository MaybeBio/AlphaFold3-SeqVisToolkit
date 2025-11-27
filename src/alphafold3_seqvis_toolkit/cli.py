import typer # For building CLI applications, command-line interfaces
from typing import List, Optional
import os
from alphafold3_seqvis_toolkit.modules.confidence_metrics_plot import plot_local_confidence, plot_global_confidence
from alphafold3_seqvis_toolkit.modules.contact_map_comparison import contact_map_diff
from alphafold3_seqvis_toolkit.modules.custom_metric_track import track_metrics

app = typer.Typer(help="AlphaFold3 SeqVis Toolkit", no_args_is_help=True)

@app.command(
    "confidence",
    no_args_is_help=True
)
def confidence_cmd(
    global_json: Optional[str] = typer.Option(None, "--global-json", help="fold_{YOUR-JOB-NAME}_summary_confidences_X.json", rich_help_panel="Input"),
    full_json: Optional[str] = typer.Option(None, "--full-json", help="fold_{YOUR-JOB-NAME}_full_data_X.json", rich_help_panel="Input"),
    output_path: str = typer.Option(".", "--output-path", "-o", help="Directory for outputs", rich_help_panel="Output"),
    chains: Optional[List[str]] = typer.Option(None, "--chains", "-c", help="Repeatable: chain IDs for local subset, this option only works when --full-json is set "
    "and is designed for plotting local metrics (like PAE, contact probs, or atom pLDDT) for specific chains. By default all chains are plotted.", rich_help_panel="Local Plot Options"),
    tick_step: int = typer.Option(100, "--tick-step", help="Residue tick step for local plots", rich_help_panel="Local Plot Options"),
):
    """
    Plot global (ipTM/pTM etc.) and/or local (PAE/contact/atom pLDDT) confidence metrics.

    \b
    Notes:
    - 1, Designed for plotting confidence metrics from AlphaFold3 JSON outputs.
    - 2, [--chains] option is only required in following cases:
        - a, when plotting local metrics from --full-json (input)
        - b, when you want to focus on specific chains in a multi-chain complex ———— this is specifically designed for visualizing PAE/contact/atom pLDDT for selected chains.
        - c, by default (if --chains is not set), all chains are plotted for local metrics.
    - 3, There may be NULL values in the above metrics produced by AlphaFold3, so we will convert them to NaN (these values may appear as NA when output, and this handling also applies to plotting). Therefore, if you are confused about the output, it is recommended to first check your original data.
        
    \b
    Examples:
    - 1, Only plot global metrics:
    af3-vis confidence --global-json fold_MYJOB_summary_confidences_0.json -o out_global_path
    - 2, Only plot local metrics for chains A and B:
    af3-vis confidence --full-json fold_MYJOB_full_data_0.json -o out_local_path -c A -c B
    - 3, Plot both global and local metrics for chain M:
    af3-vis confidence \
--global-json fold_MYJOB_summary_confidences_0.json \
--full-json fold_MYJOB_full_data_0.json \
-o out_all_path -c M
    """
    if not global_json and not full_json:
        raise typer.Exit("Error: need at least one of --global-json or --full-json")

    import os
    os.makedirs(output_path, exist_ok=True)

    if global_json:
        plot_global_confidence(confid_json_file_path=global_json, output_path=output_path)
    if full_json:
        plot_local_confidence(
            full_json_file_path=full_json,
            output_path=output_path,
            chains=chains,
            tick_step=tick_step,
        )





@app.command(
    "contact-map-diff",
    no_args_is_help=True
)
def contact_map_diff_cmd(
    mmcif_a: str = typer.Option(..., "--mmcif-a", help="Path to mmCIF file A", rich_help_panel="Inputs"),
    mmcif_b: str = typer.Option(..., "--mmcif-b", help="Path to mmCIF file B", rich_help_panel="Inputs"),

    region_1: Optional[str] = typer.Option(None, "--region-1", help="(Legacy) single region", rich_help_panel="Legacy"),
    region_2: Optional[str] = typer.Option(None, "--region-2", help="(Legacy) second region", rich_help_panel="Legacy"),
    region_pair: List[str] = typer.Option(None, "--region-pair", help="Repeatable region pair: 'a:b,c:d' or 'a-b,c-d'. Use multiple --region-pair to add more.", rich_help_panel="Regions"),
    vmax: Optional[float] = typer.Option(None, help="vmax for distance heatmap", rich_help_panel="Color scaling"),
    vmax_percentile: float = typer.Option(95.0, help="Percentile used if vmax is not set", rich_help_panel="Color scaling"),
    vdiff: Optional[float] = typer.Option(None, help="Max abs value for diff heatmap (0-centered)", rich_help_panel="Color scaling"),
    vdiff_percentile: float = typer.Option(95.0, help="Percentile used if vdiff is not set", rich_help_panel="Color scaling"),
    include_nonstandard_residue: bool = typer.Option(False, "--include-nonstandard/--no-include-nonstandard", help="Include non-standard amino acid residues"),
    out_file: Optional[str] = typer.Option(None, "--out-file", help="Save figure to file (png/pdf)", rich_help_panel="Output"),
):
    """
    Compare contact maps between two AlphaFold3 mmCIF structures for the same protein sequence, and plot the distance/diff matrices.

    \b
    Notes:
    - 1, Designed for comparing the same protein sequence under different conditions.
    - 2, Regions are 0-based and inclusive: 'start:end' means [start, end].
    - 3, If region_2 is omitted, region_1 is applied symmetrically.
    - 4, We strongly recommend using --region-pair for multiple region comparisons, as it is more flexible and clearer. 
    And you can also use --region-pair to replace the legacy --region-1/--region-2 options. 

    \b
    Examples:
    - 1, Basic diff on one region:
    af3-vis contact-map-diff \
--mmcif-a A.cif --mmcif-b B.cif \
--region-1 0:200 --out-file diff.pdf
    - 2, Diff between two regions:
    af3-vis contact-map-diff \
--mmcif-a A.cif --mmcif-b B.cif \
--region-1 0-200 --region-2 300-500 \
--out-file diff.pdf
    - 3, Multiple region pairs: (RECOMMENDED IN ANY CASE!)
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
--out-file diff.pdf
    """
    
    contact_map_diff(
        mmcif_file_a=mmcif_a,
        mmcif_file_b=mmcif_b,
        region_1=region_1,
        region_2=region_2,
        region_pairs=region_pair if region_pair else None,
        vmax=vmax,
        vmax_percentile=vmax_percentile,
        vdiff=vdiff,
        vdiff_percentile=vdiff_percentile,
        include_nonstandard_residue=include_nonstandard_residue,
        out_file=out_file,
    )


@app.command(
    "track-metrics",
    help="Aggregate metric JSON files in a folder into a single CSV (optional).",
    no_args_is_help=True,
)
def track_metrics_cmd(
    folder: str = typer.Option(..., "--folder", help="Directory of prediction runs", rich_help_panel="Inputs"),
    pattern: str = typer.Option("*.json", "--pattern", help="Glob pattern for metric files"),
    out_csv: Optional[str] = typer.Option(None, "--out-csv", help="Write merged metrics to CSV", rich_help_panel="Output"),
):
    """
    Scan a folder for metric JSON files and merge them for quick tracking.

    Examples:
      af3-vis track-metrics --folder runs/ --pattern 'metrics*.json' --out-csv merged.csv
    """
    data = track_metrics(folder=folder, pattern=pattern, out_csv=out_csv)
    typer.echo(f"Loaded {len(data)} metric records.")

if __name__ == "__main__":
    app()