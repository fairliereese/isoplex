"""Console script for iso_perplexity."""

import typer
from rich.console import Console
from iso_perplexity import utils

app = typer.Typer(help='Command line for running iso_perplexity')
console = Console()


@app.command("global-metrics")
def global_metrics(
    input_file: str = typer.Argument(..., help="Filename to input expression table (CSV or TSV)."),
    output_file: str = typer.Argument(..., help="Filename save the output file."),
    gene_col: str = typer.Option("gene_id", help="Column name for gene IDs."),
    feature_col: str = typer.Option("transcript_id", help="Column name for isoform/feature IDs."),
    expression_type: str = typer.Option("counts", help="Expression type in table: 'counts' or 'tpm'."),
    sep: str = typer.Option("\t", help="Delimiter for input/output files. Use ',' for CSV.")
):
    """
    Compute global isoform (or other feature)-level diversity metrics for a single-sample dataset.
    """

    console.print(f"[bold cyan]Loading input from:[/bold cyan] {input_file}")
    df = pd.read_csv(input_file, sep=sep)

    console.print(f"[bold cyan]Running isoform metrics computation...[/bold cyan]")
    try:
        df_out = compute_global_isoform_metrics(
            df,
            gene_col=gene_col,
            feature_col=feature_col,
            expression_type=expression_type
        )
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(code=1)

    console.print(f"[bold cyan]Saving results to:[/bold cyan] {output_file}")
    df_out.to_csv(output_file, sep=sep, index=False)

    console.print("[bold green]Done![/bold green]")

@app.command("global-metrics")
def multi_sample_metrics(
    input_file: str = typer.Argument(..., help="Filename to input expression table (CSV or TSV)."),
    output_file: str = typer.Argument(..., help="Filename save the output file."),
    gene_col: str = typer.Option("gene_id", help="Column name for gene IDs."),
    feature_col: str = typer.Option("transcript_id", help="Column name for isoform/feature IDs."),
    expression_type: str = typer.Option("counts", help="Expression type in table: 'counts' or 'tpm'."),
    sep: str = typer.Option("\t", help="Delimiter for input/output files. Use ',' for CSV.")
):
    """
    Compute global isoform (or other feature)-level diversity metrics for a single-sample dataset.
    """

    console.print(f"[bold cyan]Loading input from:[/bold cyan] {input_file}")
    df = pd.read_csv(input_file, sep=sep)

    console.print(f"[bold cyan]Running isoform metrics computation...[/bold cyan]")
    try:
        df_out = compute_multi_sample_isoform_metrics(
            df,
            gene_col=gene_col,
            feature_col=feature_col,
            expression_type=expression_type
        )
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(code=1)

    console.print(f"[bold cyan]Saving results to:[/bold cyan] {output_file}")
    df_out.to_csv(output_file, sep=sep, index=False)

    console.print("[bold green]Done![/bold green]")

if __name__ == "__main__":
    app()
