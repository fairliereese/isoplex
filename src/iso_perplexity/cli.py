"""Console script for iso_perplexity."""

import typer
from rich.console import Console

from iso_perplexity import utils

app = typer.Typer()
console = Console()


@app.command()
def main():
    """Console script for iso_perplexity."""
    console.print("Replace this message by putting your code into "
               "iso_perplexity.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")
    utils.do_something_useful()


if __name__ == "__main__":
    app()
