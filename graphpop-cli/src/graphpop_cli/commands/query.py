"""graphpop query — run arbitrary Cypher and format as TSV."""
import click
from ..cli import pass_ctx
from ..formatters import format_output


@click.command()
@click.argument("cypher")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def query(ctx, cypher, output_path, fmt):
    """Run an arbitrary Cypher statement and format the results."""
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "query", {"cypher": cypher})
