# Usage

`iso-perplexity` allows you to perform the computations either from within Python, or using the command line.

## Python

Compute single-sample statistics:

```python
import pandas as pd
import isoplex

df = pd.read_csv('my_isoform_expression_matrix.tsv', sep='\t')

perplexity = isoplex.compute_global_isoform_metrics(df,
                              gene_col='gene_id',
                              feature_col='transcript_id',
                              expression_type='tpm')
```

Compute multi-sample statistics, as well as global statistics:

```python
import pandas as pd
import isoplex

df = pd.read_csv('my_isoform_expression_matrix.tsv', sep='\t')

perplexity = isoplex.compute_multi_sample_isoform_metrics(df,
                              gene_col='gene_id',
                              feature_col='transcript_id',
                              expression_type='tpm')
```

See more details about the Python functions in the [API docs](api.md).

## CLI

Compute single-sample statistics:

```bash
Usage: isoplex global-metrics [OPTIONS] INPUT_FILE OUTPUT_FILE

Compute global isoform (or other feature) diversity metrics for a single-sample dataset.

╭─ Arguments ────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *    input_file       TEXT  Filename to input expression table (CSV or TSV). [required]                    │
│ *    output_file      TEXT  Filename save the output file. [required]                                      │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --gene-col               TEXT  Column name for gene IDs. [default: gene_id]                                │
│ --feature-col            TEXT  Column name for isoform/feature IDs. [default: transcript_id]               │
│ --expression-type        TEXT  Expression type in table: 'counts' or 'tpm'. [default: counts]              │
│ --sep                    TEXT  Delimiter for input/output files. Use ',' for CSV. [default: \t]            │
│ --help                         Show this message and exit.                                                 │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

Compute multi-sample statistics, as well as global statistics:

```bash
Usage: isoplex multi-sample-metrics [OPTIONS] INPUT_FILE OUTPUT_FILE

Compute sample-level and global isoform (or other feature) diversity metrics for a single-sample dataset.

╭─ Arguments ────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *    input_file       TEXT  Filename to input expression table (CSV or TSV). [required]                    │
│ *    output_file      TEXT  Filename save the output file. [required]                                      │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --gene-col               TEXT  Column name for gene IDs. [default: gene_id]                                │
│ --feature-col            TEXT  Column name for isoform/feature IDs. [default: transcript_id]               │
│ --expression-type        TEXT  Expression type in table: 'counts' or 'tpm'. [default: counts]              │
│ --sep                    TEXT  Delimiter for input/output files. Use ',' for CSV. [default: \t]            │
│ --help                         Show this message and exit.                                                 │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```
