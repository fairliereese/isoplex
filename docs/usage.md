# Usage

`isoplex` allows you to perform the computations either from within Python, or using the command line. There are two modes: *single-sample* (global) or *multi-sample*. If you have a multi-sample expression matrix, please use the corresponding option; global-level metrics will also be provided!

See also the [input data format specifications](input_data.md).

## Python

Compute single-sample statistics:

```python
import pandas as pd
import isoplex

df = pd.read_csv('my_isoform_expression_matrix.tsv', sep='\t')

global_perplexity = isoplex.compute_global_isoform_metrics(df,
                              gene_col='gene_id',
                              feature_col='transcript_id',
                              expression_type='tpm')
```

Compute multi-sample statistics, as well as global statistics:

```python
import pandas as pd
import isoplex

df = pd.read_csv('my_isoform_expression_matrix.tsv', sep='\t')

sample_perplexity, global_perplexity = isoplex.compute_multi_sample_isoform_metrics(df,
                              gene_col='gene_id',
                              feature_col='transcript_id',
                              expression_type='tpm')
```

See more details about the Python functions in the [API docs](api.md).

## CLI

Compute single-sample statistics:

```text
Usage: isoplex global-metrics [OPTIONS] INPUT_FILE OUTPUT_FILE

Compute global isoform (or other feature) diversity metrics for a single-sample dataset.

Arguments:
  INPUT_FILE        Filename to input expression table (CSV or TSV).  [required]
  OUTPUT_FILE       Filename to save the output file.                [required]

Options:
  --gene-col TEXT          Column name for gene IDs. [default: gene_id]
  --feature-col TEXT       Column name for isoform/feature IDs. [default: transcript_id]
  --expression-type TEXT   Expression type in table: 'counts' or 'tpm'. [default: counts]
  --sep TEXT               Delimiter for input/output files. Use ',' for CSV. [default: \t]
  --help                   Show this message and exit.
```

Compute multi-sample statistics, as well as global statistics:

```text
Usage: isoplex multi-sample-metrics [OPTIONS] INPUT_FILE OUTPUT_SAMPLE_FILE
                                    OUTPUT_GLOBAL_FILE

Compute sample-level and global isoform (or other feature) diversity
metrics for a single-sample dataset.

Arguments:
  INPUT_FILE               Filename to input expression table (CSV or TSV).
  OUTPUT_SAMPLE_FILE       Filename to save the sample-level output file.
  OUTPUT_GLOBAL_FILE       Filename to save the global-level output file.

Options:
  --gene-col TEXT          Column name for gene IDs. [default: gene_id]
  --feature-col TEXT       Column name for isoform/feature IDs. [default:
                           transcript_id]
  --expression-type TEXT   Expression type in table: 'counts' or 'tpm'.
                           [default: counts]
  --sep TEXT               Delimiter for input/output files. Use ',' for CSV.
                           [default: \t]
  --help                   Show this message and exit.
```
