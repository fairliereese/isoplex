## Input Data Format

Both the CLI and Python API work with the **same underlying table structure**.  
The CLI reads from TSV/CSV files, while the API can take a `pandas.DataFrame` directly.

Every table must include:

1. **Gene column**  
   Contains the gene identifier for each row.  
   The column name can be anything; it is specified when running the tool.

2. **Feature column**  
   Contains the feature identifier (for example, isoform or transcript ID).  
   The column name is also user-specified.

3. **Expression data**  
   One or more numeric columns containing either raw **counts** or normalized **TPM** values.  
   The type (counts vs TPM) is specified at runtime.

---

### Single-Sample Format
For datasets representing **one sample**, include exactly **one numeric expression column**.

Example:

```text
gene_id    transcript_id    sample1
G1         T1               100
G1         T2                50
G2         T3                30
G2         T4                 0
```

### Multi-Sample Format

For datasets with **multiple samples**, include **one column per sample**, and all must be numeric.

Example:

```text
gene_id    transcript_id    sample1    sample2    sample3
G1         T1               100        80         90
G1         T2                50        20         30
G2         T3                10         5          8
G2         T4                 0         2          4
```

### Notes

Column names for gene, feature, and samples are flexible — you specify them when using either the CLI or the API.

All expression columns must contain only non-negative, numeric values.

For API usage, simply pass a properly structured pandas.DataFrame; no file format or separator concerns apply.
