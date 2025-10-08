# Welcome to isoplex's documentation!

This library implements the basic computations described by Schertzer et al. in their [manuscript](https://www.biorxiv.org/content/10.1101/2025.07.02.662769v1). In brief, from a counts or TPM expression matrix of transcriptome features (transcript isoforms, ORF ids, protein IDs, etc) and their associated genes, `isoplex` will compute the gene potential, entropy, perplexity, and mark effective features based on the aforementioned metrics.

These metrics are designed to provide a more intuitive description of isoform diversity as well as to provide a less rigid framework for filtering isoforms, as the distribution of expression values are taken into account for each gene to mark isoforms as effective or not, rather than applying a uniform filter across the entire dataset.


## Contents

- [Installation](installation.md)
- [Usage](usage.md)
- [Input data formatting](input_data.md)
- [API](api.md)
- [Contributing](contributing.md)
