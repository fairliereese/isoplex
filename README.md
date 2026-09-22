# Isoform Perplexity

![PyPI version](https://img.shields.io/pypi/v/isoplex.svg)

<!-- [![Documentation Status](https://readthedocs.org/projects/isoplex/badge/?version=latest)](https://isoplex.readthedocs.io/en/latest/?version=latest) -->

Compute and filter isoforms based on perplexity

* PyPI package: https://pypi.org/project/isoplex/
* Free software: MIT License
* Documentation: https://fairliereese.github.io/isoplex/

## Features

This library implements the basic computations described in our [Genome Biology paper](https://doi.org/10.1186/s13059-026-04260-2) (Schertzer et al., 2026). In brief, from a counts or TPM expression matrix of transcriptome features (transcript isoforms, ORF ids, protein IDs, etc) and their associated genes, `isoplex` will compute the gene potential, entropy, perplexity, and mark effective features based on the aforementioned metrics.

These metrics are designed to provide a more intuitive description of isoform diversity as well as to provide a less rigid framework for filtering isoforms, as the distribution of expression values are taken into account for each gene to mark isoforms as effective or not, rather than applying a uniform filter across the entire dataset.

## Citation

If you use `isoplex`, please cite:

Schertzer MD\*, Park SH\*, Su J, Reese F, Sheynkman GM, Knowles DA. Perplexity as a metric for isoform diversity in the human transcriptome. *Genome Biology*. 2026. https://doi.org/10.1186/s13059-026-04260-2

## Credits

This package was created with [Cookiecutter](https://github.com/audreyfeldroy/cookiecutter) and the [audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage) project template.
