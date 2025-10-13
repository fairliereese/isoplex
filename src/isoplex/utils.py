import pandas as pd
import numpy as np

EXP_COL = 'counts'
FEATURE_COL = TRANSCRIPT_COL = 'transcript_id'
GENE_COL = 'gene_id'
SAMPLE_COL = 'sample'

def validate_counts_input(
    df: pd.DataFrame,
    gene_col: str = GENE_COL,
    feature_col: str = TRANSCRIPT_COL
) -> bool:
    """
    Validate a wide-format expression DataFrame.

    In wide-format input, rows correspond to features (e.g. transcripts),
    columns (other than `gene_col` and `feature_col`) correspond to samples.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format expression table.
    gene_col : str
        Column name identifying genes.
    feature_col : str
        Column name identifying features (e.g., transcripts).

    Raises
    ------
    KeyError
        If either `gene_col` or `feature_col` is missing.
    ValueError
        If:
          - any sample column is non-numeric,
          - there are missing gene or feature IDs,
          - expression values contain negatives,
          - all expression values sum to zero.
    """

    # Check that gene and feature columns exist
    required_cols = [gene_col, feature_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # Identify all putative sample columns
    sample_cols = [c for c in df.columns if c not in required_cols]
    if not sample_cols:
        raise ValueError("No sample columns found in input.")

    # All sample columns must be numeric
    non_numeric = [c for c in sample_cols
                   if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        raise ValueError(
            f"Found non-numeric columns (likely metadata) in input: {non_numeric}"
        )

    # Check for missing identifiers
    if df[gene_col].isna().any():
        raise ValueError(f"Missing values in `{gene_col}` column.")
    if df[feature_col].isna().any():
        raise ValueError(f"Missing values in `{feature_col}` column.")

    # Check that expression values are non-negative and not all zero
    expr_df = df[sample_cols]
    if (expr_df < 0).any().any():
        raise ValueError("Expression values contain negative numbers.")
    if expr_df.sum().sum() <= 0:
        raise ValueError("Total expression counts sum to zero or less.")

    return True

def rename_sample_col( df: pd.DataFrame,
                       gene_col: str = GENE_COL,
                       feature_col: str = TRANSCRIPT_COL,
                       expression_type: str = EXP_COL):
    """
    """
    sample_cols = [c for c in df.columns if c not in [gene_col,
                                                      feature_col]]
    if len(sample_cols)>1:
        raise ValueError('DF has more than one sample')

    sample = sample_cols[0]
    df.rename({sample:expression_type}, axis=1, inplace=True)
    return df

def collapse_counts_by_feature(df,
                               feature_col=TRANSCRIPT_COL,
                               expression_type=EXP_COL,
                               gene_col=GENE_COL,
                               sample_col=None):
    """
    Collapse counts by a feature (e.g. ORF, transcript) instead of transcript.

    Parameters
    ----------
    df : pd.DataFrame
        Input table with counts at transcript level.
    feature_col : str
        Alternative feature column to collapse to (e.g. 'orf_id').
    expression_type : str
        Name of expression col to collapse
    gene_col : str
        Column identifying genes.
    sample_col : str, optional
        Sample column. If None, assumes single-sample bulk.

    Returns
    -------
    pd.DataFrame
        Collapsed df with summed counts per feature.
    """
    group_cols = list(dict.fromkeys([gene_col, feature_col]))

    if sample_col is not None:
        group_cols.append(sample_col)

    # sum counts for all transcripts mapping to the same feature
    out = (df.groupby(group_cols, as_index=False)[expression_type]
          .sum())

    return out

def compute_tpm(df):
    """
    Calculate TPM values from counts for a single-sample (bulk) dataframe
    and add as a new column to the dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing 'counts'.

    Returns
    -------
    pd.DataFrame
        DataFrame with a new column 'tpm' added.
    """
    df['tpm'] = df['counts'] / df['counts'].sum() * 1e6

    return df

def compute_pi(df, gene_col=GENE_COL):
    """
    Generate pi values (isoform ratios) from input expression column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing at least <gene_col> and 'tpm'.
    gene_col : str
        Column representing the gene (default: gene_id)

    Returns
    -------
    pd.DataFrame
        DataFrame with an additional 'pi' column.
    """
    df['gene_tpm'] = df.groupby(gene_col)['tpm'].transform('sum')
    df['pi'] = df['tpm'] / df['gene_tpm']
    df['pi'] = df['pi'].fillna(0)
    df = df.drop(columns='gene_tpm')
    return df

def compute_n_detected_features(df, gene_col=GENE_COL, feature_col=TRANSCRIPT_COL):
    """
    Compute gene potential based on the number of unique expressed features per gene.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing gene and feature columns.
    gene_col : str
        Column representing the gene (default: 'gene_id')
    feature_col : str
        Column representing the feature (default: 'transcript_id')

    Returns
    -------
    pd.DataFrame
        DataFrame with an additional column:
        - 'n_detected_features': number of unique features per gene
    """
    # count unique expressed features per gene
    temp = df.loc[df.tpm>0].copy(deep=True)
    n_detected_features_df = (
        temp[[gene_col, feature_col]]
        .groupby(gene_col)
        .nunique()
        .reset_index()
        .rename({feature_col: 'n_detected_features'}, axis=1)
    )

    # merge back to original df
    df = df.merge(n_detected_features_df, how='left', on=gene_col)

    df['n_detected_features'] = df['n_detected_features'].fillna(0)

    return df

def compute_entropy(df, gene_col=GENE_COL):
    """
    Compute Shannon entropy per gene based on isoform proportions.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing gene and 'pi' columns.
    gene_col : str
        Column representing the gene (default: 'gene_id')

    Returns
    -------
    pd.DataFrame
        DataFrame with an additional column:
        - 'entropy': Shannon entropy per gene
    """
    # compute plogp for each feature; avoid log2(0)
    df['plogp'] = df['pi'] * np.log2(df['pi'].replace(0, np.nan))

    # sum plogp per gene
    entropy_df = (
        df[[gene_col, 'plogp']]
        .groupby(gene_col)
        .sum()
        .reset_index()
        .rename({'plogp': 'entropy'}, axis=1)
    )

    # multiply by -1 for positive entropy
    entropy_df['entropy'] = -1 * entropy_df['entropy']

    # merge back to original df
    df = df.merge(entropy_df, how='left', on=gene_col)

    # drop intermediate column
    df = df.drop(columns='plogp')

    return df

def compute_perplexity(df):
    """
    Compute perplexity per gene based on Shannon entropy.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'gene_id' and 'entropy' columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with an additional column:
        - 'perplexity': effective number of isoforms per gene
    """
    # compute perplexity as 2^entropy
    df['perplexity'] = 2 ** df['entropy']
    return df

def call_effective(df, gene_col=GENE_COL, feature_col=TRANSCRIPT_COL):
    """
    Identify effective features per gene based on gene perplexity.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing gene, feature, 'pi', and 'perplexity' columns.
    gene_col : str
        Column representing the gene (default: 'gene_id')
    feature_col : str
        Column representing the feature (default: 'transcript_id')

    Returns
    -------
    pd.DataFrame
        DataFrame with additional columns:
        - 'n_effective': perplexity rounded to nearest integer
        - 'feature_rank': rank of each feature within its gene (by pi)
        - 'effective': boolean indicating if feature is effective
    """
    # round perplexity to nearest integer
    df['n_effective'] = df['perplexity'].round(0)

    # rank features within each gene by pi
    df['feature_rank'] = (
        df.groupby(gene_col)['pi']
        .rank(method='first', ascending=False)
        .astype(int)
    )

    # mark features as effective if rank <= rounded perplexity
    df['effective'] = df['feature_rank'] <= df['n_effective']

    return df

def compute_expression_breadth(df, sample_col, feature_col=TRANSCRIPT_COL):
    """
    Compute percentage of samples in which each feature is effective.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing feature, sample, and 'effective' columns.
    sample_col : str
        Column representing sample IDs.
    feature_col : str
        Column representing the feature (default: 'transcript_id')

    Returns
    -------
    pd.DataFrame
        DataFrame with additional columns:
        - 'n_samples_effective': number of samples where the feature is effective
        - 'expression_breadth': percentage of samples where the feature is effective
    """
    n_samples = df[sample_col].nunique()

    temp = (
        df.loc[df['effective'], [feature_col, sample_col]]
        .groupby(feature_col)
        .nunique()
        .reset_index()
        .rename({sample_col: 'n_samples_effective'}, axis=1)
    )

    df = df.merge(temp, how='left', on=feature_col)
    df['expression_breadth'] = df['n_samples_effective'].fillna(0) / n_samples * 100

    return df

def compute_expression_var(df, sample_col, feature_col=TRANSCRIPT_COL):
    """
    Compute number of samples expressing each feature and pi standard deviation.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing feature, sample, and 'pi' columns.
    sample_col : str
        Column representing sample IDs.
    feature_col : str
        Column representing the feature (default: 'transcript_id')

    Returns
    -------
    pd.DataFrame
        DataFrame with additional columns:
        - 'n_exp_samples': number of samples where feature is expressed
        - 'expression_var': standard deviation of pi across samples
    """
    # number of samples where feature is expressed
    df['n_exp_samples'] = df.groupby(feature_col)[sample_col].transform('nunique')

    # standard deviation of pi across samples
    df['expression_var'] = df.groupby(feature_col)['pi'].transform(lambda x: x.std(ddof=1, skipna=True))

    return df

def compute_max_expression(df, sample_col, feature_col=TRANSCRIPT_COL):
    """
    Compute maximum expression across samples for each feature.
    Only considers samples where feature is expressed!

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing feature, sample, and counts columns.
    sample_col : str
        Column representing sample IDs.
    feature_col : str
        Column representing the feature (default: 'transcript_id')

    Returns
    -------
    pd.DataFrame
        DataFrame with additional column:
        - 'avg_<feature_col>_tpm': mean tpm across samples for the feature
    """
    temp = df[[feature_col, sample_col, 'tpm']]

    # filter out unexpressed features so mean denominator is correct
    temp = temp.loc[temp['tpm'] > 0]

    # sum counts per feature per sample (handles cases with multiple rows per feature)
    temp = (
        temp.groupby([feature_col, sample_col])
        .sum()
        .reset_index()
        .rename({'tpm': f'{feature_col}_tpm'}, axis=1)
    )
    temp.drop(sample_col, axis=1, inplace=True)

    # then take mean across samples
    temp = (
        temp.groupby(feature_col)
        .max()
        .reset_index()
        .rename({f'{feature_col}_tpm': f'max_{feature_col}_tpm'}, axis=1)
    )

    # merge back to original df
    df = df.merge(temp, how='left', on=[feature_col])

    return df

def compute_avg_expression(df, sample_col, feature_col=TRANSCRIPT_COL):
    """
    Compute average expression across samples for each feature.
    Only considers samples where feature is expressed!

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing feature, sample, and counts columns.
    sample_col : str
        Column representing sample IDs.
    feature_col : str
        Column representing the feature (default: 'transcript_id')

    Returns
    -------
    pd.DataFrame
        DataFrame with additional column:
        - 'avg_<feature_col>_tpm': mean tpm across samples for the feature
    """
    temp = df[[feature_col, sample_col, 'tpm']]

    # filter out unexpressed features so mean denominator is correct
    temp = temp.loc[temp['tpm'] > 0]

    # sum counts per feature per sample (handles cases with multiple rows per feature)
    temp = (
        temp.groupby([feature_col, sample_col])
        .sum()
        .reset_index()
        .rename({'tpm': f'{feature_col}_tpm'}, axis=1)
    )
    temp.drop(sample_col, axis=1, inplace=True)

    # then take mean across samples
    temp = (
        temp.groupby(feature_col)
        .mean()
        .reset_index()
        .rename({f'{feature_col}_tpm': f'avg_{feature_col}_tpm'}, axis=1)
    )

    # merge back to original df
    df = df.merge(temp, how='left', on=[feature_col])

    return df

def compute_global_isoform_metrics(df,
                                   gene_col=GENE_COL,
                                   feature_col=TRANSCRIPT_COL,
                                   expression_type='counts'):
    """
    Compute isoform or other feature diversity metrics for a single-sample (bulk) dataframe.
    Either provide counts or TPMs; if counts, will automatically convert to TPM.
    Optionally, collapse counts to a different feature or compute TPM.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing counts at transcript level.
    gene_col : str
        Column representing the gene (default: 'gene_id')
    feature_col : str
        Column representing the feature (default: 'transcript_id')
    expression_type : str
        Type of expression values in table {'counts' | 'tpm'}

    Returns
    -------
    pd.DataFrame
        DataFrame with computed metrics.
    """
    # validate input
    validate_counts_input(df,
                          gene_col=gene_col,
                          feature_col=feature_col)

    # convert from wide to long
    df = rename_sample_col(df,
                      gene_col=gene_col,
                      feature_col=feature_col,
                      expression_type=expression_type)

    # collapse counts if feature_col is not unique
    df = collapse_counts_by_feature(df,
                                    feature_col=feature_col,
                                    expression_type=expression_type,
                                    gene_col=gene_col,
                                    sample_col=None)

    # compute TPM if required
    df = df if expression_type == 'tpm' else compute_tpm(df)

    # compute isoform ratios
    df = compute_pi(df, gene_col=gene_col)

    # compute gene-level metrics
    df = compute_n_detected_features(df, gene_col=gene_col, feature_col=feature_col)
    df = compute_entropy(df, gene_col=gene_col)
    df = compute_perplexity(df)

    # mark effective features
    df = call_effective(df, gene_col=gene_col, feature_col=feature_col)

    return df

def compute_multi_sample_isoform_metrics(
    df: pd.DataFrame,
    gene_col: str = GENE_COL,
    feature_col: str = TRANSCRIPT_COL,
    expression_type: str = EXP_COL):
    """
    Compute isoform metrics across multiple samples as well as global
    metrics.

    Parameters
    ----------
    df : pd.DataFrame
        Input table with counts or TPMs for all samples.
    gene_col : str
        Column identifying genes.
    feature_col : str
        Column identifying isoforms or other features (e.g. ORFs).
    expression_type : {'counts', 'tpm'}
        Type of expression values

    Returns
    -------
    big_df : pd.DataFrame
        DataFrame with:
          • per-sample metrics (gene potential, entropy, etc.)
          • cross-sample metrics (breadth, variance, average expression)
    global_df : pd.DataFrame
        DataFrame with:
          • global entropy, detected features, and perplexity per gene,
            based on summing tpms across all samples
    """

    # validate input
    validate_counts_input(df,
                          gene_col=gene_col,
                          feature_col=feature_col)

    # prepare kwargs for downstream functions
    col_kwargs = dict(
        gene_col=gene_col,
        feature_col=feature_col,
        expression_type=expression_type,
    )

    # loop over samples and compute single-sample metrics
    samples = [c for c in df.columns if c not in [gene_col, feature_col]]
    dfs = []

    for i, s in enumerate(samples):
        s_df = df[[gene_col, feature_col, s]].copy()
        df.drop(s, axis=1, inplace=True) # remove from original df to save memory
        # if i % 25 == 0: print(i)
        s_df = compute_global_isoform_metrics(s_df,
                                              gene_col=gene_col,
                                              feature_col=feature_col,
                                              expression_type=expression_type)
        s_df['sample'] = s  # re-attach sample ID
        dfs.append(s_df)

    big_df = pd.concat(dfs, axis=0, ignore_index=True)

    # compute cross-sample metrics on the combined table
    big_df = compute_expression_breadth(big_df,
                                        sample_col='sample',
                                        feature_col=feature_col)
    big_df = compute_expression_var(big_df,
                                    sample_col='sample',
                                    feature_col=feature_col)

    for fc in [feature_col, gene_col]:
        big_df = compute_avg_expression(big_df,
                                        sample_col='sample',
                                        feature_col=fc)
        big_df = compute_max_expression(big_df,
                                        sample_col='sample',
                                        feature_col=fc)

    # sum up tpms across samples
    global_df = (big_df[[gene_col, feature_col, 'tpm']]
          .groupby([gene_col, feature_col])
          .sum()
          .reset_index()
         )

    global_df.rename({'tpm': 'global'}, axis=1, inplace=True)



    global_df = compute_global_isoform_metrics(global_df,
                                               gene_col=gene_col,
                                               feature_col=feature_col,
                                               expression_type='tpm')

    return big_df, global_df

def flatten_list(l):
    """
    Flatten a list into 1 dimension.

    Parameters
    ----------
    l : list
    """
    return [j for i in l for j in i]
