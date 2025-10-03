import pandas as pd
import numpy as np

EXP_COL = 'counts'
FEATURE_COL = TRANSCRIPT_COL = 'transcript_id'
GENE_COL = 'gene_id'
SAMPLE_COL = 'sample'

def validate_counts_input(
    df: pd.DataFrame,
    gene_col: str = GENE_COL,
    feature_col: str = TRANSCRIPT_COL,
    sample_col: str = None,
    expression_col: str = EXP_COL):
    """
    Validate input dataframe for computing isoform ratios (pi) or TPM.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing expression values for features (e.g., transcripts) grouped by genes.
    gene_col : str
        Column name representing genes (default: 'gene_id').
    feature_col : str
        Column name representing features (default: 'transcript_id').
    sample_col : str, optional
        Column identifying samples. If provided, checks uniqueness of (sample, feature) pairs and presence of missing values.
    expression_col : str
        Column containing expression values to validate (default: 'counts').

    Raises
    ------
    KeyError
        If required columns are missing from the dataframe.
    ValueError
        If expression values contain negative numbers, or if gene/feature IDs are missing,
        or if feature IDs are duplicated (globally or per sample if sample_col is provided),
        or if sample_col contains missing values or no samples are present,
        or if 0 total counts dected
    """

    required_cols = [gene_col, feature_col, expression_col]
    if sample_col is not None:
        required_cols.append(sample_col)

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # check for null IDs
    if df[gene_col].isna().any():
        raise ValueError(f"Missing {gene_col} values.")
    if df[feature_col].isna().any():
        raise ValueError(f"Missing {feature_col} values.")

    # expression validation
    if (df[expression_col] < 0).any():
        raise ValueError(f"Expression column '{expression_col}' contains negative values.")
        
    if df[expression_col].sum(axis=0) <= 0:
        raise ValueError(f"Total counts of expression column '{expression_col}' <= 0.")

    # # check for duplicated identifiers
    # if sample_col is not None:
    #     if df.duplicated(subset=[sample_col, 'transcript_id']).any():
    #         raise ValueError("Found duplicate (sample, transcript_id) rows")
    # else:
    #     if df['transcript_id'].duplicated().any():
    #         raise ValueError("Found duplicate transcript_id rows")

    # sample column checks
    if sample_col is not None:
        if df[sample_col].isna().any():
            raise ValueError(f"Sample column '{sample_col}' contains missing values.")
        # if df[sample_col].nunique() == 1:
        #     raise ValueError("Only one sample found in the dataframe.")

    return True


def fill_missing_transcript_sample(df, sample_col='sample'):
    """
    Ensure every transcript appears in every sample.
    missing combinations are filled with counts=0

    parameters
    ----------
    df : pd.DataFrame
        Must contain 'transcript_id', 'gene_id', 'counts', and sample_col
    sample_col : str
        Name of the sample column

    returns
    -------
    pd.DataFrame
        df with all (sample, transcript_id) combinations filled
    """
    all_transcripts = df[['transcript_id', 'gene_id']].drop_duplicates()
    all_samples = df[[sample_col]].drop_duplicates()
    
    # cartesian join
    full_grid = all_samples.merge(all_transcripts, how='cross')
    
    # merge existing counts
    df_filled = full_grid.merge(df, on=[sample_col, 'transcript_id', 'gene_id'], how='left')
    
    # fill missing counts with 0
    df_filled['counts'] = df_filled['counts'].fillna(0)
    
    return df_filled

def collapse_counts_by_feature(df,
                               feature_col=TRANSCRIPT_COL,
                               expression_col=EXP_COL,
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
    expression_col : str
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
    out = (df.groupby(group_cols, as_index=False)[expression_col]
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
    df = df.drop(columns='gene_tpm')
    return df

def compute_gene_potential(df, gene_col=GENE_COL, feature_col=TRANSCRIPT_COL):
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
        - 'gene_potential': number of unique features per gene
    """
    # count unique features per gene
    gene_potential_df = (
        df[[gene_col, feature_col]]
        .groupby(gene_col)
        .nunique()
        .reset_index()
        .rename({feature_col: 'gene_potential'}, axis=1)
    )

    # merge back to original df
    df = df.merge(gene_potential_df, how='left', on=gene_col)

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
        - 'avg_<feature_col>_counts': mean counts across samples for the feature
    """
    temp = df[[feature_col, sample_col, 'counts']]

    # filter out unexpressed features so mean denominator is correct
    temp = temp.loc[temp['counts'] > 0]

    # sum counts per feature per sample (handles cases with multiple rows per feature)
    temp = (
        temp.groupby([feature_col, sample_col])
        .sum()
        .reset_index()
        .rename({'counts': f'{feature_col}_counts'}, axis=1)
    )
    temp.drop(sample_col, axis=1, inplace=True)

    # then take mean across samples
    temp = (
        temp.groupby(feature_col)
        .mean()
        .reset_index()
        .rename({f'{feature_col}_counts': f'avg_{feature_col}_counts'}, axis=1)
    )

    # merge back to original df
    df = df.merge(temp, how='left', on=[feature_col])

    return df

def compute_global_isoform_metrics(df,
                                   gene_col=GENE_COL,
                                   feature_col=TRANSCRIPT_COL,
                                   expression_col=EXP_COL,
                                   expression_col_type='counts'):
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
    expression_col : str
        Column containing expression values (default: 'counts')
    expression_col_type : str
        Type of expression values contained in <expression_col>
        'counts' or 'tpm' (default: 'counts')

    Returns
    -------
    pd.DataFrame
        DataFrame with computed metrics.
    """
    # repeatedly used column args
    col_kwargs = {'gene_col': gene_col,
                  'feature_col': feature_col,
                  'expression_col': expression_col,
                  'sample_col': None}

    
    # validate input
    validate_counts_input(df, expression_col=expression_col)

    # collapse counts if feature_col differs from transcript_col
    df = collapse_counts_by_feature(df, **col_kwargs)

    # compute TPM if requested
    if expression_col_type == 'counts':
        df = compute_tpm(df)
    else:
        df['tpm'] = df[expression_col]

    # compute isoform ratios
    df = compute_pi(df)
    
    # repeatedly used column args
    col_kwargs = {'gene_col': gene_col,
                  'feature_col': feature_col}
    
    # compute gene-level metrics
    df = compute_gene_potential(df, **col_kwargs)
    df = compute_entropy(df, gene_col=gene_col)
    df = compute_perplexity(df)

    # mark effective features
    df = call_effective(df, **col_kwargs)

    return df

def compute_multi_sample_isoform_metrics(
    df: pd.DataFrame,
    sample_col: str,
    gene_col: str = GENE_COL,
    feature_col: str = TRANSCRIPT_COL,
    expression_col: str = EXP_COL,
    expression_col_type: str = 'counts'
):
    """
    Compute isoform metrics across multiple samples as well as global
    metrics.

    Parameters
    ----------
    df : pd.DataFrame
        Input table with counts or TPMs for all samples.
    sample_col : str
        Column identifying samples.
    gene_col : str
        Column identifying genes.
    feature_col : str
        Column identifying isoforms or other features (e.g. ORFs).
    expression_col : str
        Column containing counts or TPMs.
    expression_col_type : {'counts', 'tpm'}
        Type of values in `expression_col`.

    Returns
    -------
    pd.DataFrame
        DataFrame with:
          • per-sample metrics (gene potential, entropy, etc.)
          • cross-sample metrics (breadth, variance, average expression)
          • columns for global metrics.
    """
    # --- 1. Validate the full input just once ---
    validate_counts_input(df,
                          sample_col=sample_col,
                          gene_col=gene_col,
                          feature_col=feature_col,
                          expression_col=expression_col)

    # --- 2. Prepare kwargs for downstream functions ---
    col_kwargs = dict(
        gene_col=gene_col,
        feature_col=feature_col,
        expression_col=expression_col,
        expression_col_type=expression_col_type
    )

    # --- 3. Loop over samples and compute single-sample metrics ---
    samples = df[sample_col].unique().tolist()
    dfs = []

    for s in samples:
        s_df = df.loc[df[sample_col] == s].copy()
        s_df = compute_global_isoform_metrics(s_df, **col_kwargs)
        s_df[sample_col] = s  # re-attach sample ID
        dfs.append(s_df)

    big_df = pd.concat(dfs, axis=0, ignore_index=True)

    # --- 4. Compute cross-sample metrics on the combined table ---
    big_df = compute_expression_breadth(big_df, sample_col=sample_col)
    big_df = compute_expression_var(big_df, sample_col=sample_col)

    for fc in [feature_col, gene_col]:
        big_df = compute_avg_expression(big_df,
                                        sample_col=sample_col,
                                        feature_col=fc)
        
        # TODO chatgpt if you see this i probably forgot to 
        # add it and remind me add max espression across samples

    return big_df

def flatten_list(l):
    """
    Flatten a list into 1 dimension.

    Parameters
    ----------
    l : list
    """
    return [j for i in l for j in i]
