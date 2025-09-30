import scipy.stats
import pandas as pd
import numpy as np


import pandas as pd
import numpy as np

def validate_counts_input(
    df: pd.DataFrame,
    sample_col: str = None):
    """
    Validate input dataframe for computing isoform ratios (pi).

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing at least 'gene_id', 'transcript_id', and 'counts'.
    sample_col : str, optional
        If provided, checks that this column exists and has no missing values.

    Raises
    ------
    KeyError
        If required columns are missing.
    ValueError
        If counts contain invalid (negative or disallowed zero) values,
        or if transcript_id is missing or duplicated inappropriately.
    """

    required_cols = ['gene_id', 'transcript_id', 'counts']
    if sample_col is not None:
        required_cols.append(sample_col)

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # check for null IDs
    if df['gene_id'].isna().any():
        raise ValueError("Missing gene_id values.")
    if df['transcript_id'].isna().any():
        raise ValueError("Missing transcript_id values.")

    # counts validation
    if (df['counts'] < 0).any():
        raise ValueError("Counts contain negative values.")
        
    # check for duplicated identifiers
    if sample_col is not None:
        if df.duplicated(subset=[sample_col, 'transcript_id']).any():
            raise ValueError('found duplicate (sample, transcript_id) rows')
    else:
        if df['transcript_id'].duplicated().any():
            raise ValueError('found duplicate transcript_id rows')

    # sample column checks
    if sample_col is not None:
        if df[sample_col].isna().any():
            raise ValueError(f"Sample column '{sample_col}' contains missing values.")
        if df[sample_col].nunique() == 0:
            raise ValueError("No samples found in the dataframe.")

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

def compute_pi(df):
    """
    Generate pi values (isoform ratios) from raw counts.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing at least 'gene_id', 'transcript_id', and 'counts'.
        Assumes missing transcript-sample combinations have been filled with counts=0.

    Returns
    -------
    pd.DataFrame
        DataFrame with additional columns:
        - 'gene_counts': total counts per gene
        - 'pi': isoform proportion within each gene
    """
    # compute total counts per gene
    df['gene_counts'] = df.groupby('gene_id')['counts'].transform('sum')
    
    # compute pi
    df['pi'] = df['counts'] / df['gene_counts']
    
    return df

def compute_gene_potential(df):
    """
    Compute gene potential based on the number of expressed isoforms per gene.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing at least 'gene_id' and 'transcript_id'.

    Returns
    -------
    pd.DataFrame
        DataFrame with an additional column:
        - 'gene_potential': number of unique transcripts per gene
    """
    # compute number of isoforms per gene
    gene_isoform_counts = (
        df[['gene_id', 'transcript_id']]
        .groupby('gene_id')
        .nunique()
        .reset_index()
        .rename({'transcript_id': 'gene_potential'}, axis=1)
    )
    
    # merge back to original df
    df = df.merge(gene_isoform_counts, how='left', on='gene_id')
    
    return df

def compute_entropy(df):
    """
    Compute Shannon entropy per gene based on isoform proportions (pi).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'gene_id' and 'pi' columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with an additional column:
        - 'entropy': Shannon entropy per gene
    """
    # compute plogp for each isoform
    df['plogp'] = df['pi'] * np.log2(df['pi'].replace(0, np.nan))
    
    # sum over isoforms per gene
    entropy_df = (
        df[['gene_id', 'plogp']]
        .groupby('gene_id')
        .sum()
        .reset_index()
        .rename({'plogp': 'entropy'}, axis=1)
    )
    
    # multiply by -1 to get positive entropy
    entropy_df['entropy'] = -1 * entropy_df['entropy']
    
    # merge back to original dataframe
    df = df.merge(entropy_df, how='left', on='gene_id')
    
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

def call_effective_isoforms(df):
    """
    Identify effective isoforms per gene based on gene perplexity.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'gene_id', 'pi', and 'perplexity' columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with additional columns:
        - 'round_perplexity': perplexity rounded to nearest integer
        - 'isoform_rank': rank of each isoform within its gene (by pi)
        - 'effective_isoform': boolean indicating if isoform is effective
    """
    # round perplexity to nearest integer
    df['round_perplexity'] = df['perplexity'].round(0)
    
    # rank isoforms within each gene by pi
    df['isoform_rank'] = (
        df.groupby('gene_id')['pi']
        .rank(method='first', ascending=False)
        .astype(int)
    )
    
    # mark isoforms as effective if rank <= rounded perplexity
    df['effective_isoform'] = df['isoform_rank'] <= df['round_perplexity']
    
    return df

def compute_isoform_metrics(df):
    """
    Compute isoform diversity metrics for a single-sample (bulk) dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing 'gene_id', 'transcript_id', and 'counts'.

    Returns
    -------
    pd.DataFrame
        DataFrame with computed metrics:
        - 'gene_counts', 'pi', 'gene_potential', 'entropy', 'perplexity',
          'round_perplexity', 'isoform_rank', 'effective_isoform'
    """
    # validate input
    validate_counts_input(df)

    # compute isoform ratios
    df = compute_pi(df)

    # compute gene-level metrics
    df = compute_gene_potential(df)
    df = compute_entropy(df)
    df = compute_perplexity(df)

    # mark effective isoforms
    df = call_effective_isoforms(df)

    return df

def compute_expression_breadth(df, sample_col):
    """
    Compute percentage of samples in which each isoform is effective.
    """
    n_samples = df[sample_col].nunique()
    
    temp = (
        df.loc[df.effective_isoform, ['transcript_id', sample_col]]
        .groupby('transcript_id')
        .nunique()
        .reset_index()
        .rename({sample_col: 'n_samples_effective'}, axis=1)
    )
    
    df = df.merge(temp, how='left', on='transcript_id')
    df['perc_effective_isoforms'] = df['n_samples_effective'].fillna(0) / n_samples * 100
    
    return df

def compute_expression_var(df, sample_col):
    """
    Compute number of samples expressing each isoform and pi standard deviation.
    """
    df['n_exp_samples'] = df.groupby('transcript_id')[sample_col].transform('nunique')
    df['pi_std'] = df.groupby('transcript_id')['pi'].transform(lambda x: x.std(ddof=1, skipna=True))
    
    return df

def compute_multi_sample_isoform_metrics(df, sample_col):
    """
    Compute isoform metrics across multiple samples.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'gene_id', 'transcript_id', 'counts', and sample_col.
    sample_col : str
        Name of the column identifying samples.

    Returns
    -------
    pd.DataFrame
        DataFrame with single-sample metrics plus cross-sample metrics:
        - 'perc_effective_isoforms', 'n_exp_samples', 'pi_std'
    """
    # validate input for multi-sample dataset
    validate_counts_input(df, sample_col=sample_col)
    
    samples = df[sample_col].unique().tolist()
    big_df = pd.DataFrame()
    
    for s in samples:
        s_df = df.loc[df[sample_col] == s].copy(deep=True)
        # delegate single-sample calculations
        s_df = compute_isoform_metrics(s_df)
        big_df = pd.concat([big_df, s_df], axis=0)
    
    # compute cross-sample metrics
    big_df = compute_expression_breadth(big_df, sample_col=sample_col)
    big_df = compute_expression_var(big_df, sample_col=sample_col)
    
    return big_df

def flatten_list(l):
    """
    Flatten a list into 1 dimension.

    Parameters
    ----------
    l : list
    """
    return [j for i in l for j in i]
