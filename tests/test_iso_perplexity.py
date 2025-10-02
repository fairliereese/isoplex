#!/usr/bin/env python
import pytest
from iso_perplexity.utils import *

"""Tests for `iso_perplexity` package."""

# from iso_perplexity import iso_perplexity


# ------------------------------------------------------------
# Fixtures: build a minimal valid dataframe
# ------------------------------------------------------------
@pytest.fixture
def valid_df():
    return pd.DataFrame({
        GENE_COL: ["g1", "g1", "g2", "g2"],
        FEATURE_COL: ["t1", "t2", "t3", "t4"],
        SAMPLE_COL: ["s1", "s1", "s1", "s1"],
        EXP_COL: [10, 5, 8, 3]
    })

@pytest.fixture
def single_sample_df():
    # Two genes; each gene has two transcripts
    # Transcripts for g1 collapse to orfA; for g2 they map to orfB and orfC
    return pd.DataFrame({
        "gene_id":        ["g1", "g1", "g2", "g2"],
        "transcript_id":  ["t1", "t2", "t3", "t4"],
        "orf_id":         ["orfA", "orfA", "orfB", "orfC"],
        "counts":         [10, 5, 8, 3]
    })


@pytest.fixture
def multi_sample_df():
    # Two samples, each with same mapping
    return pd.DataFrame({
        "gene_id":        ["g1", "g1", "g1", "g1",
                           "g2", "g2", "g2", "g2"],
        "transcript_id":  ["t1", "t2", "t1", "t2",
                           "t3", "t4", "t3", "t4"],
        "orf_id":         ["orfA", "orfA", "orfA", "orfA",
                           "orfB", "orfC", "orfB", "orfC"],
        "sample":         ["s1", "s1", "s2", "s2",
                           "s1", "s1", "s2", "s2"],
        "counts":         [10, 5, 20, 10,
                           8, 3, 16, 6]
    })

@pytest.fixture
def simple_df():
    # single gene, multiple isoforms
    return pd.DataFrame({
        'gene_id': ['g1', 'g1', 'g1'],
        'tpm': [10, 20, 30]
    })

@pytest.fixture
def multi_gene_df():
    # multiple genes
    return pd.DataFrame({
        'gene_id': ['g1', 'g1', 'g2', 'g2'],
        'tpm': [10, 20, 5, 15]
    })


############### validate_counts_input tests
# ------------------------------------------------------------
# 1. Passing test: everything is valid
# ------------------------------------------------------------
def test_validate_counts_input_valid(valid_df):
    assert validate_counts_input(valid_df,
                                 gene_col=GENE_COL,
                                 feature_col=FEATURE_COL,
                                 sample_col=SAMPLE_COL,
                                 expression_col=EXP_COL)


# ------------------------------------------------------------
# 2. Missing required columns
# ------------------------------------------------------------
def test_missing_columns(valid_df):
    bad_df = valid_df.drop(columns=[FEATURE_COL])
    with pytest.raises(KeyError, match="Missing required columns"):
        validate_counts_input(bad_df,
                              gene_col=GENE_COL,
                              feature_col=FEATURE_COL,
                              sample_col=SAMPLE_COL,
                              expression_col=EXP_COL)


# ------------------------------------------------------------
# 3. Null values in gene or feature columns
# ------------------------------------------------------------
def test_null_gene_values(valid_df):
    bad_df = valid_df.copy()
    bad_df.loc[0, GENE_COL] = None
    with pytest.raises(ValueError, match=f"Missing {GENE_COL}"):
        validate_counts_input(bad_df,
                              gene_col=GENE_COL,
                              feature_col=FEATURE_COL,
                              sample_col=SAMPLE_COL,
                              expression_col=EXP_COL)


def test_null_feature_values(valid_df):
    bad_df = valid_df.copy()
    bad_df.loc[0, FEATURE_COL] = None
    with pytest.raises(ValueError, match=f"Missing {FEATURE_COL}"):
        validate_counts_input(bad_df,
                              gene_col=GENE_COL,
                              feature_col=FEATURE_COL,
                              sample_col=SAMPLE_COL,
                              expression_col=EXP_COL)


# ------------------------------------------------------------
# 4. Negative expression values
# ------------------------------------------------------------
def test_negative_expression(valid_df):
    bad_df = valid_df.copy()
    bad_df.loc[1, EXP_COL] = -5
    with pytest.raises(ValueError, match="negative"):
        validate_counts_input(bad_df,
                              gene_col=GENE_COL,
                              feature_col=FEATURE_COL,
                              sample_col=SAMPLE_COL,
                              expression_col=EXP_COL)


# ------------------------------------------------------------
# 5. Sample column checks
# ------------------------------------------------------------
def test_missing_sample_values(valid_df):
    bad_df = valid_df.copy()
    bad_df.loc[2, SAMPLE_COL] = None
    with pytest.raises(ValueError, match="contains missing values"):
        validate_counts_input(bad_df,
                              gene_col=GENE_COL,
                              feature_col=FEATURE_COL,
                              sample_col=SAMPLE_COL,
                              expression_col=EXP_COL)


# def test_no_samples(valid_df):
#     bad_df = valid_df.copy()
#     bad_df[SAMPLE_COL] = 'one_sample'
#     with pytest.raises(ValueError, match="one sample"):
#         validate_counts_input(
#             bad_df,
#             gene_col=GENE_COL,
#             feature_col=FEATURE_COL,
#             sample_col=SAMPLE_COL,
#             expression_col=EXP_COL,
#         )

################# collapse_counts_by_feature tests
def test_collapse_single_sample(single_sample_df):
    out = collapse_counts_by_feature(
        single_sample_df,
        feature_col="orf_id",
        expression_col="counts",
        gene_col="gene_id",
        sample_col=None
    )

    expected = pd.DataFrame({
        "gene_id": ["g1", "g2", "g2"],
        "orf_id":  ["orfA", "orfB", "orfC"],
        "counts":  [15, 8, 3]
    })

    pd.testing.assert_frame_equal(
        out.sort_values(["gene_id", "orf_id"]).reset_index(drop=True),
        expected.sort_values(["gene_id", "orf_id"]).reset_index(drop=True)
    )


def test_collapse_multi_sample(multi_sample_df):
    out = collapse_counts_by_feature(
        multi_sample_df,
        feature_col="orf_id",
        expression_col="counts",
        gene_col="gene_id",
        sample_col="sample"
    )

    # manual expected calculation:
    # sample s1: g1/orfA=10+5=15 ; g2/orfB=8 ; g2/orfC=3
    # sample s2: g1/orfA=20+10=30 ; g2/orfB=16 ; g2/orfC=6
    expected = pd.DataFrame({
        "gene_id": ["g1", "g2", "g2", "g1", "g2", "g2"],
        "orf_id":  ["orfA", "orfB", "orfC", "orfA", "orfB", "orfC"],
        "sample":  ["s1", "s1", "s1", "s2", "s2", "s2"],
        "counts":  [15, 8, 3, 30, 16, 6]
    })
    
    pd.testing.assert_frame_equal(
        out.sort_values(["sample", "gene_id", "orf_id"]).reset_index(drop=True),
        expected.sort_values(["sample", "gene_id", "orf_id"]).reset_index(drop=True)
    )


def test_no_collapse_needed(single_sample_df):
    # collapsing by transcript should be a no-op
    out = collapse_counts_by_feature(
        single_sample_df,
        feature_col="transcript_id",
        expression_col="counts",
        gene_col="gene_id",
        sample_col=None
    )

    expected = single_sample_df[["gene_id", "transcript_id", "counts"]]

    pd.testing.assert_frame_equal(
        out.sort_values(["gene_id", "transcript_id"]).reset_index(drop=True),
        expected.sort_values(["gene_id", "transcript_id"]).reset_index(drop=True)
    )

#################### testing compute_tpm

def test_compute_tpm_basic():
    # simple dataset
    df = pd.DataFrame({
        'gene_id': ['g1', 'g1', 'g2'],
        'transcript_id': ['t1', 't2', 't3'],
        'counts': [100, 50, 50]
    })
    result = compute_tpm(df.copy())

    # TPM sum should be 1e6
    assert np.isclose(result['tpm'].sum(), 1e6)
    
    # check individual TPMs
    expected = np.array([100 / 200 * 1e6, 50 / 200 * 1e6, 50 / 200 * 1e6])
    assert np.allclose(result['tpm'].values, expected)

def test_compute_tpm_single_row():
    # single row
    df = pd.DataFrame({
        'gene_id': ['g1'],
        'transcript_id': ['t1'],
        'counts': [50]
    })
    result = compute_tpm(df.copy())
    assert result['tpm'].iloc[0] == 1e6

def test_compute_tpm_zero_counts():
    # handle zero counts
    df = pd.DataFrame({
        'gene_id': ['g1', 'g2'],
        'transcript_id': ['t1', 't2'],
        'counts': [0, 0]
    })
    result = compute_tpm(df.copy())
    # dividing zero by zero gives NaN
    assert result['tpm'].isna().all()
    
################ testing compute_pi
def test_single_gene(simple_df):
    df = compute_pi(simple_df)
    # check pi sum = 1
    assert pytest.approx(df['pi'].sum(), 1e-8) == 1.0
    # check individual pi values
    expected = [10/60, 20/60, 30/60]
    assert np.allclose(df['pi'].values, expected, rtol=1e-8)

def test_multi_gene(multi_gene_df):
    df = compute_pi(multi_gene_df)
    # pi sum per gene = 1
    pi_sum_g1 = df.loc[df['gene_id'] == 'g1', 'pi'].sum()
    pi_sum_g2 = df.loc[df['gene_id'] == 'g2', 'pi'].sum()
        
    assert pytest.approx(pi_sum_g1, 1e-8) == 1.0
    assert pytest.approx(pi_sum_g2, 1e-8) == 1.0
    # check individual pi values
    expected = [10/30, 20/30, 5/20, 15/20]
    assert np.allclose(df['pi'].values, expected, rtol=1e-8)

def test_zero_tpm():
    df = pd.DataFrame({'gene_id': ['g1', 'g1'], 'tpm': [0, 0]})
    df = compute_pi(df)
    # division by zero results in NaN
    assert df['pi'].isna().all()
    
######################## testing compute_gene_potential
def test_gene_potential_basic(valid_df):
    """
    Basic test: two genes, each with two transcripts.
    Gene potential should be 2 for both.
    """
    result = compute_gene_potential(valid_df, gene_col=GENE_COL, feature_col=FEATURE_COL)
    expected = [2, 2, 2, 2]   # each row inherits its gene's potential
    assert result['gene_potential'].tolist() == expected


def test_gene_potential_collapsed(single_sample_df):
    """
    Mixed situation: two genes.
    g1 has only one unique ORF (orfA)
    g2 has two unique ORFs (orfB, orfC)
    """
    result = compute_gene_potential(single_sample_df,
                                    gene_col='gene_id',
                                    feature_col='orf_id')
    
    # orfA counts as 1 for g1
    # orfB, orfC count as 2 for g2
    expected = [1, 1, 2, 2]
    assert result['gene_potential'].tolist() == expected
    assert result['gene_id'].tolist() == ['g1', 'g1', 'g2', 'g2']
    


def test_gene_potential_single_gene(simple_df):
    """
    Single gene with 3 transcripts.
    Potential should be 3 for all rows.
    """
    # fabricate a feature_col (transcript_id) just for test
    df = simple_df.copy()
    df['transcript_id'] = ['t1', 't2', 't3']
    result = compute_gene_potential(df,
                                    gene_col='gene_id',
                                    feature_col='transcript_id')
    assert result['gene_potential'].tolist() == [3, 3, 3]


def test_gene_potential_empty_df():
    """
    Edge case: empty dataframe should return empty result.
    """
    df = pd.DataFrame(columns=[GENE_COL, FEATURE_COL])
    result = compute_gene_potential(df, gene_col=GENE_COL, feature_col=FEATURE_COL)
    # should still have the gene_potential column, but be empty
    assert 'gene_potential' in result.columns
    assert result.empty