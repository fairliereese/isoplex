#!/usr/bin/env python
import pytest
from iso_perplexity.utils import *

"""Tests for `iso_perplexity` package."""

# from iso_perplexity import iso_perplexity


# ------------------------------------------------------------
# Fixtures: build a minimal valid dataframe
# ------------------------------------------------------------
@pytest.fixture
def long_df():
    # 2 features × 3 samples
    return pd.DataFrame({
        "transcript_id": ["T1", "T1", "T1", "T2", "T2", "T2"],
        "sample":        ["S1", "S2", "S3", "S1", "S2", "S3"],
        "tpm":           [10, 20,  0,  5,  0,  15],
    })

@pytest.fixture
def multi_sample_counts_df():
    # Simple table with two genes, two isoforms, two samples
    return pd.DataFrame({
        "gene_id": ["G1", "G1", "G2", "G2"],
        "transcript_id": ["T1", "T2", "T3", "T4"],
        "S1": [100, 50, 10, 0],
        "S2": [80, 20, 5, 2],
    })


@pytest.fixture
def multi_sample_tpm_df():
    # Pretend these are TPMs
    return pd.DataFrame({
        "gene_id": ["G1", "G1", "G2", "G2"],
        "transcript_id": ["T1", "T2", "T3", "T4"],
        "S1": [60.0, 40.0, 80.0, 0.0],
        "S2": [50.0, 50.0, 70.0, 10.0],
    })



@pytest.fixture
def wide_df():
    return pd.DataFrame({
        "gene_id": ["g1", "g1", "g2"],
        "transcript_id": ["t1", "t2", "t3"],
        "s1": [10, 5, 3],
        "s2": [20, 7, 8],
    })

@pytest.fixture
def wide_valid_df():
    # 2 genes, 2 transcripts, 2 samples
    return pd.DataFrame({
        "gene_id": ["g1", "g1", "g2", "g2"],
        "transcript_id": ["t1", "t2", "t3", "t4"],
        "s1": [10, 5, 8, 3],
        "s2": [20, 0, 12, 6],
    })

@pytest.fixture
def manuscript_sample_df():
    df = pd.DataFrame()
    df['gene_id'] = ['A' for i in range(8*4)]
    df['transcript_id'] = flatten_list([[f'A_{j+1}' for i in range(4)] for j in range(8)])
    df['sample'] = flatten_list([['heart', 'brain', 'lungs', 'kidney'] for i in range(8)])
    df['pi'] = [0.5, 0.45, 0.2, 0.2,
                0, 0, 0.1, 0.2,
                0, 0, 0.08, 0.2,
                0.5, 0.1, 0.02, 0.2,
                0, 0, 0.5, 0,
                0, 0, 0.04, 0.2,
                0, 0.45, 0, 0,
                0, 0, 0.06, 0]

    # long formitze it
    df = df.pivot(index=['gene_id', 'transcript_id'], columns='sample', values='pi')
    df = df.reset_index()
    df.columns.name = ''

    return df

@pytest.fixture
def manuscript_df():
    df = pd.DataFrame()
    df['gene_id'] = ['A' for i in range(2)]+\
                    ['B' for i in range(7)]+\
                    ['C' for i in range(8)]
    df['transcript_id'] = [f'A_{i+1}' for i in range(2)]+\
                          [f'B_{i+1}' for i in range(7)]+\
                          [f'C_{i+1}' for i in range(8)]
    df['sample'] = [0.5, 0.5,
                 0.4, 0.35, 0.10, 0.07, 0.04, 0.03, 0.01]+\
                 [0.125 for i in range(8)]
    df['sample'] = [i*100 for i in df['sample'].tolist()]

    return df

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

@pytest.fixture
def simple_counts_df():
    """2 genes, 2 isoforms each, with varying counts."""
    return pd.DataFrame({
        "gene_id": ["G1", "G1", "G2", "G2"],
        "transcript_id": ["T1", "T2", "T3", "T4"],
        "global_counts": [100, 50, 30, 0],
    })


@pytest.fixture
def counts_with_feature():
    """A dataset where the feature column differs from transcript_id."""
    return pd.DataFrame({
        "gene_id": ["G1", "G1", "G1", "G1"],
        "transcript_id": ["T1", "T2", "T3", "T4"],
        "alt_feature": ["F1", "F1", "F2", "F2"],
        "global_counts": [10, 20, 0, 5],
    })


@pytest.fixture
def tpm_df(simple_counts_df):
    """Fake TPM dataset corresponding to the simple counts."""
    df = simple_counts_df.copy()
    df.drop('global_counts', axis=1, inplace=True)
    df["tpm"] = [60.0, 40.0, 80.0, 0.0]
    return df



############### validate_counts_input tests
# ------------------------------------------------------------
# 1. Passing test: everything is valid
# ------------------------------------------------------------
def test_validate_counts_input_valid(wide_valid_df):
    assert validate_counts_input(wide_valid_df,
                                 gene_col="gene_id",
                                 feature_col="transcript_id")


# ------------------------------------------------------------
# 2. Missing required columns
# ------------------------------------------------------------
def test_missing_columns(wide_valid_df):
    bad_df = wide_valid_df.drop(columns=["transcript_id"])
    with pytest.raises(KeyError, match="Missing required columns"):
        validate_counts_input(bad_df,
                              gene_col="gene_id",
                              feature_col="transcript_id")


# ------------------------------------------------------------
# 3. Null values in gene or feature columns
# ------------------------------------------------------------
def test_null_gene_values(wide_valid_df):
    bad_df = wide_valid_df.copy()
    bad_df.loc[0, "gene_id"] = None
    with pytest.raises(ValueError, match="Missing values in `gene_id`"):
        validate_counts_input(bad_df,
                              gene_col="gene_id",
                              feature_col="transcript_id")


def test_null_feature_values(wide_valid_df):
    bad_df = wide_valid_df.copy()
    bad_df.loc[0, "transcript_id"] = None
    with pytest.raises(ValueError, match="Missing values in `transcript_id`"):
        validate_counts_input(bad_df,
                              gene_col="gene_id",
                              feature_col="transcript_id")

# ------------------------------------------------------------
# 4. Negative expression values
# ------------------------------------------------------------
def test_negative_expression(wide_valid_df):
    bad_df = wide_valid_df.copy()
    bad_df.loc[0, "s1"] = -5
    with pytest.raises(ValueError, match="negative"):
        validate_counts_input(bad_df,
                              gene_col="gene_id",
                              feature_col="transcript_id")

def test_all_zero_expression(wide_valid_df):
    bad_df = wide_valid_df.copy()
    bad_df[["s1", "s2"]] = 0
    with pytest.raises(ValueError, match="sum to zero"):
        validate_counts_input(bad_df,
                              gene_col="gene_id",
                              feature_col="transcript_id")

# ------------------------------------------------------------
# 4.5. Negative expression values
# ------------------------------------------------------------
def test_negative_expression(wide_valid_df):
    bad_df = wide_valid_df.copy()
    bad_df.loc[0, "s1"] = -5
    with pytest.raises(ValueError, match="negative"):
        validate_counts_input(bad_df,
                              gene_col="gene_id",
                              feature_col="transcript_id")


# ------------------------------------------------------------
# 6. Extra non-numeric metadata column
# ------------------------------------------------------------
def test_extra_non_numeric_column(wide_valid_df):
    bad_df = wide_valid_df.copy()
    bad_df["metadata"] = ["meta1", "meta2", "meta3", "meta4"]
    with pytest.raises(ValueError, match="non-numeric columns"):
        validate_counts_input(bad_df,
                              gene_col="gene_id",
                              feature_col="transcript_id")

# ########### wide to long tests
# def test_wide_to_long_basic(wide_df):
#     out = wide_to_long(wide_df,
#                        gene_col="gene_id",
#                        feature_col="transcript_id",
#                        expression_type="counts")

#     # expected shape: rows = original_rows * num_samples
#     assert out.shape == (wide_df.shape[0] * 2, 4)
#     assert set(out.columns) == {"gene_id", "transcript_id", "sample", "counts"}

#     # check a few known values
#     first_row = out.iloc[0]
#     assert first_row["gene_id"] == "g1"
#     assert first_row["transcript_id"] == "t1"
#     assert first_row["sample"] in ["s1", "s2"]
#     assert first_row["counts"] in [10, 20]

# # ------------------------------------------------------------
# # 2. Works with different expression column name
# # ------------------------------------------------------------
# def test_wide_to_long_with_tpm(wide_df):
#     out = wide_to_long(wide_df,
#                        gene_col="gene_id",
#                        feature_col="transcript_id",
#                        expression_type="tpm")
#     assert "tpm" in out.columns
#     assert "counts" not in out.columns

################# collapse_counts_by_feature tests
def test_collapse_single_sample(single_sample_df):
    out = collapse_counts_by_feature(
        single_sample_df,
        feature_col="orf_id",
        expression_type="counts",
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
        expression_type="counts",
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
        expression_type="counts",
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
    valid_df.rename({EXP_COL:'tpm'}, axis=1, inplace=True) # filter here requires tpm col
    result = compute_gene_potential(valid_df, gene_col=GENE_COL, feature_col=FEATURE_COL)
    expected = [2, 2, 2, 2]   # each row inherits its gene's potential
    assert result['gene_potential'].tolist() == expected


def test_gene_potential_collapsed(single_sample_df):
    """
    Mixed situation: two genes.
    g1 has only one unique ORF (orfA)
    g2 has two unique ORFs (orfB, orfC)
    """
    single_sample_df.rename({'counts':'tpm'}, axis=1, inplace=True)
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
    df = pd.DataFrame(columns=[GENE_COL, FEATURE_COL, 'tpm'])
    result = compute_gene_potential(df, gene_col=GENE_COL, feature_col=FEATURE_COL)
    # should still have the gene_potential column, but be empty
    assert 'gene_potential' in result.columns
    assert result.empty

########### test compute_entropy
def test_entropy_single_gene_uniform(simple_df):
    """
    For a single gene with equal isoform proportions,
    entropy should be log2(n), where n is number of isoforms.
    """
    # uniform proportions: each pi = 1/3
    df = simple_df.copy()
    df['pi'] = [1/3, 1/3, 1/3]

    result = compute_entropy(df, gene_col=GENE_COL)

    # expected entropy = log2(3)
    expected_entropy = np.log2(3)

    # all rows should inherit the same gene-level entropy
    assert np.allclose(result['entropy'].values, expected_entropy, rtol=1e-8)


def test_entropy_single_gene_skewed(simple_df):
    """
    Skewed distribution:
    pi = [0.5, 0.25, 0.25]
    Compute entropy explicitly.
    """
    df = simple_df.copy()
    df['pi'] = [0.5, 0.25, 0.25]

    result = compute_entropy(df, gene_col=GENE_COL)

    expected = -(
        0.5 * np.log2(0.5) +
        0.25 * np.log2(0.25) +
        0.25 * np.log2(0.25)
    )

    assert np.allclose(result['entropy'].values, expected, rtol=1e-8)


def test_entropy_multiple_genes(multi_gene_df):
    """
    Multiple genes with different isoform distributions.
    Each gene's entropy should be computed separately.
    """
    df = multi_gene_df.copy()
    # For gene g1: proportions 10, 20 → pi = [1/3, 2/3]
    # For gene g2: proportions 5, 15 → pi = [0.25, 0.75]
    df['pi'] = [1/3, 2/3, 0.25, 0.75]

    result = compute_entropy(df, gene_col=GENE_COL)

    # expected entropy per gene
    H_g1 = -(1/3 * np.log2(1/3) + 2/3 * np.log2(2/3))
    H_g2 = -(0.25 * np.log2(0.25) + 0.75 * np.log2(0.75))

    expected = [H_g1, H_g1, H_g2, H_g2]
    assert np.allclose(result['entropy'].values, expected, rtol=1e-8)


def test_entropy_with_zero_pi():
    """
    Zero pi values should not cause math errors.
    Entropy contribution for zero pi is defined as zero.
    """
    df = pd.DataFrame({
        GENE_COL: ['g1', 'g1', 'g1'],
        'pi': [0.0, 0.5, 0.5]
    })

    result = compute_entropy(df, gene_col=GENE_COL)

    expected = -(0.5 * np.log2(0.5) + 0.5 * np.log2(0.5))

    assert np.allclose(result['entropy'].values, expected, rtol=1e-8)


############# test compute_perplexity
def test_perplexity_basic():
    """
    Perplexity is 2^entropy.
    Check for simple, known values.
    """
    df = pd.DataFrame({
        GENE_COL: ['g1', 'g1', 'g2', 'g2'],
        'entropy': [1.0, 1.0, 2.0, 2.0]
    })

    result = compute_perplexity(df)

    # expected: g1 -> 2^1 = 2; g2 -> 2^2 = 4
    expected = [2.0, 2.0, 4.0, 4.0]

    assert np.allclose(result['perplexity'].values, expected, rtol=1e-8)


def test_perplexity_zero_entropy():
    """
    When entropy is zero, perplexity should be 1.
    """
    df = pd.DataFrame({
        GENE_COL: ['g1', 'g1', 'g2', 'g2'],
        'entropy': [0.0, 0.0, 0.0, 0.0]
    })

    result = compute_perplexity(df)

    assert np.allclose(result['perplexity'].values, 1.0, rtol=1e-8)


def test_perplexity_floats():
    """
    Works with non-integer entropy.
    """
    df = pd.DataFrame({
        GENE_COL: ['g1', 'g1', 'g2', 'g2'],
        'entropy': [0.5, 0.5, 1.5, 1.5]
    })

    result = compute_perplexity(df)

    expected = [2 ** 0.5, 2 ** 0.5, 2 ** 1.5, 2 ** 1.5]

    assert np.allclose(result['perplexity'].values, expected, rtol=1e-8)

################### test call_effective
def test_call_effective_basic():
    """
    Basic case:
    - perplexity = 2 → two top-ranked features marked as effective.
    """
    df = pd.DataFrame({
        GENE_COL: ['g1', 'g1', 'g1', 'g2', 'g2'],
        TRANSCRIPT_COL: ['t1', 't2', 't3', 't1', 't2'],
        'pi': [0.5, 0.3, 0.2, 0.7, 0.3],
        'perplexity': [2.0, 2.0, 2.0, 1.0, 1.0],
    })

    result = call_effective(df)

    # perplexity = 2 → top 2 isoforms effective for g1
    g1 = result[result[GENE_COL] == 'g1']
    assert set(g1[g1['effective']][TRANSCRIPT_COL]) == {'t1', 't2'}

    # perplexity = 1 → only top isoform effective for g2
    g2 = result[result[GENE_COL] == 'g2']
    assert set(g2[g2['effective']][TRANSCRIPT_COL]) == {'t1'}

    # check rank ordering
    assert all(g1.sort_values('pi', ascending=False)['feature_rank'].values == [1, 2, 3])
    assert all(g2.sort_values('pi', ascending=False)['feature_rank'].values == [1, 2])


def test_call_effective_rounding():
    """
    Check that perplexity is rounded to the nearest integer.
    """
    df = pd.DataFrame({
        GENE_COL: ['g1', 'g1', 'g1'],
        TRANSCRIPT_COL: ['t1', 't2', 't3'],
        'pi': [0.5, 0.3, 0.2],
        'perplexity': [1.6, 1.6, 1.6],  # rounds to 2
    })

    result = call_effective(df)

    # expect 2 effective features due to rounding
    assert (result['n_effective'].unique() == [2]).all()
    assert sum(result['effective']) == 2


def test_call_effective_ties():
    """
    Ties in pi are handled by 'first' method.
    """
    df = pd.DataFrame({
        GENE_COL: ['g1', 'g1', 'g1'],
        TRANSCRIPT_COL: ['t1', 't2', 't3'],
        'pi': [0.4, 0.4, 0.4],
        'perplexity': [2.0, 2.0, 2.0],
    })

    result = call_effective(df)

    # Top 2 features by order of appearance get marked effective
    assert sum(result['effective']) == 2
    assert set(result[result['effective']][TRANSCRIPT_COL]) <= {'t1', 't2'}


def test_call_effective_all_effective():
    """
    If perplexity >= number of features, all features should be effective.
    """
    df = pd.DataFrame({
        GENE_COL: ['g1', 'g1', 'g1'],
        TRANSCRIPT_COL: ['t1', 't2', 't3'],
        'pi': [0.5, 0.3, 0.2],
        'perplexity': [3.5, 3.5, 3.5],  # rounds to 4
    })

    result = call_effective(df)

    assert all(result['effective'])

################## test calc_expressino_breadth

def test_expression_breadth_basic():
    """
    Basic case: some features effective in all samples, some in subset, some in none.
    """
    df = pd.DataFrame({
        TRANSCRIPT_COL: ['t1', 't1', 't2', 't2', 't3', 't3'],
        SAMPLE_COL:     ['s1', 's2', 's1', 's2', 's1', 's2'],
        'effective': [True, True, True, False, False, False]
    })

    result = compute_expression_breadth(df, SAMPLE_COL, TRANSCRIPT_COL)

    # There are 2 samples total
    # t1 effective in both samples → breadth = 100%
    # t2 effective in only 1 sample → breadth = 50%
    # t3 effective in 0 samples    → breadth = 0%
    breadth = dict(zip(result[TRANSCRIPT_COL], result['expression_breadth']))

    assert breadth['t1'] == 100
    assert breadth['t2'] == 50
    assert breadth['t3'] == 0

    n_eff = dict(zip(result[TRANSCRIPT_COL], result['n_samples_effective']))
    assert n_eff['t1'] == 2
    assert n_eff['t2'] == 1
    assert pd.isna(n_eff['t3']) or n_eff['t3'] == 0  # after fillna(0)


def test_expression_breadth_all_none():
    """
    All features ineffective → all breadth = 0%.
    """
    df = pd.DataFrame({
        TRANSCRIPT_COL: ['t1', 't1', 't2', 't2'],
        SAMPLE_COL:     ['s1', 's2', 's1', 's2'],
        'effective': [False, False, False, False]
    })

    result = compute_expression_breadth(df, SAMPLE_COL, TRANSCRIPT_COL)

    assert (result['expression_breadth'] == 0).all()
    assert (result['n_samples_effective'].fillna(0) == 0).all()


def test_expression_breadth_all_effective():
    """
    All features effective in all samples → all breadth = 100%.
    """
    df = pd.DataFrame({
        TRANSCRIPT_COL: ['t1', 't1', 't2', 't2'],
        SAMPLE_COL:     ['s1', 's2', 's1', 's2'],
        'effective': [True, True, True, True]
    })

    result = compute_expression_breadth(df, SAMPLE_COL, TRANSCRIPT_COL)

    assert (result['expression_breadth'] == 100).all()
    assert (result['n_samples_effective'] == 2).all()


def test_expression_breadth_single_sample():
    """
    Edge case: only one sample → breadth should be 0% or 100%.
    """
    df = pd.DataFrame({
        TRANSCRIPT_COL: ['t1', 't2', 't3'],
        SAMPLE_COL:     ['s1', 's1', 's1'],
        'effective': [True, False, True]
    })

    result = compute_expression_breadth(df, SAMPLE_COL, TRANSCRIPT_COL)

    # Only one sample: breadth is either 0 or 100
    breadth = dict(zip(result[TRANSCRIPT_COL], result['expression_breadth']))
    assert breadth['t1'] == 100
    assert breadth['t2'] == 0
    assert breadth['t3'] == 100

############ expression variance tests

def test_expression_var_basic():
    """
    Basic case with two features across two samples.
    """
    df = pd.DataFrame({
        TRANSCRIPT_COL: ['t1', 't1', 't2', 't2'],
        SAMPLE_COL:     ['s1', 's2', 's1', 's2'],
        'pi':           [0.4, 0.6, 0.1, 0.3]
    })

    result = compute_expression_var(df, SAMPLE_COL, TRANSCRIPT_COL)

    # Each feature appears in 2 samples
    assert (result['n_exp_samples'] == 2).all()

    # Expected standard deviations
    # t1: std of [0.4, 0.6]
    expected_t1_std = np.std([0.4, 0.6], ddof=1)
    # t2: std of [0.1, 0.3]
    expected_t2_std = np.std([0.1, 0.3], ddof=1)

    # Check that computed std matches expected
    std_by_feature = dict(zip(result[TRANSCRIPT_COL], result['expression_var']))
    assert std_by_feature['t1'] == pytest.approx(expected_t1_std, rel=1e-8)
    assert std_by_feature['t2'] == pytest.approx(expected_t2_std, rel=1e-8)


def test_expression_var_single_sample():
    """
    If a feature is observed in only one sample, std should be NaN.
    """
    df = pd.DataFrame({
        TRANSCRIPT_COL: ['t1', 't2'],
        SAMPLE_COL:     ['s1', 's1'],
        'pi':           [0.5, 0.7]
    })

    result = compute_expression_var(df, SAMPLE_COL, TRANSCRIPT_COL)

    # Each feature is in 1 sample → std = NaN
    assert (result['n_exp_samples'] == 1).all()
    assert result['expression_var'].isna().all()


def test_expression_var_some_missing_pi():
    """
    Missing pi values should be skipped in std calculation.
    """
    df = pd.DataFrame({
        TRANSCRIPT_COL: ['t1', 't1', 't1'],
        SAMPLE_COL:     ['s1', 's2', 's3'],
        'pi':           [0.2, np.nan, 0.4]
    })

    result = compute_expression_var(df, SAMPLE_COL, TRANSCRIPT_COL)

    # Feature appears in 3 samples
    assert (result['n_exp_samples'] == 3).all()

    # std should be computed only from non-NaN values [0.2, 0.4]
    expected_std = np.std([0.2, 0.4], ddof=1)
    assert result['expression_var'].iloc[0] == pytest.approx(expected_std, rel=1e-8)


def test_expression_var_feature_with_constant_pi():
    """
    Zero variance case → std should be 0.
    """
    df = pd.DataFrame({
        TRANSCRIPT_COL: ['t1', 't1', 't1'],
        SAMPLE_COL:     ['s1', 's2', 's3'],
        'pi':           [0.5, 0.5, 0.5]
    })

    result = compute_expression_var(df, SAMPLE_COL, TRANSCRIPT_COL)

    assert (result['n_exp_samples'] == 3).all()
    assert result['expression_var'].iloc[0] == pytest.approx(0.0, rel=1e-8)

############ tests for compute_average-expression

@pytest.fixture
def simple_multi_sample_df():
    """
    Two features across three samples.
    t1 is expressed in all samples.
    t2 is expressed only in one sample.
    """
    return pd.DataFrame({
        'transcript_id': ['t1', 't1', 't1', 't2', 't2', 't2'],
        'sample':        ['s1', 's2', 's3', 's1', 's2', 's3'],
        'tpm':        [10, 20, 30, 5, 0, 0]
    })


@pytest.fixture
def duplicated_rows_df():
    """
    A feature appears twice in the same sample (e.g. technical replicates)
    which should be summed before averaging.
    """
    return pd.DataFrame({
        'transcript_id': ['t1', 't1', 't1', 't1', 't2', 't2'],
        'sample':        ['s1', 's1', 's2', 's3', 's1', 's2'],
        'tpm':        [5, 5, 20, 30, 10, 15]  # t1 has duplicate rows in s1
    })


@pytest.fixture
def all_unexpressed_df():
    """
    All counts are zero; average expression should be NaN.
    """
    return pd.DataFrame({
        'transcript_id': ['t1', 't1', 't2', 't2'],
        'sample':        ['s1', 's2', 's1', 's2'],
        'tpm':        [0, 0, 0, 0]
    })

# --- Tests ---

def test_basic_avg_expression(simple_multi_sample_df):
    """
    Average expression excludes zero-count samples.
    """

    df = compute_avg_expression(simple_multi_sample_df, sample_col='sample', feature_col='transcript_id')

    # For t1: (10+20+30)/3 = 20
    assert df.loc[df['transcript_id'] == 't1', 'avg_transcript_id_tpm'].iloc[0] == 20

    # For t2: only s1 is expressed, so mean = 5 (not 5/3)
    assert df.loc[df['transcript_id'] == 't2', 'avg_transcript_id_tpm'].iloc[0] == 5

def test_all_unexpressed(all_unexpressed_df):
    """
    If a feature is never expressed, the average should be NaN.
    """
    df = compute_avg_expression(all_unexpressed_df, sample_col='sample', feature_col='transcript_id')

    # Both features are unexpressed -> NaN
    assert df['avg_transcript_id_tpm'].isna().all()


def test_merge_preserves_rows(simple_multi_sample_df):
    """
    Ensure no rows are lost after merging back.
    """
    df = compute_avg_expression(simple_multi_sample_df, sample_col='sample', feature_col='transcript_id')

    assert len(df) == len(simple_multi_sample_df), "Row count changed after merging"
    assert 'avg_transcript_id_tpm' in df.columns

########## max expression tests
def test_basic_max_expression(simple_multi_sample_df):
    """
    The max expression for each feature should be the maximum
    TPM among non-zero samples.
    """
    df = compute_max_expression(simple_multi_sample_df, sample_col='sample', feature_col='transcript_id')

    # For t1: max across s1=10, s2=20, s3=30 → 30
    assert df.loc[df['transcript_id'] == 't1', 'max_transcript_id_tpm'].iloc[0] == 30

    # For t2: max across s1=5, s2=0, s3=0 → 5
    assert df.loc[df['transcript_id'] == 't2', 'max_transcript_id_tpm'].iloc[0] == 5

def test_duplicate_rows_handled(duplicated_rows_df):
    """
    Duplicated entries for the same feature+sample
    should be summed before taking the max.
    """

    df = compute_max_expression(duplicated_rows_df, sample_col='sample', feature_col='transcript_id')

    # t1 has s1: 5+5=10, s2:20, s3:30 → max=30
    assert df.loc[df['transcript_id'] == 't1', 'max_transcript_id_tpm'].iloc[0] == 30

    # t2 has s1=10, s2=15 → max=15
    assert df.loc[df['transcript_id'] == 't2', 'max_transcript_id_tpm'].iloc[0] == 15


def test_all_unexpressed(all_unexpressed_df):
    """
    If a feature is never expressed (>0), max should be NaN.
    """
    df = compute_max_expression(all_unexpressed_df, sample_col='sample', feature_col='transcript_id')

    assert df['max_transcript_id_tpm'].isna().all()


def test_merge_preserves_rows(simple_multi_sample_df):
    """
    Ensure the row count and order are preserved after merging back.
    """
    df = compute_max_expression(simple_multi_sample_df, sample_col='sample', feature_col='transcript_id')

    assert len(df) == len(simple_multi_sample_df)
    assert 'max_transcript_id_tpm' in df.columns


############# testing compute_global_isoform_metrics
def test_happy_path_counts(simple_counts_df):
    df_out = compute_global_isoform_metrics(simple_counts_df,
                                            expression_type="counts")

    # Expect same rows as input
    assert len(df_out) == len(simple_counts_df)

    # Columns we expect to be present
    expected_cols = {"gene_id", "transcript_id", "counts", "tpm",
                     "pi", "gene_potential", "entropy", "perplexity", "effective"}
    assert expected_cols.issubset(df_out.columns)

    # Check that TPM was computed
    assert not df_out["tpm"].isna().any()


def test_happy_path_tpm(tpm_df):

    df_out = compute_global_isoform_metrics(tpm_df,
                                            expression_type='tpm')

    # Should not overwrite supplied TPM
    assert np.allclose(df_out["tpm"], tpm_df["tpm"], equal_nan=True)


def test_invalid_input_missing_expression(simple_counts_df):
    df_bad = simple_counts_df.drop(columns=["global_counts"])
    with pytest.raises(ValueError, match="No sample"):
        compute_global_isoform_metrics(df_bad,
                                       expression_type="counts")

def test_invalid_input_negative_counts(simple_counts_df):
    df_bad = simple_counts_df.copy()
    df_bad.loc[0, "global_counts"] = -10
    with pytest.raises(ValueError, match="negative"):
        compute_global_isoform_metrics(df_bad,
                                       expression_type="counts")

# def test_zero_counts_behavior(simple_counts_df):
#     df = simple_counts_df.copy()
#     df["counts"] = [0, 0, 0, 0]
#     df_out = compute_global_isoform_metrics(df,
#                                             expression_type="counts",
#                                             expression_type_type="counts")
#     # All TPM and pi should be zero or NaN-safe
#     assert (df_out["tpm"] == 0).all() or (df_out["tpm"].fillna(0) == 0).all()


def test_collapsing_feature(counts_with_feature):
    counts_with_feature = counts_with_feature.drop('transcript_id', axis=1)
    df_out = compute_global_isoform_metrics(counts_with_feature,
                                            gene_col="gene_id",
                                            feature_col="alt_feature",
                                            expression_type="counts")
    # Expect output to have unique alt_feature per gene
    assert df_out["alt_feature"].nunique() == 2


# # def test_unexpected_sample_column(simple_counts_df):
# #     df = simple_counts_df.copy()
# #     df["sample"] = ["s1", "s1", "s1", "s1"]
# #     # Should error because this function is single-sample only
# #     with pytest.raises(ValueError, match="single-sample"):
# #         compute_global_isoform_metrics(df,
# #                                        expression_type="counts",
# #                                        expression_type_type="counts")


def test_order_independence(simple_counts_df):
    shuffled = simple_counts_df.sample(frac=1, random_state=42)
    out1 = compute_global_isoform_metrics(simple_counts_df,
                                          expression_type="counts")
    out2 = compute_global_isoform_metrics(shuffled,
                                          expression_type="counts")

    # Sort both outputs for comparison
    out1_sorted = out1.sort_values(["gene_id", "transcript_id"]).reset_index(drop=True)
    out2_sorted = out2.sort_values(["gene_id", "transcript_id"]).reset_index(drop=True)

    pd.testing.assert_frame_equal(out1_sorted[out1_sorted.columns],
                                  out2_sorted[out2_sorted.columns])

def test_manuscript_global(manuscript_df):
    df = compute_global_isoform_metrics(manuscript_df,
                                          expression_type="counts")

    # gene-level metrics
    assert df.set_index('gene_id')['gene_potential'].to_dict() == {'A': 2, 'B': 7, 'C': 8}
    assert df.set_index('gene_id')['entropy'].to_dict() == pytest.approx({'A': 1, 'B': 2.06, 'C': 3}, rel=1e6)
    assert df.set_index('gene_id')['perplexity'].to_dict() == pytest.approx({'A': 2, 'B': 4.18, 'C': 8}, rel=1e6)

    # efective isoforms
    eff_isos = df.loc[df.effective].transcript_id.tolist()
    truth = ['A_1', 'A_2', 'B_1', 'B_2', 'B_3', 'B_4',
             'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6', 'C_7', 'C_8']
    assert eff_isos == truth

############## testing multi sample workflow
def test_multi_sample_happy_path(multi_sample_counts_df):
    out = compute_multi_sample_isoform_metrics(
        multi_sample_counts_df,
        gene_col="gene_id",
        feature_col="transcript_id",
        expression_type="counts"
    )

    # Expect 2 samples × 4 isoforms = 8 rows
    assert len(out) == 8

    # There should be a column for sample ID
    assert "sample" in out.columns
    assert set(out["sample"]) == {"S1", "S2"}

    # Per-sample metrics should exist
    for col in ["tpm", "pi", "gene_potential", "entropy", "perplexity", "effective"]:
        assert col in out.columns

    # Cross-sample metrics should exist
    assert any("breadth" in c for c in out.columns)
    assert any("expression_var" in c for c in out.columns)
    assert any("avg_transcript_id_tpm" in c for c in out.columns)
    assert any("avg_gene_id_tpm" in c for c in out.columns)
    assert any("max_gene_id_tpm" in c for c in out.columns)
    assert any("max_transcript_id_tpm" in c for c in out.columns)



def test_multi_sample_tpm(multi_sample_tpm_df):

    out = compute_multi_sample_isoform_metrics(
        multi_sample_tpm_df,
        gene_col="gene_id",
        feature_col="transcript_id",
        expression_type="tpm"
    )


    # Should use provided TPM directly, not recompute
    assert np.allclose(
        out[out["sample"] == "S1"].sort_values("transcript_id")["tpm"].values,
        multi_sample_tpm_df["S1"].values
    )

def test_multi_sample_missing_sample(multi_sample_counts_df):
    bad = multi_sample_counts_df.drop(columns=["S1", "S2"])
    with pytest.raises(Exception):   # ValueError/KeyError depending on validator
        compute_multi_sample_isoform_metrics(
            bad,
            gene_col="gene_id",
            feature_col="transcript_id",
            expression_type="counts"
        )


def test_multi_sample_negative_counts(multi_sample_counts_df):
    bad = multi_sample_counts_df.copy()
    bad.loc[0, "S1"] = -5
    with pytest.raises(ValueError, match="negative"):
        compute_multi_sample_isoform_metrics(
            bad,
            gene_col="gene_id",
            feature_col="transcript_id",
            expression_type="counts"
        )

def test_multi_sample_order_independence(multi_sample_counts_df):
    shuffled = multi_sample_counts_df.sample(frac=1, random_state=42)

    out1 = compute_multi_sample_isoform_metrics(
        multi_sample_counts_df,
        gene_col="gene_id",
        feature_col="transcript_id",
        expression_type="counts"
    )
    out2 = compute_multi_sample_isoform_metrics(
        shuffled,
        gene_col="gene_id",
        feature_col="transcript_id",
        expression_type="counts"
    )

    # Sort for comparison
    cols_to_check = ["gene_id", "transcript_id", "sample", "tpm"]
    out1s = out1[cols_to_check].sort_values(cols_to_check).reset_index(drop=True)
    out2s = out2[cols_to_check].sort_values(cols_to_check).reset_index(drop=True)

    pd.testing.assert_frame_equal(out1s, out2s)

# now the most important test is if the manuscript output is correct
def test_manuscript_sample(manuscript_sample_df):
    df = compute_multi_sample_isoform_metrics(manuscript_sample_df,
                                              gene_col='gene_id',
                                              feature_col='transcript_id',
                                              expression_type='counts')

    # gene-level metrics
    assert df.set_index('sample')['gene_potential'].to_dict() == {'heart': 2, 'brain': 3, 'lungs': 7, 'kidney': 5}

    # effective isoforms

    # sample level metrics
    assert df.set_index('sample')['perplexity'].to_dict() == pytest.approx({'heart':2, 'brain':2.58, 'lungs':4.38, 'kidney':5}, rel=1e6)

    # tissue-level metrics
    assert df.set_index('transcript_id')['expression_breadth'].to_dict() == pytest.approx({'A_1':100, 'A_2':50, 'A_3':50, 'A_4':75, 'A_5':25, 'A_6':25, 'A_7':25, 'A_8':0}, rel=1e6)
    assert df.set_index('transcript_id')['expression_breadth'].to_dict() == pytest.approx({'A_1':0.160, 'A_2':0.096, 'A_3':0.094, 'A_4':0.210, 'A_5':0.250, 'A_6':0.095, 'A_7':0.225, 'A_8':0.030}, rel=1e6)
