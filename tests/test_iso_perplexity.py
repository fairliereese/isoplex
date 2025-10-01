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
