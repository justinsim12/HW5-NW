# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    # Example usage with BLOSUM62 matrix and affine gaps
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    # Run the alignment
    score, aligned_seqA, aligned_seqB = nw.align(seq1, seq2)

    # Expected Alignment Matrix
    expected_align_matrix = np.array([
        [  0., -10., -11., -12.,],
        [-10.,   5.,  -6.,  -7.,],
        [-11.,  -6.,   4.,  -7.,],
        [-12.,  -7.,  -1.,   5.,],
        [-13.,  -8.,  -6.,   4.,]])

    # Expected Gap Matrices
    expected_gapA_matrix = np.array([
        [  0, -10, -11, -12],
        [-10,  -1, -11, -12],
        [-11,  -2,  -2, -12],
        [-12,  -3,  -3,  -3],
        [-13,  -4,  -4,  -4]
    ])

    expected_gapB_matrix = np.array([
        [  0, -10, -11, -12],
        [-10,  -1,  -2,  -3],
        [-11, -11,  -2,  -3],
        [-12, -12,  -3,  -4],
        [-13, -13,  -4,  -4]
    ])

    # Expected Backtrace Matrix
    expected_backtrace = np.array([
        ['' , 'L', 'L', 'L'],
        ['U', 'D', 'L', 'L'],
        ['U', 'U', 'D', 'L'],
        ['U', 'U', 'D', 'D'],
        ['U', 'U', 'U', 'D']
    ])

    # Assert Matrices
    np.testing.assert_array_almost_equal(nw._align_matrix, expected_align_matrix, decimal=3)
    np.testing.assert_array_almost_equal(nw._gapA_matrix, expected_gapA_matrix, decimal=3)
    np.testing.assert_array_almost_equal(nw._gapB_matrix, expected_gapB_matrix, decimal=3)
    assert np.array_equal(nw._back, expected_backtrace)

    # Assert Final Score
    assert score == 4, f"Expected alignment score 4, but got {score}"

    # Assert Correct Alignment
    assert aligned_seqA == "MYQR"
    assert aligned_seqB == "M-QR"

    
def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    score, aligned_seqA, aligned_seqB = nw.align(seq3, seq4)

    expected_score = 17
    expected_aligned_seqA = 'MAVHQLIRRP'
    expected_aligned_seqB = 'M---QLIRHP'

    assert score == expected_score, f"Expected {expected_score}, but got {score}"
    assert aligned_seqA == expected_aligned_seqA, f"Expected {expected_aligned_seqA}, but got {aligned_seqA}"
    assert aligned_seqB == expected_aligned_seqB, f"Expected {expected_aligned_seqB}, but got {aligned_seqB}"