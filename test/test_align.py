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
    nw = NeedlemanWunsch("BLOSUM62.txt", gap_open=-10, gap_extend=-1)

    score, aligned_seqA, aligned_seqB = nw.align(seq1, seq2)
    # expected_score
    # expected_aligned_seqA
    # expected_aligned_seqB

    # assert score == expected_score, f"Expected {expected_score}, but got {score}"
    # assert aligned_seqA == expected_aligned_seqA, f"Expected {expected_aligned_seqA}, but got {aligned_seqA}"
    # assert aligned_seqB == expected_aligned_seqB, f"Expected {expected_aligned_seqB}, but got {aligned_seqB}"


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
    
    nw = NeedlemanWunsch("BLOSUM62.txt", gap_open=-10, gap_extend=-1)

    score, aligned_seqA, aligned_seqB = nw.align(seq3, seq4)

    expected_score = 17
    expected_aligned_seqA = 'MAVHQLIRRP'
    expected_aligned_seqB = 'M---QLIRHP'

    assert score == expected_score, f"Expected {expected_score}, but got {score}"
    assert aligned_seqA == expected_aligned_seqA, f"Expected {expected_aligned_seqA}, but got {aligned_seqA}"
    assert aligned_seqB == expected_aligned_seqB, f"Expected {expected_aligned_seqB}, but got {aligned_seqB}"