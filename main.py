# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    species = [
        ("Gallus gallus", gg_seq),
        ("Mus musculus", mm_seq),
        ("Balaeniceps rex", br_seq),
        ("Tursiops truncatus", tt_seq)
    ]
    
    # Align each species to human BRD2 and store scores
    alignment_results = []
    for name, seq in species:
        aligner = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
        score, aligned_hs, aligned_seq = aligner.align(hs_seq, seq)
        alignment_results.append((name, score, aligned_hs, aligned_seq))
    
    # Sort species by alignment score (most similar first)
    alignment_results.sort(key=lambda x: x[1], reverse=True)
    
    print("Species in order of most similar to human BRD2:")
    for name, score, _, _ in alignment_results:
        print(f"{name}: {score}")
    
    print("\nPairwise alignment scores and sequences:")
    for name, score, aligned_hs, aligned_seq in alignment_results:
        print(f"Human BRD2 vs {name}: {score}")
        print(f"Human BRD2 Aligned:   {aligned_hs}")
        print(f"{name} Aligned:       {aligned_seq}")
        print()  # For better readability
    

if __name__ == "__main__":
    main()
