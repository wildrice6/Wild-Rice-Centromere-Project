#!/usr/bin/env python3
import argparse
from itertools import combinations
import hashlib
import math


def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = ""
            else:
                seqs[name] += line.upper()
    return seqs


def hash_kmer(kmer):
    """Stable integer hash for a k-mer."""
    return int(hashlib.md5(kmer.encode()).hexdigest(), 16)


def get_modimizers(seq, k=7, s=10):
    """
    Mod-sampling: keep k-mers whose hash % s == 0.
    """
    mods = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if hash_kmer(kmer) % s == 0:
            mods.add(kmer)
    return mods


def cmod_similarity(seq1, seq2, k=7, s=10):
    """
    Compute ModDotPlot Cmod:
        Cmod = (|A ∩ B| / |A|) / (1 - (1 - 1/s)^(|A_kmers|))
    """
    A = get_modimizers(seq1, k, s)
    B = get_modimizers(seq2, k, s)

    if len(A) == 0:
        return 0.0

    raw_containment = len(A & B) / len(A)

    # number of k-mers in seq1
    nA = len(seq1) - k + 1
    if nA <= 0:
        return 0.0

    # bias correction term from ModDotPlot (formula 4)
    correction = 1 - (1 - 1/s) ** nA

    if correction == 0:
        return 0.0

    return raw_containment / correction


def cmod_to_ani(Cmod, k):
    """
    Convert Cmod to ANI using Mash equation:
        ANI = 1 + ln(Cmod) / k
    """
    if Cmod <= 0:
        return 0.0
    return 1 + math.log(Cmod) / k


def main():
    parser = argparse.ArgumentParser(description="Compute ModDotPlot Cmod + ANI")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--k", type=int, default=7, help="k-mer size (default: 7)")
    parser.add_argument("--s", type=int, default=10, help="mod-sampling sparsity (default: 10)")
    args = parser.parse_args()

    seqs = read_fasta(args.input)
    names = list(seqs.keys())

    with open(args.output, "w") as out:
        out.write("seq1\tseq2\tCmod\tANI\n")
        for a, b in combinations(names, 2):
            Cmod = cmod_similarity(seqs[a], seqs[b], k=args.k, s=args.s)
            ANI = cmod_to_ani(Cmod, args.k)
            out.write(f"{a}\t{b}\t{Cmod:.6f}\t{ANI:.6f}\n")


if __name__ == "__main__":
    main()


