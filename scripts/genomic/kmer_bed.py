import pandas as pd
import numpy as np
import dinopy
import pyBigWig

from tqdm import tqdm

from collections import deque

K_MER_LENGTH = 10

def iter_kmers(iterator, k=10):

    kmers = deque(maxlen=k)

    for letter in iterator:
        kmers.append(letter)

        if len(kmers) == k:
            yield ''.join(kmers)

def get_kmer_dict(kmer_score_h5, motif_id, score_column, nlargest):
    """
        Reads kmer scores from hdf5
    """
    kmer_scores = pd.read_hdf(kmer_score_h5, score_column)[motif_id]
    kmer_scores = kmer_scores.nlargest(nlargest)
    kmer_scores = kmer_scores.sort_values(ascending=False)

    ans = {}
    for i, kmer in enumerate(kmer_scores.index, start=1):
        ans[kmer] = f'{i}_{kmer}'

    return ans

def make_bigwig(genomic_fasta,
                kmer_score_h5,
                motif_id,
                score_column,
                nlargest,
                output_file):

    kmer_dict = get_kmer_dict(kmer_score_h5, motif_id,
                              score_column, nlargest)

    far = dinopy.FastaReader(genomic_fasta)

    with open(output_file, 'w') as bed_file:
        for sequence, chromosome, length, interval in far.entries():
            chromosome = chromosome.decode()
            sequence = sequence.decode()

            for start, kmer in enumerate(iter_kmers(sequence, k=K_MER_LENGTH)):
                try:
                    name = kmer_dict[kmer]
                except KeyError:
                    try:
                        kmer = dinopy.reverse_complement(kmer)
                        name = kmer_dict[kmer]
                    except KeyError:
                        continue

                end = start+K_MER_LENGTH

                bed_file.write(f'{chromosome}\t{start}\t{end}\t{name}\n')

make_bigwig(snakemake.input.genomic_fasta,
            snakemake.input.kmer_score_h5,
            snakemake.wildcards.motif_id,
            snakemake.params.score_column,
            snakemake.params.nlargest,
            snakemake.output.bed)
