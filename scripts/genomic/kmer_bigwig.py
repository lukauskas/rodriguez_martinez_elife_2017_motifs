import pandas as pd
import numpy as np
import dinopy
import pyBigWig

from kmer_helpers.kmer_encoding import *
from tqdm import tqdm

K_MER_LENGTH = 10


def read_score_lookup_array(kmer_score_h5, motif_id, score_column, nlargest):
    """
        Reads kmer scores from hdf5
    """
    kmer_scores = pd.read_hdf(kmer_score_h5, score_column)[motif_id]
    kmer_scores = kmer_scores.nlargest(nlargest)

    kmer_scores = kmer_scores.rank(pct=True)
    # Change scale
    kmer_scores = kmer_scores * 1000

    kmer_score_lookup = create_lookup_array(
        kmer_scores.index.astype(bytes),
        kmer_scores.values)

    # Fillna with 0.0
    return np.nan_to_num(kmer_score_lookup, nan=0.0)

def make_bigwig(genomic_fasta,
                kmer_score_h5,
                motif_id,
                score_column,
                nlargest,
                output_file):

    print('Reading scores')
    kmer_score_lookup = read_score_lookup_array(kmer_score_h5, motif_id,
                                                score_column, nlargest)

    shp = dinopy.shape.Shape(K_MER_LENGTH)

    far = dinopy.FastaReader(genomic_fasta)


    header = []
    for __, chromosome, length, __ in far.entries():
        header.append((chromosome.decode(), length))

    with pyBigWig.open(output_file, 'w') as bw:
        bw.addHeader(header)

        for sequence, chromosome, length, interval in far.entries():
            chromosome = chromosome.decode()
            values = consume_sequence(kmer_score_lookup, sequence)
            # Cast to int to get more repeating values
            # values = np.asarray(values, dtype=int)


            pos = 0
            bw.addEntries(chromosome,
                          0,
                          values=values,
                          span=1,
                          step=1)




make_bigwig(snakemake.input.genomic_fasta,
            snakemake.input.kmer_score_h5,
            snakemake.wildcards.motif_id,
            snakemake.params.score_column,
            snakemake.params.nlargest,
            snakemake.output.bw)
