# Import this from main snakefile

rule make_bigwig:
    input:
        genomic_fasta = rules.genomic_sequence_download.output.fa_gz,
        kmer_score_h5 = rules.combine_kmers_to_dataframe.output.h5,
    output:
        bw=DIR_OUTPUT / 'bw' / '{genome_version}' / '{motif_id}.ranks.bw'
    params:
        score_column = 'z_score',
        nlargest = 1000 # show only top 1000 kmers
    script:
        'kmer_bigwig.py'
