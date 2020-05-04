rule top_sequences_as_fasta:
    input:
        kmer_data = DIR_INTERIM / '{kmer_name}.h5'
    output:
        fa = temp(DIR_INTERIM / 'fasta' / '{kmer_name}.best_sequences.fa')
    params:
        column_to_use='z_score',
        n = config['N_KMERS_FOR_MEME']
    run:
        import pandas as pd
        n = int(params.n)
        series = pd.read_hdf(input.kmer_data, params.column_to_use)

        series = series.nlargest(n)
        series = series.sort_values(ascending=False)

        with open(output.fa, 'w') as file_:
            for i, kmer in enumerate(series.index, start=1):
                file_.write(f'>{i}\n{kmer}\n')

rule top_sequences_as_fasta_gz:
    input:
        fa=rules.top_sequences_as_fasta.output.fa
    output:
        fa=DIR_OUTPUT / 'fasta' / '{kmer_name}.best_sequences.fa.gz'
    shell:
        'cat {input.fa} | gzip > {output.fa}'
