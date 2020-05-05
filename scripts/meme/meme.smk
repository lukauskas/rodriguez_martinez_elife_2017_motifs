# Expects some global variables, import this from main snakefile

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

rule meme_motifs:
    input:
        fa=rules.top_sequences_as_fasta.output.fa
    output:
        main_directory = directory(DIR_OUTPUT / 'meme' / '{kmer_name}')
    params:
        # These are params from the publication
        # meme_params = '-dna -mod anr -nmotifs 10 -minw 8 -maxw 18 -time 7200 -maxsize 60000 -revcomp'
        # Slight adjustment below
        meme_params = '-dna -mod anr -nmotifs 10 -minw 8 -maxw 10 -time 7200 -revcomp -norand -minsites 10 -noendgaps -searchsize 1000'
    threads:
        2
    log: DIR_LOG / 'meme' / '{kmer_name}.log'
    shell:
        """
        meme {input.fa} {params.meme_params} -oc {output.main_directory} -p {threads} &> {log}
        """
