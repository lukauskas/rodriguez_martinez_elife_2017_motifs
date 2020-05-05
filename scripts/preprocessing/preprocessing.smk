# Expects some global variables, import this from main snakefile

rule parse_kmers:
    """
        Parses kmers from the format provided
    """
    input:
        kmers_txt = DIR_KMER_INPUT / '{kmer_name}.10mer.txt.gz'
    output:
        kmers_filtered_txt = temp(DIR_INTERIM / '{kmer_name}.h5')
    run:
        import pandas as pd

        df = pd.read_csv(input.kmers_txt, sep='\t', header=None,
                         names=['kmer', 'kmer_rc', 'z_score', 'z_score_rank'])

        df = df.set_index('kmer')
        df.index.name = 'kmer'

        with pd.HDFStore(output.kmers_filtered_txt, 'w') as store:
            for col in ['z_score', 'z_score_rank']:
                series = df[col]
                series.name = wildcards.kmer_name
                store[col] = series


rule combine_kmers_to_dataframe:
    input:
        kmer_files =[ DIR_INTERIM / f'{kmer_name}.h5' for kmer_name in KMERS ]
    output:
        h5 = DIR_OUTPUT / 'kmer_z_score_dataframe.h5'
    threads: workflow.cores # Uses one CPU but a lot of memory, better not use any more
    run:
        import pandas as pd

        with pd.HDFStore(output.h5, 'w', complevel=1) as store:
            for col in ['z_score', 'z_score_rank']:
                ans = []

                for filename in input.kmer_files:
                    ans.append(pd.read_hdf(filename, col))

                ans = pd.concat(ans, axis=1)
                ans.sort_index(axis=0, inplace=True)
                ans.sort_index(axis=1, inplace=True)

                ans.index.name = 'kmer'

                store.put(col, ans, format='table')
