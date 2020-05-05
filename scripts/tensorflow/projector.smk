# Expects some global variables, import this from main snakefile

rule projector_data:
    input:
        h5=rules.combine_kmers_to_dataframe.output.h5
    output:
        vectors = DIR_OUTPUT / 'tensorflow_projector' / 'vectors.tsv',
        metadata = DIR_OUTPUT / 'tensorflow_projector' / 'metadata.tsv'
    params:
        col = 'z_score_rank'
    run:
        import pandas as pd

        df = pd.read_hdf(input.h5, params.col)
        df.to_csv(output.vectors, sep='\t', index=False, header=False)
        df.to_csv(output.metadata, sep='\t', index=True, header=True)
