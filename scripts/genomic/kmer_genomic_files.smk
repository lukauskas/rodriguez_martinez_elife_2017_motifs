# Import this from main snakefile

rule make_bigwig:
    input:
        genomic_fasta = rules.genomic_sequence_download.output.fa_gz,
        kmer_score_h5 = rules.combine_kmers_to_dataframe.output.h5,
    output:
        bw=DIR_OUTPUT / 'bw' / '{genome_version}' / '{motif_id}.ranks.bw'
    params:
        score_column = 'z_score',
        nlargest = 200 # show only top 200 kmers
    threads:
        8  # large memory usage
    script:
        'kmer_bigwig.py'

rule make_bed:
    input:
        genomic_fasta = rules.genomic_sequence_download.output.fa_gz,
        kmer_score_h5 = rules.combine_kmers_to_dataframe.output.h5,
    output:
        bed=temp(DIR_INTERIM / 'bed' / '{genome_version}' / '{motif_id}.kmers.bed')
    params:
        score_column = 'z_score',
        nlargest = 50
    script:
        'kmer_bed.py'

rule gzip_bed:
    input:
        bed = rules.make_bed.output.bed
    output:
        bed_gz = DIR_OUTPUT / 'bed' / '{genome_version}' / '{motif_id}.kmers.bed.gz'
    shell:
        """
        cat {input.bed} | gzip > {output.bed_gz}
        """
