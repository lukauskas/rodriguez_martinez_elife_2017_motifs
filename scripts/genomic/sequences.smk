# Import this from main snakefile

rule genomic_sequence_download:
    output:
        fa_gz=DIR_INTERIM / 'genomes' / '{genome_version}.sequence.fa.gz'
    shell:
        '''
        wget --quiet https://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.genome_version}/bigZips/{wildcards.genome_version}.fa.gz -O {output.fa_gz}
        '''
