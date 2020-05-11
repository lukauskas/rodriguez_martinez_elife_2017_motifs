# Rodriguez-Martinez, eLife 2017 bZIP family motif kmer bigwigs

Some further post-processing of data from Rodriguez-Martinez et al. eLife 2017;6:e19272

What these scripts do, briefly, is to parse the kmer tracks from their supplementary data,
and create a bed file of occurences of top 50 motifs.

Alternatively, if the appropriate line is uncommented, the
script can generate a bigwig of top 200 motifs:

The resulting bigwigs show signal values of the level 0-200 where the scores
indicate the kmer ranks in the following fashion:

1. 0 - all kmers up to 200 best kmers
2. 1-200 top 200 kmers ranked in increasing score

## Usage

1. Download kmer files from (Rodríguez-Martínez, 2007) publication.
   These are available at the [Ansari Lab webpage](https://ansarilab.biochem.wisc.edu/computation.html).
   Only 10bp motifs are needed, and the txt files should be gzipped.
   The data should be put into `data/rodriguez_matrinez_elife_10bp_data/`.
   Check the checksums in the corresponding `checksums.txt`
2. Create a new virtual environment (python 3)
3. navigate to `scripts` folder and install requirements.txt `pip install -r requirements.txt`
4. install the kmer counter helper library: `pip install -e .`
5. run the snakemake pipeline `snakemake -j 8`
6. the output is in `output`. Note that each bigwig is about 1Gb.

## References

* Rodríguez-Martínez, J., Reinke, A., Bhimsaria, D., Keating, A., Ansari, A. (2017). Combinatorial bZIP dimers display complex DNA-binding specificity landscapes eLife  6(), e19272. https://dx.doi.org/10.7554/elife.19272
