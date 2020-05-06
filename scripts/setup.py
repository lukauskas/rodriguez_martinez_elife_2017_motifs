from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name='kmer helper scripts',
    ext_modules=cythonize("kmer_helpers/kmer_encoding.pyx",
                          include_path = [numpy.get_include()]),
    zip_safe=False,
    include_dirs=[numpy.get_include()]

)
