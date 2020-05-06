import cython
cimport cython

import numpy as np
cimport numpy as np


cdef enum NucleotideEncoding:

  A = 0 # 0b00
  C = 1 # 0b01
  G = 2 # 0b10
  T = 3 # 0b11

# K-mer length (code assumes all kmers will be of length 10)
# this is hardcoded because it's simpler and this library is designed for one thing only
cdef int K = 10


cpdef int encode_kmer (char* sequence, bint reverse_complement = 0) except -1:
    """
    Encodes a kmer (where K is a constant defined above)
    using bit representation in the `NucleotideEncoding` enum.

    Can also do a reverse complement encoding. It is simpler to do it here
    than to play with bits.

    :param sequence: sequence to encode (python bytes string)
    :param reverse_complement: boolean whether to encode reverse complement of sequence.
    """

    cdef int i = 0
    cdef int pos = 0

    cdef int letter_code = 0

    cdef int answer = 0

    for i in range(K):

        if not reverse_complement:
            pos = K-i-1
            if sequence[pos] == 'A':
                letter_code = NucleotideEncoding.A
            elif sequence[pos] == 'C':
                letter_code = NucleotideEncoding.C
            elif sequence[pos] == 'G':
                letter_code = NucleotideEncoding.G
            elif sequence[pos] == 'T':
                letter_code = NucleotideEncoding.T
            else:
                return -1
        else:
            pos = i
            if sequence[pos] == 'A':
                letter_code = NucleotideEncoding.T
            elif sequence[pos] == 'C':
                letter_code = NucleotideEncoding.G
            elif sequence[pos] == 'G':
                letter_code = NucleotideEncoding.C
            elif sequence[pos] == 'T':
                letter_code = NucleotideEncoding.A
            else:
                return -1

        answer = (answer << 2) + letter_code

    return answer


cpdef bytes decode_kmer (int encoded):
    """
    Decodes kmer encoded by `encode_kmer`
    :param encoded: the KMER encoded into integer
    """

    cdef bytes answer = b''

    cdef int mask
    cdef int value

    for i in range(K):
        # 0b11 << 2-n
        mask = 3 << (2*i)

        value = encoded & mask
        value = value >> (2*i)

        if value == NucleotideEncoding.A:
            answer += b'A'
        elif value == NucleotideEncoding.C:
            answer += b'C'
        elif value == NucleotideEncoding.G:
            answer += b'G'
        else:
            answer += b'T'

    return answer

cpdef int advance_kmer (int encoded, char next_char) except -1:
    """
        Shifts encoded kmer to the right by one nucleotide.
        E.g.
        ACGTACGTAC kmer shifted by T will become CGTACGTACT

        :param encoded: encoded kmer
        :param char: next character to add
    """

    cdef int ans = (encoded >> 2)

    cdef int letter_code;

    if next_char == 'A':
      letter_code = NucleotideEncoding.A
    elif next_char == 'C':
      letter_code = NucleotideEncoding.C
    elif next_char == 'G':
      letter_code = NucleotideEncoding.G
    elif next_char == 'T':
      letter_code = NucleotideEncoding.T
    else:
      return -1

    letter_code = (letter_code << 2*(K-1))
    return ans + letter_code


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)
cpdef np.ndarray[double, ndim=1] consume_sequence(
    np.ndarray[double, ndim=1] lookup,
    bytes sequence):
    """
        Consumes a DNA sequence provided in `sequence`
        and returns a numpy (float) array which contains the values
        from `lookup` corresponding to encoded kmers.

        Kmer positions are assumed to be centered:

        ACGTACGTAC
            ^ pos


        Lookup is expected to be a numpy array of dimension 4^K
        where each entry in position i corresponds to a value for
        kmer which is encoded to i.

        :params lookup: kmer value lookup array
        :params sequence: sequence to consume
    """

    # upercase to avoid problems
    sequence = sequence.upper()

    # Center sequences at Kmer
    cdef int i = -(K / 2)
    cdef int length = len(sequence)

    assert lookup.shape[0] == 4**K

    cdef np.ndarray[double, ndim=1] ans = np.full(length, np.nan)

    cdef bytes kmer = b''
    cdef int kmer_encoded = -1

    for c in sequence:
        i += 1
        # if c == b'N':
        if c == 78:
            kmer = b''
            kmer_encoded = -1
        elif kmer_encoded == -1:
            # Hacky but I cannot figure out the types now
            kmer += chr(c).encode()
            if len(kmer) == K:
                kmer_encoded = encode_kmer(kmer)
                ans[i] = lookup[kmer_encoded] # lookup
            else:
                continue
        else:
            kmer_encoded = advance_kmer(kmer_encoded, c)
            ans[i] = lookup[kmer_encoded]  # lookup

    return ans

cpdef np.ndarray[double, ndim=1] create_lookup_array(
    index,
    double[:] values,
):
    """
        Takes an index array and values array and creates appropriate
        lookup array

        Assumes kmers are in the index of series.
        Assumes that values for reverse-complement are the same and
        series contains only +, or - strand, but not both.

        Missing values will have a np.nan in the lookup

        :param index: kmers to which values correspond to
        :param values: values to be encoded

    """

    assert len(index) == len(values)

    cdef int n = len(index)

    cdef np.ndarray[double, ndim=1] ans = np.full(4**K, np.nan)

    cdef int i;

    cdef int encoding;
    cdef int encodeing_rc;

    # There's no good way to work with byte arrays anyway,
    # So we might as well take the performance hit:
    for i, sequence in enumerate(index):
        encoding = encode_kmer(sequence)
        encoding_rc = encode_kmer(sequence, reverse_complement=1)

        ans[encoding] = values[i]
        ans[encoding_rc] = values[i]

    return ans
