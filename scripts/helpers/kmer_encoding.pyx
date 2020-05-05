import numpy as np
cimport numpy as np

cdef enum NucleotideEncoding:

  A = 0 # 0b00
  C = 1 # 0b01
  G = 2 # 0b10
  T = 3 # 0b11

# K-mer length (code assumes all kmers will be of length 10)
# Everything will break if they're not..

cdef int K = 10


cpdef int encode_kmer (char* sequence) except -1:

  cdef int i = 0

  cdef int letter_code = 0

  cdef int answer = 0

  for i in range(K):

    if sequence[K-i-1] == 'A':
      letter_code = NucleotideEncoding.A
    elif sequence[K-i-1] == 'C':
      letter_code = NucleotideEncoding.C
    elif sequence[K-i-1] == 'G':
      letter_code = NucleotideEncoding.G
    elif sequence[K-i-1] == 'T':
      letter_code = NucleotideEncoding.T
    else:
      return -1


    answer = (answer << 2) + letter_code

  return answer


cpdef bytes decode_kmer (int encoded):

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


# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)
cpdef np.ndarray[double, ndim=1] consume_sequence(
    np.ndarray[double, ndim=1] lookup,
    bytes sequence):

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
