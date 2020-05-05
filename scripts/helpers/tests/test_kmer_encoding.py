from helpers.kmers_to_bigwig import *
from nose2.tools import params
from hypothesis import given
from hypothesis.strategies import text

import unittest
import numpy as np
from numpy.testing import assert_array_equal

def window(seq, n=2):
    # Based on sliding window iterator in Stackoverflow (originally from python docs)
    result = seq[:n]
    if len(result) == n:
        yield result
    for elem in seq[n:]:
        result = result[1:] + chr(elem).encode()
        yield result


class TestKmerEncoding(unittest.TestCase):

    @given(text(alphabet='ACGT', min_size=10, max_size=10))
    def test_roundtrip(self, sequence):
        sequence = sequence.encode() # bytes

        encoded = encode_kmer(sequence)
        decoded = decode_kmer(encoded)

        self.assertEquals(sequence, decoded)

    @given(text(alphabet='ACGT', min_size=11, max_size=11))
    def test_advance(self, sequence):
        sequence = sequence.encode() # bytes

        sequence_start = sequence[:10]
        next_char = sequence[10]
        sequence_next = sequence[1:]

        encoded_start = encode_kmer(sequence_start)
        encoded_next = advance_kmer(encoded_start, next_char)
        decoded_next = decode_kmer(encoded_next)

        self.assertEquals(sequence_next, decoded_next)


    def test_consume_sequence_manual(self):

        self_reporting_lookup = np.arange(4**10, dtype=np.float)
        some_sequence = b'AAAAAAAAAACGT'

        kmers = [None] * 5 + [b'AAAAAAAAAA', b'AAAAAAAAAC', b'AAAAAAAACG', b'AAAAAAACGT'] + [None] * 4
        expected_values = np.array([np.nan if x is None else encode_kmer(x) for x in kmers])

        actual_values = consume_sequence(self_reporting_lookup, some_sequence)
        assert_array_equal(expected_values, actual_values)

    @given(text(alphabet='ACGTacgtN', min_size=11, max_size=100))
    def test_consume_sequence_automatic(self, sequence):

        # To bytes
        sequence = sequence.encode()

        self_reporting_lookup = np.arange(4**10, dtype=np.float)
        # To make sure we're using lookup not just index
        offset = 13
        self_reporting_lookup += offset

        expected_ans = np.full(len(sequence), np.nan)

        for i, kmer in enumerate(window(sequence.upper(), n=10)):
            if b'N' in kmer:
                continue

            encoded = encode_kmer(kmer)
            expected_ans[i+5] = encoded + offset

        actual_ans = consume_sequence(self_reporting_lookup, sequence)
        assert_array_equal(expected_ans, actual_ans)
