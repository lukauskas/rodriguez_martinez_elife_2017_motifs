from kmer_helpers.kmer_encoding import *
from nose2.tools import params
from hypothesis import given, settings
from hypothesis.strategies import text

import unittest
import numpy as np
from numpy.testing import assert_array_equal

import dinopy
import itertools

N_HYPOTHESIS_TESTS = 5000


def window(seq, n=2):
    # Based on sliding window iterator in Stackoverflow (originally from python docs)
    result = seq[:n]
    if len(result) == n:
        yield result
    for elem in seq[n:]:
        result = result[1:] + chr(elem).encode()
        yield result


class TestKmerEncoding(unittest.TestCase):

    @settings(max_examples=N_HYPOTHESIS_TESTS)
    @given(text(alphabet='ACGT', min_size=10, max_size=10))
    def test_roundtrip(self, sequence):
        sequence = sequence.encode() # bytes

        encoded = encode_kmer(sequence)
        decoded = decode_kmer(encoded)

        self.assertEquals(sequence, decoded)

    @settings(max_examples=N_HYPOTHESIS_TESTS)
    @given(text(alphabet='ACGT', min_size=10, max_size=10))
    def test_reverse_complement(self, sequence):

        sequence = sequence.encode() # bytes
        complement = dinopy.processors.reverse_complement(sequence)

        encoded = encode_kmer(sequence, reverse_complement=True)
        encoded_complement = encode_kmer(complement, reverse_complement=False)

        self.assertEquals(encoded, encoded_complement)

    @settings(max_examples=N_HYPOTHESIS_TESTS)
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

    @settings(max_examples=N_HYPOTHESIS_TESTS)
    @given(text(alphabet='ACGTacgtN', min_size=1, max_size=200))
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

    def test_creating_lookup_array_smallscale(self):

        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        index = np.array([b'AAAAAAAAAA',
                          b'AAAAAAAAAC',
                          b'AAAAAAAAAT',
                          b'CCCCCCCCCC',
                          b'ACACACACAT',
                          b'CCGCTCACCT'], dtype='a10')

        expected_array = np.full(4**10, np.nan)

        for ix, v in zip(index, values):
            expected_array[encode_kmer(ix)] = v
            expected_array[encode_kmer(ix, reverse_complement=True)] = v

        actual_array = create_lookup_array(index, values)

        assert_array_equal(expected_array, actual_array)


    def test_creating_lookup_array_largescale(self):

        seen = set()
        index = []
        values = []

        expected_array = np.full(4**10, np.nan)

        for i, kmer in enumerate(itertools.product('ACGT', repeat=10)):
            kmer = ''.join(kmer).encode()
            rc_kmer = dinopy.processors.reverse_complement(kmer)

            if kmer in seen or rc_kmer in seen:
                continue

            seen.add(kmer)
            seen.add(rc_kmer)

            value = 42.0 + i
            index.append(kmer)
            values.append(value)

            expected_array[encode_kmer(kmer)] = value
            expected_array[encode_kmer(kmer, reverse_complement=True)] = value

        index = np.array(index)
        values = np.array(values)

        actual_array = create_lookup_array(index, values)
        assert_array_equal(expected_array, actual_array)
