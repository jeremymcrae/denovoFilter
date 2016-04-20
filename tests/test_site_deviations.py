'''
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import unittest
import math

from pandas import DataFrame, Series

from denovoFilter.site_deviations import site_strand_bias, test_sites, test_genes

class TestSiteDeviations(unittest.TestCase):
    
    def setUp(self):
        
        self.counts = DataFrame({'person_stable_id': ['a', 'b'],
            'chrom': ['1', '1'],
            'pos': [1, 2],
            'ref': ['A', 'G'],
            'alt': ['C', 'T'],
            'symbol': ['TEST1', 'TEST2'],
            'child_ref_F': [40, 15], 'child_ref_R': [15, 15],
            'child_alt_F': [20, 15], 'child_alt_R': [25, 15],
            'mother_ref_F': [40, 20], 'mother_ref_R': [15, 20],
            'mother_alt_F': [0, 0], 'mother_alt_R': [1, 0],
            'father_ref_F': [60, 30], 'father_ref_R': [30, 30],
            'father_alt_F': [0, 0], 'father_alt_R': [1, 1],
            'min_parent_alt': [1, 0]
            })
    
    def check_series(self, a, b):
        '''
        '''
        
        for x, y in zip(a, b):
            if math.isnan(x):
                self.assertTrue(math.isnan(x))
                self.assertTrue(math.isnan(y))
            else:
                self.assertAlmostEqual(x, y, 14)
    
    def test_site_strand_bias(self):
        ''' check that site_strand_bias works correctly
        '''
        
        site = {'ref_F': 5, 'ref_R': 30, 'alt_F': 10, 'alt_R': 10}
        self.assertAlmostEqual(site_strand_bias(site), 0.010024722592, places=11)
        
        site = {'ref_F': 30, 'ref_R': 30, 'alt_F': 10, 'alt_R': 10}
        self.assertEqual(site_strand_bias(site), 1.0)
        
        # zero counts give a p-value of 1.0
        site = {'ref_F': 0, 'ref_R': 0, 'alt_F': 0, 'alt_R': 0}
        self.assertEqual(site_strand_bias(site), 1.0)
        
        # check that values which would ordinarily give out of bounds errors
        # instead are converted to a p-value of 1.0. Some versions of scipy have
        # fixed this bug, and give a correct value, which we need to check too.
        site = {'ref_F': 1, 'ref_R': 2, 'alt_F': 9, 'alt_R': 84419233}
        self.assertIn(site_strand_bias(site), (1.0, 3.5536923140874242e-07))
    
    def test_test_sites(self):
        ''' check p-values from tests of strand and parental alt bias.
        '''
        
        expected = [[0.00061560815415820467, 1.0], [0.035457115371929658, 0.18307032892094907]]
        
        for x, y in zip(test_sites(self.counts), expected):
            self.check_series(x, y)
        
        # check when we mask some variants due to failing earlier variants
        expected = [[float('nan'), 1.0], [float('nan'), 0.18307032892094907]]
        for x, y in zip(test_sites(self.counts, pass_status=[False, True]), expected):
            self.check_series(x, y)
    
    def test_test_genes(self):
        ''' check p-values from test of parental alt bias within genes.
        '''
        
        # one of the sites fails the strand bias filter, so this gets dropped
        # for checking gene-based paternal alt biases
        sb, pa = test_sites(self.counts)
        expected = [float('nan'), 0.18307032892094907]
        
        self.check_series(test_genes(self.counts, sb), expected)
        
        # check when we mask some variants due to failing earlier variants
        expected = [float('nan'), float('nan')]
        self.check_series(test_genes(self.counts, sb, pass_status=[True, False]), expected)
