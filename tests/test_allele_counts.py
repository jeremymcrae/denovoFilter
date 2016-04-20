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

from pandas import DataFrame, Series

from denovoFilter.allele_counts import extract_alt_and_ref_counts, \
    get_recurrent_genes, get_allele_counts, get_depths_and_proportions
from tests.compare_dataframes import CompareTables

class TestAlleleCounts(CompareTables):
    
    def setUp(self):
        
        self.variants = DataFrame({'person_stable_id': ['a', 'b'],
            'chrom': ['1', '1'],
            'pos': [1, 2],
            'ref': ['A', 'G'],
            'alt': ['C', 'T'],
            'symbol': ['TEST1', 'TEST2'],
            })
        
        self.variants['dp4_child'] = ['20,15,30,25', '15,15,15,15']
        self.variants['dp4_mother'] = ['20,15,0,1', '20,20,0,0']
        self.variants['dp4_father'] = ['30,30,0,1', '30,30,0,1']
        
        self.counts = DataFrame({'person_stable_id': ['a', 'b'],
            'chrom': ['1', '1'],
            'pos': [1, 2],
            'ref': ['A', 'G'],
            'alt': ['C', 'T'],
            'symbol': ['TEST1', 'TEST2'],
            'child_ref_F': [20, 15], 'child_ref_R': [15, 15],
            'child_alt_F': [30, 15], 'child_alt_R': [25, 15],
            'mother_ref_F': [20, 20], 'mother_ref_R': [15, 20],
            'mother_alt_F': [0, 0], 'mother_alt_R': [1, 0],
            'father_ref_F': [30, 30], 'father_ref_R': [30, 30],
            'father_alt_F': [0, 0], 'father_alt_R': [1, 1],
            'min_parent_alt': [1, 0]
            })
    
    def test_extract_alt_and_ref_counts(self):
        ''' check that counting alleles from DP4 entries works correctly
        '''
        
        self.compare_tables(extract_alt_and_ref_counts(self.variants), self.counts)
    
    def test_get_recurrent_genes(self):
        ''' check that we identify the recurrently mutated genes correctly
        '''
        
        variants = DataFrame({'symbol': ['TEST1', 'TEST1', 'TEST2', 'TEST3']})
        
        self.assertEqual(get_recurrent_genes(variants), Series(['TEST1']))
    
    def test_get_allele_counts(self):
        ''' check that counting alleles works correctly
        '''
        
        # check the counts when the variants are aggregated for a gene
        expected = {'gene_alt': 3, 'gene_ref': 195}
        self.assertEqual(get_allele_counts(self.counts, gene=True), expected)
        
        # check the counts when the variants are aggregated for a site
        expected = {'ref_F': 135, 'ref_R': 125, 'alt_F': 45, 'alt_R': 43,
            'parent_alt': 3, 'parent_ref': 195}
        self.assertEqual(get_allele_counts(self.counts, gene=False), expected)
    
    def test_get_depths_and_proportions(self):
        ''' check that counting the depths and proportions works correctly
        '''
        
        expected = DataFrame({
            'child_depth': [90, 60], 'dad_depth': [61, 61], 'mom_depth': [36, 40],
            'child_alts': [55, 30], 'dad_alts': [1, 1], 'mom_alts': [1, 0],
            'child_prp': [55/90.0, 30/60.0], 'dad_prp': [1/61.0, 1/61.0],
            'mom_prp': [1/36.0, 0/40.0]}, index=self.counts.index)
        
        self.compare_tables(get_depths_and_proportions(self.counts), expected)
