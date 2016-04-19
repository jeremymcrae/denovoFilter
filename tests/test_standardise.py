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

from pandas import DataFrame

from denovoFilter.standardise import standardise_columns
from tests.compare_dataframes import CompareTables

class TestStandardise(CompareTables):
    
    def setUp(self):
        self.initial = DataFrame({
            'person_stable_id': ['sample_01', 'sample_02'],
            'chrom': ['6', '6'],
            'pos': [1, 2],
            'ref': ['A', 'G'],
            'alt': ['C', 'T'],
            'symbol': ['TEST1', 'TEST2'],
            'consequence': ['missense_variant', 'stop_gained'],
            'max_af': [0.0001, 0],
            'pp_dnm': [1, 0.95],
            })
        
        self.expected = self.initial.copy()
    
    def test_standardise_columns(self):
        ''' test that standardising the columns works
        '''
        
        self.initial['child_ref_F'] = [23, 35]
        self.compare_tables(standardise_columns(self.initial), self.expected)
        
        # if the dataframes are different, expect an error
        self.expected['pp_dnm'] = [0.0, 0.0]
        with self.assertRaises(AssertionError):
            self.compare_tables(standardise_columns(self.initial), self.expected)
    
    def test_standardise_columns_with_pass(self):
        ''' test that standardising the columns works when a 'pass' column exists
        '''
        
        self.initial['pass'] = [True, False]
        self.initial['child_ref_F'] = [23, 35]
        
        self.expected['pass'] = [True, False]
        self.compare_tables(standardise_columns(self.initial), self.expected)
    
    def test_standardise_error(self):
        ''' test if a required column is missing, we raise an error
        '''
        
        self.initial = self.initial.drop('pp_dnm', axis=1)
        
        with self.assertRaises(KeyError):
            standardise_columns(self.initial)
