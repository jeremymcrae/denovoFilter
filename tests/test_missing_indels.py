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

from denovoFilter.missing_indels import filter_missing_indels

class TestMissingIndels(unittest.TestCase):
    
    def setUp(self):
        
        self.variants = DataFrame({'person_stable_id': ['a', 'b'],
            'chrom': ['1', '1'],
            'pos': [1, 1],
            'ref': ['A', 'A'],
            'alt': ['C', 'TC'],
            'symbol': ['TEST1', 'TEST1'],
            })
        
        self.variants['dp4_child'] = ['20,15,30,25', '20,15,30,25']
        self.variants['dp4_mother'] = ['20,20,0,1', '20,20,0,1']
        self.variants['dp4_father'] = ['30,30,0,1', '30,30,0,1']
    
    def test_filter_missing_indels(self):
        ''' test that filter_missing_indels() passes the default variants
        '''
        
        status = filter_missing_indels(self.variants)
        self.assertTrue(all(status == Series([True, True])))
    
    def test_filter_missing_indels_child_depth(self):
        ''' test that filter_missing_indels() screens low child alt depth
        '''
        
        self.variants['dp4_child'] = ['20,15,30,25', '2,2,1,1']
        status = filter_missing_indels(self.variants)
        self.assertTrue(all(status == Series([True, False])))
    
    def test_filter_missing_indels_parental_alt(self):
        ''' test that filter_missing_indels() screens high parental alts
        '''
        
        self.variants['dp4_mother'] = ['20,20,0,1', '20,20,0,2']
        self.variants['dp4_father'] = ['30,30,0,1', '30,30,0,2']
        status = filter_missing_indels(self.variants)
        self.assertTrue(all(status == Series([True, False])))
    
    def test_filter_missing_indels_parental_depth(self):
        ''' test that filter_missing_indels() screens low parental depth
        '''
        
        self.variants['dp4_mother'] = ['20,20,0,1', '3,3,0,0']
        status = filter_missing_indels(self.variants)
        self.assertTrue(all(status == Series([True, False])))
    
    def test_filter_missing_indels_parental_proportion(self):
        ''' test that filter_missing_indels() screens high parental proportion
        '''
        
        self.variants['dp4_mother'] = ['20,20,0,1', '4,4,0,1']
        status = filter_missing_indels(self.variants)
        self.assertTrue(all(status == Series([True, False])))
    
    def test_filter_missing_indels_child_proportion(self):
        ''' test that filter_missing_indels() screens low child proportion
        '''
        
        self.variants['dp4_child'] = ['20,15,30,25', '30,30,7,7']
        status = filter_missing_indels(self.variants)
        self.assertTrue(all(status == Series([True, False])))
    
