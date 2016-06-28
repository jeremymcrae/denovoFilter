'''
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the 'Software'), to deal in
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
import time

from pandas import DataFrame, Series

from denovoFilter.preliminary_filtering import fix_maf, preliminary_filtering, \
    check_coding

class TestPreliminaryFiltering(unittest.TestCase):
    
    def setUp(self):
        self.variants = DataFrame({'person_stable_id': ['a', 'b'],
            'max_af': ['0.001', 'missing'],
            'in_child_vcf': [1, 1],
            'in_father_vcf': [0, 0],
            'in_mother_vcf': [0, 0],
            })
    
    def test_fix_maf(self):
        ''' test that fix_maf works correctly
        '''
        
        values = Series(['0.001', '0.2,0.3', '', '.', 'missing', None])
        
        self.assertTrue(all(fix_maf(values) == Series([0.001, 0.2, 0, 0, 0, 0])))
    
    def test_preliminary_filtering(self):
        ''' test that preliminary_filtering works correctly.
        '''
        
        # for a clean table, all variants should pass
        status = preliminary_filtering(self.variants)
        self.assertTrue(all(status == Series([True, True])))
        
        # when we define a set of samples that fail, their variants will not pass
        status = preliminary_filtering(self.variants, sample_fails=['b'])
        self.assertTrue(all(status == Series([True, False])))
        
        # when we adjust the MAF threshold, candidates above will faik
        status = preliminary_filtering(self.variants, maf_cutoff=0.0001)
        self.assertTrue(all(status == Series([False, True])))
        
        # if we set a parent to have been called, that site will fail
        self.variants['in_father_vcf'] = [1, 0]
        status = preliminary_filtering(self.variants)
        self.assertTrue(all(status == Series([False, True])))
    
    def test_check_coding(self):
        ''' check that check_coding works correctly
        '''
        
        # checking whether a site is coding or not works correctly
        self.variants['consequence'] = ['synonymous_variant', 'intergenic_variant']
        status = check_coding(self.variants)
        self.assertTrue(all(status == Series([True, False])))
        
        # make sure we can easily use a different column name
        self.variants = self.variants.drop('consequence', axis=1)
        self.variants['cq'] = ['synonymous_variant', 'intergenic_variant']
        status = check_coding(self.variants, cq_name='cq')
        self.assertTrue(all(status == Series([True, False])))
        
        # raise an error if we try a missing column.
        with self.assertRaises(KeyError):
            check_coding(self.variants, cq_name='UNKNOWN')
    
