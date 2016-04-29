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

from denovoFilter.check_independence import check_independence, \
    person_recurrence, family_recurrence

class TestCheckIndependence(unittest.TestCase):
    
    def setUp(self):
        
        self.variants = DataFrame({'person_stable_id': ['a', 'b'],
            'chrom': ['1', '1'],
            'pos': [1, 2],
            'ref': ['A', 'G'],
            'alt': ['C', 'T'],
            'symbol': ['TEST1', 'TEST2'],
            'consequence': ['missense_variant', 'splice_donor_variant']
            })
    
    def test_person_recurrence(self):
        ''' test that person_recurrence() works correctly
        '''
        
        # in two variants for different samples, neither can be recurrent
        status = person_recurrence(self.variants)
        self.assertTrue(all(status == Series([False, False])))
        
        # if the two variants are in the same individual and gene, then one is
        # recurrent, the one with the les severe consequence.
        self.variants['person_stable_id'] = ['a', 'a']
        self.variants['symbol'] = ['TEST1', 'TEST1']
        status = person_recurrence(self.variants)
        self.assertTrue(all(status == Series([True, False])))
    
    def test_family_recurrence(self):
        ''' test that family_recurrence() works correctly
        '''
        
        family_ids = {'a': 'fam1', 'b': 'fam1'}
        status = family_recurrence(self.variants, family_ids)
        self.assertTrue(all(status == Series([False, False])))
        
        # set up some variants that are at the same position
        self.variants['symbol'] = ['TEST1', 'TEST1']
        self.variants['chrom'] = ['1', '1']
        self.variants['pos'] = [100, 100]
        self.variants['ref'] = ['C', 'C']
        self.variants['alt'] = ['T', 'T']
        
        # if the sites are at the same position, and in the ssame family, then
        # we'll fail all but the first occurrence
        status = family_recurrence(self.variants, family_ids)
        self.assertTrue(all(status == Series([False, True])))
        
        # if the sites are not in the same family, then they can't be recurrent
        family_ids = {'a': 'fam1', 'b': 'fam2'}
        status = family_recurrence(self.variants, family_ids)
        self.assertTrue(all(status == Series([False, False])))
    
    def test_check_independence(self):
        ''' check that check_independence works correctly
        '''
        
        family_ids = {'a': 'fam1', 'b': 'fam2'}
        status = check_independence(self.variants, family_ids)
        self.assertTrue(all(status == Series([True, True])))
        
        # now check sites with have intra-individual gene recurrence, the status
        # should be inverted with respect to the person_recurrence and
        # family_recurrence functions
        self.variants['person_stable_id'] = ['a', 'a']
        self.variants['symbol'] = ['TEST1', 'TEST1']
        status = check_independence(self.variants, family_ids)
        self.assertTrue(all(status == Series([False, True])))
