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

from denovoFilter.filter_denovogear_sites import filter_denovogear_sites

class TestFilterDenovogearSites(unittest.TestCase):
    
    def setUp(self):
        
        self.variants = DataFrame({'person_stable_id': ['a', 'b'],
            'chrom': ['1', '1'],
            'pos': [1, 1],
            'ref': ['A', 'A'],
            'alt': ['C', 'T'],
            'symbol': ['TEST1', 'TEST1'],
            })
        
        self.variants['dp4_child'] = ['20,15,30,25', '20,15,30,25']
        self.variants['dp4_mother'] = ['30,30,0,1', '30,30,0,1']
        self.variants['dp4_father'] = ['30,30,0,1', '30,30,0,1']
    
    def test_filter_denovogear_sites(self):
        ''' check that filter_denovogear_sites() works correctly
        '''
        
        initial = [True, True]
        
        status = filter_denovogear_sites(self.variants, initial)
        self.assertTrue(all(status == Series([True, True])))
    
    def test_filter_denovogear_sites_strand_bias(self):
        ''' check that filter_denovogear_sites() identifies strand biased sites
        '''
        
        initial = [True, True]
        
        # define a site with extreme strand bias
        self.variants['pos'] = [1, 2]
        self.variants['dp4_child'] = ['40,15,20,25', '20,15,30,25']
        self.variants['dp4_mother'] = ['40,15,0,1', '20,20,0,1']
        self.variants['dp4_father'] = ['60,30,0,1', '30,30,0,1']
        
        status = filter_denovogear_sites(self.variants, initial)
        self.assertTrue(all(status == Series([False, True])))
        
        # but indels don't use the same strand bias filter
        self.variants['ref'] = ['AG', 'A']
        status = filter_denovogear_sites(self.variants, initial)
        self.assertTrue(all(status == Series([True, True])))
    
    def test_filter_denovogear_sites_parental_bias(self):
        ''' check that filter_denovogear_sites() identifies parentally biased sites
        '''
        
        initial = [True, True]
        
        # define a site with extreme parental bias
        self.variants['pos'] = [1, 2]
        self.variants['dp4_child'] = ['20,15,20,25', '20,15,30,25']
        self.variants['dp4_mother'] = ['10,10,0,2', '20,20,0,1']
        self.variants['dp4_father'] = ['10,10,0,1', '30,30,0,0']
        
        status = filter_denovogear_sites(self.variants, initial)
        self.assertTrue(all(status == Series([False, True])))
        
        print('\nzzz\n')
        
        # define a site which only fails one of the filtering criteria, an
        # extreme parental alt bias at a site.
        self.variants['pos'] = [1, 2]
        self.variants['dp4_child'] = ['20,15,20,25', '20,15,30,25']
        self.variants['dp4_mother'] = ['8,8,0,3', '100,100,0,1']
        self.variants['dp4_father'] = ['10,10,0,0', '100,100,0,0']
        
        status = filter_denovogear_sites(self.variants, initial)
        print('STATUS: {}'.format(list(status)))
        self.assertTrue(all(status == Series([True, True])))
    
