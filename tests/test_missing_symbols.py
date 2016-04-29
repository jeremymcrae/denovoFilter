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

from pandas import DataFrame

from denovoFilter.missing_symbols import fix_missing_gene_symbols, open_url, \
    rate_limit_requests, get_gene_id

class TestMissingSymbols(unittest.TestCase):
    
    def setUp(self):
        self.variants = DataFrame({'person_stable_id': ['a', 'b'],
            'chrom': ['6', '2'],
            'pos': [157528051, 129119889],
            'ref': ['C', 'G'],
            'alt': ['T', 'A'],
            'symbol': ['', ''],
            'consequence': ['stop_gained', 'intergenic_variant'],
            })
    
    def test_fix_missing_gene_symbols(self):
        ''' check that get_most_severe works correctly
        '''
        
        symbols = fix_missing_gene_symbols(self.variants)
        self.assertEqual(list(symbols), ['ARID1B', 'fake_symbol.2_129119889'])
    
    def test_open_url(self):
        ''' check that open_url works correctly
        '''
        
        headers = {'content-type': 'application/json'}
        
        url = 'http://grch37.rest.ensembl.org/overlap/region/human/6:157528051-157528051?feature=gene'
        response, status_code, headers = open_url(url, headers)
        self.assertIn('external_name', response)
        self.assertEqual(status_code, 200)
        self.assertIn('content-length', headers)
        
        # a nonsense URL gives an error
        with self.assertRaises(ValueError):
            open_url('example', headers)
        
        # # The following code doesn't seem to work on travis-ci, so I have
        # # commented it out. The HTTPError part is later captured when testing
        # # the get_gene_id function, by passing an invalid ensembl REST URL.
        # response, status_code, headers = open_url('http://httpbin.org/status/500', headers)
        # self.assertIn(status_code, [500, 411])
    
    def test_rate_limit_requests(self):
        ''' test that rate_limit_requests() works correctly
        '''
        
        # run this once, otherwise the first time around it won't pause, since
        # the starting time was when the function was imported.
        rate_limit_requests()
        
        initial = time.time()
        rate_limit_requests()
        self.assertTrue(time.time() - initial > 0.0667)
        
        initial = time.time()
        rate_limit_requests(rate_limit=0.1)
        self.assertTrue(time.time() - initial > 0.1)
    
    def test_get_gene_id(self):
        ''' test that get_gene_id() works correctly
        '''
        
        self.assertEqual(get_gene_id('6', 157528051, 157528051), 'ARID1B')
        
        with self.assertRaises(ValueError):
            get_gene_id('NOTCHROM', 'AAAA', 'BBBB')
        
        with self.assertRaises(ValueError):
            get_gene_id('6', 157528051, 157528051, attempts=6)
        
    
