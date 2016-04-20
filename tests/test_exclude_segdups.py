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

from denovoFilter.exclude_segdups import check_segdups

class TestExcludeSegdups(unittest.TestCase):
    
    def test_exclude_segdups(self):
        ''' check that counting alleles from DP4 entries works correctly
        '''
        
        # define variants that lie around the boundaries of a segdup region
        variants = DataFrame({'person_stable_id': ['a', 'a', 'a', 'a', 'a', 'a'],
            'chrom': ['1', '1', '1', '1', '1', '1'],
            'pos': [1379893, 1379895, 1379894, 1384309, 1384310, 1384311],
            'ref': ['A', 'G', 'A', 'G', 'A', 'G'],
            'alt': ['C', 'T', 'C', 'T', 'C', 'T'],
            'symbol': ['TEST1', 'TEST1', 'TEST1', 'TEST1', 'TEST1', 'TEST1'],
            })
        
        expected = [True, False, False, False, False, True]
        
        self.assertEqual(check_segdups(variants), expected)
    
