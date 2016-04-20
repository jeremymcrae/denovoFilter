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

from pandas import Series

from denovoFilter.site_deviations import site_strand_bias

class TestSiteDeviations(unittest.TestCase):
    
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
    
    
