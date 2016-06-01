"""
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
"""

import numpy
import pandas
import unittest

from denovoFilter.min_depth import min_depth

class TestMinDepth(unittest.TestCase):
    
    def test_min_depth(self):
        ''' test estimating the minumum allowed depth
        '''
        
        self.assertEqual(min_depth(100, 0.01), 2)
        self.assertEqual(min_depth(150, 0.03), 7)
        
        # and try modulating the certainty threshold
        self.assertEqual(min_depth(150, 0.03, threshold=0.5), 3)
        self.assertEqual(min_depth(150, 0.03, threshold=0.999), 9)
    
    def test_min_depth_data_types(self):
        ''' test that min_depth works on numerous data types
        '''
        
        self.assertEqual(min_depth([75, 150], 0.03), 5)
        self.assertEqual(min_depth((75, 150), 0.03), 5)
        self.assertEqual(min_depth(numpy.array([75, 150]), 0.03), 5)
        self.assertEqual(min_depth(pandas.Series([75, 150]), 0.03), 5)
    
    def test_min_depth_errors(self):
        ''' test that min depth raises appropriate errors
        '''
        
        with self.assertRaises(AssertionError):
            # raise errors with non-standard lengths of depths
            min_depth([75, 150, 50], 0.03)
            min_depth([75], 0.03)
            min_depth([], 0.03)
            
            # raise an error with an unobtainable certainty threshold
            min_depth(75, 0.01, 1.0)
