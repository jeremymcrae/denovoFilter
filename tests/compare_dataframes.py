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

import unittest
import math

import numpy

class CompareTables(unittest.TestCase):
    def compare_tables(self, first, second):
        """ check if two Dataframes contain the same information.
        
        Note that the dataframes don't have to have the same column order.
        
        Args:
            first: pandas DataFrame
            second: pandas DataFrame
        """
        
        # make sure the two tables have the same columns, and number of rows
        self.assertEqual(set(first.columns), set(second.columns))
        self.assertEqual(len(first), len(second))
        
        # make sure the values in the columns are identical between tables.
        # We have to run through all the values one by one, to make sure that
        # NaN values compare correctly.
        for column in first:
            for pos in range(len(first)):
                first_val = first[column][pos]
                second_val = second[column][pos]
                
                # prepare a suitable diagnostic error message if the values are
                #  not equal
                if type(first_val) == str or first_val is None or type(first_val) == numpy.bool_:
                    msg = "{} != {} at position {} in {}".format(first_val,
                        second_val, pos, column)
                else:
                    msg = "{:.20f} != {:.20f} at position {} in {}".format( \
                        first_val, second_val, pos, column)
                
                if type(first_val) in [float, numpy.float64] and math.isnan(first_val):
                    # if one is nan, check that the other is also nan, since
                    # nan entries cannot be chacked for eqaulity.
                    self.assertTrue(math.isnan(second_val), msg)
                elif type(first_val) in [float, numpy.float64]:
                    # Only check if floats are nearly equal to get around float
                    # precision issues.
                    self.assertAlmostEqual(first_val, second_val, msg=msg)
                else:
                    self.assertEqual(first_val, second_val, msg)
