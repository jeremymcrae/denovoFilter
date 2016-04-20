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
import tempfile
import json

from pandas import DataFrame

from denovoFilter.change_last_base_sites import change_conserved_last_base_consequence
from tests.compare_dataframes import CompareTables

class TestChangeLastBaseSites(CompareTables):
    
    def test_change_last_base_sites(self):
        ''' check conversion of consequences for sites at the last bae of exons
        '''
        
        temp = tempfile.NamedTemporaryFile(mode='w')
        json.dump([['1', 10], ['2', 10]], temp)
        temp.flush()
        
        # define a table of variants. The first is a variant where the
        # consequence needs to be changed, the second isn't at an appropriate
        # location, and the third, although at an appropriate location, isn't
        # for a SNV, so does not need to be changed
        variants = DataFrame({'person_stable_id': ['a', 'b', 'c'],
            'chrom': ['1', '2', '2'],
            'pos': [10, 1, 10],
            'ref': ['A', 'G', 'GG'],
            'alt': ['C', 'T', 'G'],
            'consequence': ['synonymous_variant', 'synonymous_variant', 'frameshift_variant'],
            })
        
        # define an expected output, aside from a modified consequence
        expected = variants.copy()
        expected['consequence'] = ['conserved_exon_terminus_variant', 'synonymous_variant', 'frameshift_variant']
        
        self.compare_tables(change_conserved_last_base_consequence(variants, temp.name), expected)
