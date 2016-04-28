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

from denovoFilter.most_severe import get_most_severe

class TestMostSevere(unittest.TestCase):
    
    def test_get_most_severe(self):
        ''' check that get_most_severe works correctly
        '''
        
        cq = ['missense_variant', 'protein_altering_variant',
            'splice_region_variant', 'incomplete_terminal_codon_variant']
        self.assertEqual(get_most_severe(cq), 'missense_variant')
        
        cq = ['stop_lost', 'start_lost', 'transcript_amplification',
            'conserved_exon_terminus_variant']
        self.assertEqual(get_most_severe(cq), 'stop_lost')
        
        # an empty list raises an error
        with self.assertRaises(TypeError):
            get_most_severe([])
        
        
        
    
