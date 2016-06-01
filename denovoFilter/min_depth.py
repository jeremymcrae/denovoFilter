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

from scipy.stats import binom

def min_depth(error, depth, threshold=0.98):
    '''  determine the maximum parental depth permitted from both parents
    
    We look at the minimum alternate depth across both parents. We need this to
    be low. This function determines how low we can set this and still capture
    the vast majority of true candidate de novo mutations.
    
    Args:
        depth: depth of parental sequencing. Can either be a single value (as a
            summary value for both parents), or a list of two depths, one for
            each parent.
        error: site-specific error rate (e.g. 0.002)
        threshold: probability threshold that we are wanting to exceed. We look for
            a alt count where the probability exceeds this value. Assuming the
            value is 0.98, then 98% of cases will have a min depth less than the
            returned count.
    
    Returns:
        maximum permitted alternate depth across both parents.
    '''
    
    x = 0
    while True:
        try:
            product = 1
            for i in depth:
                product *= binom.sf(x, i, error)
            prob = 1 - product
            
            # raise an error if the length is not two. This will not be a
            # problem for ints, which have been caught already.
            assert len(depth) == 2
        except TypeError:
            prob = 1 - (binom.sf(x, depth, error))**2
        
        if prob > threshold:
            return(x)
        
        x += 1
