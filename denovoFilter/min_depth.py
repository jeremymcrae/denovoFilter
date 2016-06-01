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

def min_depth(depth, error, threshold=0.98):
    '''  determine the maximum parental depth permitted from both parents
    
    We look at the minimum alternate depth across both parents. We need this to
    be low. This function determines how low we can set this and still capture
    the vast majority of true candidate de novo mutations.
    
    There are four posible scenarios:
        1) both parents have depths below or equal to the depth
        2) the first parent exceeds the depth but the second parent does not
        3) the second parent exceeds the depth but the first parent does not
        4) both parents exceed the depth.
    
    We only need to consider the fourth scenario. The probability that this
    happens (at a given depth) is 1 - prob(first parent exceeds) * prob(second
    parent exceeds). We repeatedly increment the depth threshold upwards until
    this probability is sufficiently high.
    
    Args:
        depth: depth of parental sequencing. Can either be a single value (as a
            summary value for both parents), or a list of two depths, one for
            each parent.
        error: site-specific error rate (e.g. 0.002)
        threshold: probability threshold that we are wanting to exceed. We look
            for a alt count where the probability exceeds this value. Assuming
            the value is 0.98, then 98% of cases will have a min depth less than
            the returned count.
    
    Returns:
        maximum permitted alternate depth across both parents.
    '''
    
    # convert int variables to a list, so we don't need code specific to ints
    try:
        depth = [int(depth), int(depth)]
    except TypeError:
        pass
    
    # raise an error if the length is not two.
    assert len(depth) == 2
    
    x = 0
    while True:
        product = 1
        for i in depth:
            product *= binom.sf(x, i, error)
        prob = 1 - product
        
        if prob > threshold:
            return(x)
        
        x += 1
