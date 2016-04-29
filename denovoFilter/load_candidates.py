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

import pandas

def load_candidates(candidates_path):
    """ load the candidate dataset
    
    Args:
        candidates_path: path to candidate DNMs
    
    Returns:
        pandas dataframe of candidate de novo sites
    """
    
    candidates = pandas.read_table(candidates_path, na_filter=False)
    
    # the missing indels don't have some columns that are present in the
    # denovogear input, use mock columns so that later processing works smoothly.
    if 'in_child_vcf' not in candidates.columns:
        candidates["in_child_vcf"] = 1
        candidates["in_mother_vcf"] = 0
        candidates["in_father_vcf"] = 0
        candidates["pp_dnm"] = None
        
    return candidates
