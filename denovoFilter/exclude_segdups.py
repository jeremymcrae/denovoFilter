"""
Copyright (c) 2015 Wellcome Trust Sanger Institute

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

from pkg_resources import resource_filename
import gzip

from intervaltree import IntervalTree

def load_segdups():
    """ load all the segdup regions
    
    Returns:
        dictionary of IntervalTree based segdup regions per chromosome
    """
    
    segdup_path = resource_filename(__name__, "data/segdup_regions.gz")
    
    segdups = {}
    with gzip.open(segdup_path, "rt") as handle:
        
        header = handle.readline()
        
        for line in handle:
            line = line.strip().split("\t")
            chrom = line[0].strip("chr")
            start = int(line[1])
            end = int(line[2]) + 1
            
            if chrom not in segdups:
                segdups[chrom] = IntervalTree()
            
            segdups[chrom].addi(start, end)
    
    return segdups

def check_segdups(de_novos):
    """ identifies de novo calls within segmental duplications.
    
    Args:
        de_novos: dataframe of candidate de novo calls
    
    Returns:
        list of booleans for whether each candidate is not in a segdup region.
    """
    
    segdups = load_segdups()
    
    # check if each candidate is not in a segdup region. This uses an interval
    # tree for efficient searching.
    statuses = []
    for i, row in de_novos.iterrows():
        chrom = row["chrom"]
        pos = row["pos"]
        
        status = len(segdups[chrom].search(pos)) == 0
        statuses.append(status)
    
    return statuses
