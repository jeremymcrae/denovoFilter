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

def load_segdups():
    """ load all the segdup regions
    
    Returns:
        dictionary of sets of segdup regions per chromosome
    """
    
    segdup_path = resource_filename(__name__, "data/segdup_regions.gz")
    
    segdups = {}
    with gzip.open(segdup_path, "rt") as handle:
        
        header = handle.readline()
        
        for line in handle:
            line = line.strip().split("\t")
            chrom = line[0].strip("chr")
            start = int(line[1])
            end = int(line[2])
            
            if chrom not in segdups:
                segdups[chrom] = set([])
            
            segdups[chrom].add((start, end))
    
    return segdups

def exclude_segdups(de_novos):
    """ removes de novo calls within segmental duplications.
    
    Args:
        de_novos: dataframe of candidate de novo calls
    
    Returns:
        dataframe of candidate de novo calls with sites in segmental
        duplications removed.
    """
    
    segdups = load_segdups()
    
    # for each de novo, check if it is in a segdup region. This is definitely
    #  not the most efficient way to check this.
    de_novos["in_segdup"] = True
    for row_index, de_novo in de_novos.iterrows():
        chrom = de_novo["chrom"]
        pos = de_novo["pos"]
        
        in_segdup = any([ pos >= x[0] and pos <= x[1] for x in segdups[chrom] ])
        de_novos.loc[row_index, "in_segdup"] = in_segdup
    
    return de_novos[-de_novos["in_segdup"]]
