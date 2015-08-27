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

import json

def change_conserved_last_base_consequence(de_novos, last_base_path):
    """ reannotate the consequence of conserved sites at the end of exons
    
    Args:
        de_novos: pandas dataframe of candidate de novos
        last_base_path: path to file listing conserved sites at the end of exons
    
    Returns:
        dataframe, but with the conserved sites consequence reannotated as
        "conserved_exon_terminus_variant"
    """
    
    with open(last_base_path, "r") as handle:
        sites_json = json.load(handle)
    
    # get a dictionary of sites per chromosome
    sites = {}
    for (chrom, pos) in sites_json:
        if chrom not in sites:
            sites[chrom] = set([])
        
        sites[chrom].add(int(pos))
    
    # find the sites for each chromosome at conserved last base of exons, and
    # give them a new consequence value.
    for chrom in sites:
        positions = sites[chrom]
        
        # find which sites are for the correct chromosome, at one of the
        # identified sites, and have a single base alleles (i.e "C" or "G")
        in_chrom = de_novos["chrom"] == chrom
        in_pos = de_novos["pos"].isin(positions)
        single_base = (de_novos["ref"].str.len() == 1) & (de_novos["alt"].str.len() == 1)
        
        cq = de_novos["consequence"].copy()
        cq[in_chrom & in_pos & single_base] = "conserved_exon_terminus_variant"
        de_novos["consequence"] = cq
    
    return de_novos
