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

from __future__ import division

import pandas

def extract_alt_and_ref_counts(de_novos):
    """extract ALT and REF counts of forward and reverse reads for trio members
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        data frame of de novos, but with an extra columns for the read depths
        of forward and reverse reads for each member of the trio.
    """
    
    members = ["child", "father", "mother"]
    
    counts = de_novos[['person_stable_id', 'chrom', 'pos', 'ref', 'alt', 'symbol']].copy()
    
    for member in members:
        dp4 = de_novos["dp4_{0}".format(member)].str.split(",")
        
        # split the forward and reverse counts by ref and alt for each member
        counts["{0}_ref_F".format(member)] = [ int(x[0]) for x in dp4 ]
        counts["{0}_ref_R".format(member)] = [ int(x[1]) for x in dp4 ]
        counts["{0}_alt_F".format(member)] = [ int(x[2]) for x in dp4 ]
        counts["{0}_alt_R".format(member)] = [ int(x[3]) for x in dp4 ]
    
    # get the minimum alternate allele count from the parents
    alts = pandas.DataFrame({"a": counts[["mother_alt_F", "mother_alt_R"]].sum(axis=1),
        "b": counts[["father_alt_F", "father_alt_R"]].sum(axis=1)})
    counts["min_parent_alt"] = alts.min(axis=1)
    
    return counts

def get_recurrent_genes(de_novos):
    """ count candidate de novo mutations per gene
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        dictionary of counts indexed by HGNC symbol
    """
    
    # count of number of candidates per gene
    counts = de_novos["symbol"].value_counts()
    
    return counts.index[counts > 1]

def get_allele_counts(variants, gene=False):
    """ counts REF and ALT alleles in reads for candidate de novo events
    
    Args:
        variants: dataframe of de novo variants, all at a single site, or within
            a single gene.
        gene: whether the set of variants of all within a single gene (rather
            than at a site)
    
    Returns:
        dictionary of forward reference, forward alternate.
    """
    
    counts = pandas.to_numeric(variants.sum(), errors='coerce')
    
    # count parental alt and ref
    parent_alt = sum(counts[["mother_alt_F", "father_alt_F", "mother_alt_R", "father_alt_R"]])
    parent_ref = sum(counts[["mother_ref_F", "father_ref_F", "mother_ref_R", "father_ref_R"]])
    values = {}
    
    if not gene:
        values["ref_F"] = sum(counts[["child_ref_F", "mother_ref_F", "father_ref_F"]])
        values["ref_R"] = sum(counts[["child_ref_R", "mother_ref_R", "father_ref_R"]])
        values["alt_F"] = sum(counts[["child_alt_F", "mother_alt_F", "father_alt_F"]])
        values["alt_R"] = sum(counts[["child_alt_R", "mother_alt_R", "father_alt_R"]])
        values["parent_alt"] = parent_alt
        values["parent_ref"] = parent_ref
    else:
        values["gene_alt"] = parent_alt
        values["gene_ref"] = parent_ref
    
    return values

def get_depths_and_proportions(counts):
    """ get the read depths, alt depths and alt proportions for the trio members
    
    Args:
        counts: pandas dataframe containing colums of child_ref_F, child alt_R
            etc
    
    Returns:
        dataframe including additional alt depth, total depth, and alt
        proportion columns for each trio member.
    """
    
    # get the trio depths
    child_depth = counts[["child_ref_F", "child_ref_R", "child_alt_F", "child_alt_R"]].sum(axis=1)
    dad_depth = counts[["father_ref_F", "father_ref_R", "father_alt_F", "father_alt_R"]].sum(axis=1)
    mom_depth = counts[["mother_ref_F", "mother_ref_R", "mother_alt_F", "mother_alt_R"]].sum(axis=1)
    
    # get the depth of the alternate alleles
    child_alts = counts[["child_alt_F", "child_alt_R"]].sum(axis=1)
    dad_alts = counts[["father_alt_F", "father_alt_R"]].sum(axis=1)
    mom_alts = counts[["mother_alt_F", "mother_alt_R"]].sum(axis=1)
    
    values = pandas.DataFrame({
        'child_depth': child_depth, 'dad_depth': dad_depth, 'mom_depth': mom_depth,
        'child_alts': child_alts, 'dad_alts': dad_alts, 'mom_alts': mom_alts,
        'child_prp': child_alts/child_depth, 'dad_prp': dad_alts/dad_depth,
        'mom_prp': mom_alts/mom_depth}, index=counts.index)
    
    return values
