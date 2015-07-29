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

import pandas
import numpy

def extract_alt_and_ref_counts(de_novos):
    """extract ALT and REF counts of forward and reverse reads for the trio members
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        data frame of de novos, but with an extra columns for the read depths
        of forward and reverse reads for each member of the trio.
    """
    
    members = ["child", "father", "mother"]
    
    for member in members:
        dp4 = de_novos["dp4_{0}".format(member)].str.split(",")
        
        # split the forward and reverse counts by ref and alt for each member
        de_novos["{0}_ref_F".format(member)] = [ int(x[0]) for x in dp4 ]
        de_novos["{0}_ref_R".format(member)] = [ int(x[1]) for x in dp4 ]
        de_novos["{0}_alt_F".format(member)] = [ int(x[2]) for x in dp4 ]
        de_novos["{0}_alt_R".format(member)] = [ int(x[3]) for x in dp4 ]
    
    # get the minimum alternate allele count from the parents
    alts = pandas.DataFrame({"a": de_novos["mother_alt_F"] + de_novos["mother_alt_R"],
        "b": de_novos["father_alt_F"] + de_novos["father_alt_R"]})
    de_novos["min_parent_alt"] = alts.min(axis=1)
    
    return de_novos

def count_gene_recurrence(de_novos):
    """ count of number of sites called per gene
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        data frame of de novos, but with an extra column indicating the number
        of candidate sites in the gene.
    """
    
    de_novos["key"] = de_novos["chrom"].map(str) + de_novos["pos"].map(str) + de_novos["alt"].map(str)
    
    # count of number times that gene is called
    counts = de_novos["symbol"].value_counts()
    counts = pandas.DataFrame({"symbol": counts.index, "count.genes": counts})
    de_novos = de_novos.merge(counts)
    
    return de_novos

def get_allele_counts(vars, gene=False):
    """ counts REF and ALT alleles in reads for candidate de novo events
    
    Args:
        vars: dataframe of de novo variants, all at a single site, or within a
             single gene.
        gene: whether the set of variants of all within a single gene (rather
            than at a site)
    
    Returns:
        list of forward reference, forward alternate.
    """
    
    vars = vars.drop([vars.columns[0]], axis=1)
    counts = vars.sum()
    
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

def get_depths_and_proportions(de_novos):
    """ get the read depths, alt depths and alt proportions for the trio members
    
    Args:
        de_novos: pandas dataframe containing colums of child_ref_F, child alt_R
            etc
    
    Returns:
        dataframe including additional alt depth, total depth, and alt
        proportion columns for each trio member.
    """
    
    # get the trio depths
    de_novos["child_depth"] = de_novos[["child_ref_F", "child_ref_R", "child_alt_F", "child_alt_R"]].sum(axis=1)
    de_novos["father_depth"] = de_novos[["father_ref_F", "father_ref_R", "father_alt_F", "father_alt_R"]].sum(axis=1)
    de_novos["mother_depth"] = de_novos[["mother_ref_F", "mother_ref_R", "mother_alt_F", "mother_alt_R"]].sum(axis=1)
    
    # get the depth of the alternate alleles
    de_novos["child_alts"] = de_novos[["child_alt_F", "child_alt_R"]].sum(axis=1)
    de_novos["father_alts"] = de_novos[["father_alt_F", "father_alt_R"]].sum(axis=1)
    de_novos["mother_alts"] = de_novos[["mother_alt_F", "mother_alt_R"]].sum(axis=1)
    
    # calculate the alt proportions
    de_novos["child_prp"] = de_novos["child_alts"]/de_novos["child_depth"]
    de_novos["mother_prp"] = de_novos["mother_alts"]/de_novos["mother_depth"]
    de_novos["father_prp"] = de_novos["father_alts"]/de_novos["father_depth"]
    
    return de_novos
