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

from denovoFilter.preliminary_filtering import preliminary_filtering
from denovoFilter.exclude_segdups import check_segdups
from denovoFilter.allele_counts import extract_alt_and_ref_counts, \
    get_depths_and_proportions

def load_missing_indels(candidates_path):
    """ load the missed candidates dataset
    
    Args:
        candidates_path:
        families_path:
    
    Returns:
        pandas dataframe of candidate de novo sites
    """
    
    missed = pandas.read_table(candidates_path, na_filter=False)
    
    # rename columns so that later processing works smoothly
    missed["in_child_vcf"] = 1
    missed["in_mother_vcf"] = 0
    missed["in_father_vcf"] = 0
    missed["pp_dnm"] = None
    
    return missed

def filter_missing_indels(candidates):
    """ filter the candidate missing indels.
    
    We have a set of sites that have been called in the child, but not in the
    parents. These are the candidates for de novo mutations. Many of these sites
    have been checked by denovogear, for which we have a different filtering
    process. This function filters candidate sites which have not been examined
    by denovogear. Denovogear has a reduced sensitivity for indels, so this
    filtering only examines candidate indel sites not examined by denovogear.
    
    Args:
        candidates: pandas dataframe of de novo indel sites
    
    Returns:
        dataframe of candidate sites that pass the required criteria.
    """
    
    counts = extract_alt_and_ref_counts(candidates)
    depths = get_depths_and_proportions(counts)
    
    depths["min_parent_depth"] = depths[["mom_depth", "dad_depth"]].min(axis=1)
    depths["max_parent_proportion"] = depths[["mom_prp", "dad_prp"]].max(axis=1)
    
    # apply the filtering criteria for the missing indels
    good_depth = depths["child_alts"] > 2
    low_parental_alt = counts["min_parent_alt"] < 2
    good_parental_depth = depths["min_parent_depth"] > 7
    good_parental_proportion = depths["max_parent_proportion"] < 0.1
    good_child_proportion = depths["child_prp"] > 0.2
    
    return good_depth & low_parental_alt & good_parental_depth & \
        good_parental_proportion & good_child_proportion
