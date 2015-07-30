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

from denovoFilter.constants import P_CUTOFF

def get_filter_status(de_novos):
    """ set flags for filtering, fail samples with strand bias < threshold, or any 2 of
     (i) both parents have ALTs
     (ii) site-specific parental alts < threshold,
     (iii) gene-specific parental alts < threshold, if > 1 sites called in gene
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        vector of true/false for whether each variant passes the filters
    """
    
    overall_pass = pandas.Series([True] * len(de_novos))
    
    # fail SNVs with excessive strand bias
    overall_pass[(de_novos["SB_pval"] < P_CUTOFF) & (de_novos["var_type"] == "DENOVO-SNP")] = False
    
    # find if each de novo has passed each of three different filtering strategies
    # fail sites with gene-specific parental alts, only if >1 sites called per gene
    gene_fail = (de_novos["PA_pval_gene"] < P_CUTOFF) & (de_novos["count.genes"] > 1)
    site_fail = (de_novos["PA_pval_site"] < P_CUTOFF)
    excess_alts = de_novos["min_parent_alt"] > 0
    
    # exclude sites that fail two of three classes
    sites = pandas.DataFrame({"gene": gene_fail, "site": site_fail, "alts": excess_alts})
    overall_pass[sites.sum(axis=1) >= 2] = False
    
    return overall_pass

def subset_de_novos(de_novos):
    """ subset down to specific columns
    
    Args:
        de_novos dataframe of de novo variants
    
     Returns:
        data frame with only the pertinent columns included
    """
    
    de_novos = de_novos[["person_stable_id", "sex", "chrom", "pos", "ref", "alt",
        "symbol", "var_type", "consequence", "max_af", "pp_dnm",
        "child_ref_F", "child_ref_R", "child_alt_F", "child_alt_R",
        "mother_ref_F", "mother_ref_R", "mother_alt_F", "mother_alt_R",
        "father_ref_F", "father_ref_R", "father_alt_F", "father_alt_R"]]
    
    return de_novos
