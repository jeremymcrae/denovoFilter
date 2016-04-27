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

from denovoFilter.allele_counts import extract_alt_and_ref_counts, \
    get_recurrent_genes
from denovoFilter.site_deviations import test_sites, test_genes
from denovoFilter.constants import P_CUTOFF

def filter_denovogear_sites(de_novos, pass_status, include_metrics=False):
    """ set flags for filtering, fail samples with strand bias < threshold, or any 2 of
     (i) both parents have ALTs
     (ii) site-specific parental alts < threshold,
     (iii) gene-specific parental alts < threshold, if > 1 sites called in gene
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        vector of true/false for whether each variant passes the filters
    """
    
    counts = extract_alt_and_ref_counts(de_novos)
    recurrent = get_recurrent_genes(de_novos)
    
    # check if sites deviate from expected strand bias and parental alt depths
    strand_bias, parental_site_bias = test_sites(counts, pass_status)
    parental_gene_bias = test_genes(counts, strand_bias, pass_status)
    
    # fail SNVs with excessive strand bias. Don't check strand bias in indels.
    overall_pass = (strand_bias >= P_CUTOFF) | (de_novos["ref"].str.len() != 1) | \
        (de_novos["alt"].str.len() != 1)
    
    # find if each de novo has passed each of three different filtering strategies
    # fail sites with gene-specific parental alts, only if >1 sites called per gene
    gene_fail = (parental_gene_bias < P_CUTOFF) & de_novos["symbol"].isin(recurrent)
    site_fail = (parental_site_bias < P_CUTOFF)
    excess_alts = counts["min_parent_alt"] > 0
    
    # exclude sites that fail two of three classes
    sites = pandas.DataFrame({"gene": gene_fail, "site": site_fail, "alts": excess_alts})
    overall_pass[sites.sum(axis=1) >= 2] = False
    
    if not include_metrics:
        return overall_pass
    else:
        return overall_pass, strand_bias, parental_site_bias, parental_gene_bias, \
            counts["min_parent_alt"]
