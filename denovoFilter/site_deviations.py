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

import scipy.stats
import pandas

from denovoFilter.allele_counts import get_allele_counts
from denovoFilter.constants import P_CUTOFF, ERROR_RATE

def site_strand_bias(site):
    """ checks for strand bias per de novo site using the ref and alt counts
    
    Args:
        site: dataframe of de novo variants, all at a single site, or within a
            single gene.
    
    Returns:
        p-value for whether the forward or reverse are biased in the
        proportion of ref and alt alleles.
    """
    
    x = [[site["ref_F"], site["ref_R"]], [site["alt_F"], site["alt_R"]]]
    
    try:
        return scipy.stats.fisher_exact(x)[1]
    except ValueError:
        return float(1)

def test_sites(de_novos):
    """ tests each site for deviation from expected behaviour
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        p-value for whether the forward or reverse are biased in the proportion
        of ref and alt alleles.
    """
    
    alleles = de_novos[["key",
        "child_ref_F", "child_ref_R", "child_alt_F", "child_alt_R",
        "mother_ref_F", "mother_ref_R", "mother_alt_F", "mother_alt_R",
        "father_ref_F", "father_ref_R", "father_alt_F", "father_alt_R"]]
    
    alleles = alleles.convert_objects(convert_numeric=True)
    variants = alleles.groupby("key")
    
    # count the ref and alt alleles for each de novo site
    counts = [ get_allele_counts(site) for name, site in variants ]
    results = pandas.DataFrame(counts)
    results["key"] = [ name for name, site in variants ]
    
    # check for overabundance of parental alt alleles using binomial test
    parent_counts = pandas.DataFrame({"alt": results["parent_alt"], "ref": results["parent_ref"]})
    results["PA_pval_site"] = parent_counts.apply(scipy.stats.binom_test, axis=1, p=ERROR_RATE)
    
    # check for strand bias by fishers exact test on the allele counts
    results["SB_pval"] = results.apply(site_strand_bias, axis=1)
    
    de_novos = de_novos.merge(results, how="left")
    
    return de_novos

def test_genes(de_novos):
    """ checks if the variants in a gene have more parental ALTs than expected
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        p-value for whether the forward or reverse are biased in the proportion
        of ref and alt alleles within each gene.
    """
    
    # exclude de novo SNVs that fail the strand bias filter, otherwise these
    # skew the parental alts within genes
    sites = de_novos[~(de_novos["SB_pval"] < P_CUTOFF) & (de_novos["var_type"] == "DENOVO-SNP")]
    
    # loop to calculate gene-specific PA values after SB filtering
    sites = sites[["symbol",
        "mother_ref_F", "mother_ref_R", "mother_alt_F", "mother_alt_R",
        "father_ref_F", "father_ref_R", "father_alt_F", "father_alt_R"]]
    
    sites = sites.convert_objects(convert_numeric=True)
    
    # count the number of parental alleles within genes, and restructure data
    # for testing
    genes = sites.groupby("symbol")
    counts = [ get_allele_counts(site, gene=True) for name, site in genes ]
    results = pandas.DataFrame(counts)
    results["symbol"] = [ name for name, site in genes ]
    
    # check for overabundance of parental alt alleles using binomial test
    parent_counts = pandas.DataFrame({"alt": results["gene_alt"], "ref": results["gene_ref"]})
    results["PA_pval_gene"] = parent_counts.apply(scipy.stats.binom_test, axis=1, p=ERROR_RATE)
    
    de_novos = de_novos.merge(results, how="left")
    
    return de_novos
