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

from pandas import Series

def fix_maf(max_af):
    """ cleans up the max AF entries in the de novo dataframe
    
    Args:
        max_af: pandas Series of maximum allele frequencies from the table of
            candidate de novos.
    
    Returns:
        modified pandas Series of max_af values
    """
    
    max_af = max_af.copy()
    
    # some de novos have comma-separated lists of maf values. We select the
    # first value (which correspononds to the alternate allele, the additional
    # alleles are from denovogear selecting all possibly alternates at a
    # candidate de novo site)
    max_af = [ x.split(',')[0] if type(x) == str else None for x in max_af ]
    
    # fix the null max AF values
    missing = set(['', '.', 'missing'])
    max_af = [ x if x not in missing else 0 for x in max_af ]
    
    # replace missing maf values, and convert everything to floats
    max_af = [ float(x) if x is not None else 0.0 for x in max_af ]
    
    return Series(max_af)

def preliminary_filtering(de_novos, sample_fails=None, maf_cutoff=0.01):
    """run some preliminary filtering of de novos.
    
    We want to filter out de novos with high MAF, where they are not present in
    the child VCF, or are present in the parental VCFs
    
    Args:
        de_novos: dataframe of de novo variants
        sample_fails: list of IDs for samples who have failed QC checks (due to
            having too many de novo calls).
        maf_cutoff: sites have to have a minor allele frequency below this.
    
    Returns:
        data frame of de novos, but where we have excluded sites with high
        allele frequency and where the variant appears in one or more parents.
    """
    
    max_af = fix_maf(de_novos['max_af'])
    passes_maf = (max_af <= maf_cutoff) & max_af.notnull()
    
    # keep sites in child vcf, and not in parental vcfs
    appears_de_novo = (de_novos["in_child_vcf"] == 1) & (de_novos["in_father_vcf"] == 0) & (de_novos["in_mother_vcf"] == 0)
    
    good_samples = [True] * len(de_novos)
    if sample_fails is not None:
        # remove samples with too many candidates
        good_samples = ~de_novos["person_stable_id"].isin(sample_fails)
    
    return passes_maf & appears_de_novo & good_samples

def check_coding(de_novos, cq_name="consequence"):
    """check whether the consequence is for a coding consequence
    
    Args:
        de_novos: dataframe of de novo variants
        cq_name: name of the column for the consequence
    
    Returns:
        data frame with an additional column to show if the site is coding
    """
    
    # annotate candidates within coding exons or splice sites
    coding_splicing = ["coding_sequence_variant", "frameshift_variant",
        "inframe_deletion", "inframe_insertion", "initiator_codon_variant",
        "missense_variant", "protein_altering_variant", "splice_acceptor_variant",
        "splice_donor_variant", "splice_region_variant", "start_lost",
        "stop_gained", "stop_lost", "synonymous_variant",
        'conserved_exon_terminus_variant']
    
    return de_novos[cq_name].isin(coding_splicing)
