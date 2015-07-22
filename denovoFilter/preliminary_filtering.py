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

def fix_maf(de_novos):
    """ cleans up the max AF entries in the de novo dataframe
    
    Args:
        de_novos: DataFrame of candidate de novos, include a column of "max_af"
    
    Returns:
        column of max_af values
    """
    
    max_af = de_novos["max_af"].copy()
    
    # one de novo has a comma separated list of MAF values (both above 0.1)
    max_af[max_af.str.contains(",")] = 1
    
    # fix the blank and null max AF values
    max_af[(max_af == "") | (max_af == ".")] = 0
    max_af = max_af.convert_objects(convert_numeric=True)
    
    return max_af

def preliminary_filtering(de_novos):
    """run some preliminary filtering of de novos
    
    We want to filter out de novos with high MAF, where they are not present in
    the child VCF, or are present in the parental VCFs
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        data frame of de novos, but where we have excluded sites with high
        allele frequency and where the variant appears in one or more parents.
    """
    max_af = fix_maf(de_novos)
    de_novos = de_novos[(max_af < 0.01) & max_af.notnull()]
    
    # keep sites in child vcf, and not in parental vcfs
    de_novos = de_novos[(de_novos["in_child_vcf"] == 1) & (de_novos["in_father_vcf"] == 0) & (de_novos["in_mother_vcf"] == 0)]
    
    # remove de_novos in samples with >> too many DNMs, focus on too many DNMs at
    # high quality
    # NEED TO UPDATE WITH CAROLINE'S NEW SAMPLE FILE LIST, OR USE PRE-FILTERED SET OF DNMS
    sample_fails = [276227, 258876, 273778, 258921, 272110, 260337, 264083,
        264084, 269602, 265624]
    de_novos = de_novos[-de_novos["decipher_id"].isin(sample_fails)]
    
    # Annotate de_novos hitting coding exons or splice sites
    coding_splicing = ["coding_sequence_variant", "frameshift_variant",
        "inframe_deletion", "inframe_insertion", "initiator_codon_variant",
        "missense_variant", "protein_altering_variant", "splice_acceptor_variant",
        "splice_donor_variant", "splice_region_variant", "start_lost",
        "stop_gained", "stop_lost", "synonymous_variant"]
    
    de_novos["coding"] = de_novos["consequence"].isin(coding_splicing)
    
    return de_novos
