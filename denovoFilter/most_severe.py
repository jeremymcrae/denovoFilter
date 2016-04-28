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

consequences = ["transcript_ablation", "splice_donor_variant",
    "splice_acceptor_variant", "stop_gained", "frameshift_variant",
    "initiator_codon_variant", "stop_lost", "start_lost", "transcript_amplification",
    'conserved_exon_terminus_variant',
    "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
    "splice_region_variant", "incomplete_terminal_codon_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_exon_variant", "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation",
    "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "feature_elongation", "feature_truncation",
    "intergenic_variant"]
severity = pandas.DataFrame({"consequence": consequences, "rank": range(len(consequences))})

def get_most_severe(consequences):
    """ get the most severe consequence from a list of VEP consequences
    
    Args:
        consequences: vector of VEP consequence strings
    
    Returns:
        single string for the most severe consequence
    """
    
    best_rank = None
    
    for consequence in consequences:
        # get consequence and severity rank in the current transcript
        temp_rank = int(severity["rank"][severity["consequence"] == consequence])
        
        # check if this is the most severe consequence; prefer HGNC transcripts
        if best_rank is None or temp_rank < best_rank:
            best_rank = temp_rank
    
    return list(severity["consequence"][severity["rank"] == best_rank])[0]
