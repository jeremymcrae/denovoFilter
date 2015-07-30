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

def remove_within_person_recurrences(de_novos):
    """ remove de novos recurrent in a gene within individuals.
    
    Find the de novos that are recurrent within a single individual in a
    single gene. We shall treat these as a single de novo event. Prioritise
    including the most severe event within a gene, then take the first variant
    left after that.
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        dataframe with duplicated sites per gene removed for each individual
    """
    
    # find the variants which are recurrent within a person in a single gene
    from_start = de_novos[["person_stable_id", "symbol"]].duplicated()
    from_end = de_novos[["person_stable_id", "symbol"]].duplicated(take_last=True)
    
    person_dups = from_start | from_end
    in_person_dups = de_novos[person_dups]
    
    # split the dataset, so we can process gene by gene
    genes = in_person_dups.groupby(["person_stable_id", "symbol"])
    
    # construct a blank dataframe to add rows to
    dups = pandas.DataFrame(columns=in_person_dups.columns)
    dups["most_severe"] = False
    dups["exclude"] = False
    
    # pick a variant for each person, ie the first of the most severe consequence
    # This code is awkwardly written, perhaps there is a simpler way to express
    # this.
    for key, gene in genes:
        gene["most_severe"] = gene["consequence"] == get_most_severe(gene["consequence"])
        gene["exclude"] = True
        gene["exclude"][gene["most_severe"]] = False
        gene["exclude"][gene["most_severe"][1:].index] = True
        dups = dups.append(gene)
    
    # remove the selected de novos
    person_dups[person_dups] = dups["exclude"]
    de_novos = de_novos[-person_dups]
    
    return de_novos

def get_independent_de_novos(de_novos, trios_path):
    """ exclude de novos that would otherwise lead to double counting.
    
    Remove de novos that are shared between multiple probands of a family or
    that are recurrent in a single gene in a single person.
    
    Args:
        de_novos: dataframe of de novo variants
        trios_path
    
    Returns:
        dataframe with duplicated sites removed
    """
    
    families = pandas.read_table(trios_path, na_filter=False)
    families = de_novos.merge(families, how="left", left_on="person_stable_id", right_on="individual_id")
    
    # get the family ID, and de novo coordinates, which are sufficient to
    # identify duplicates within families
    families = families[["family_id", "chrom", "pos", "ref", "alt"]]
    families.index = de_novos.index.copy()
    
    # restrict ourselves to the non-duplicates (this retains the first de novo
    # for each family)
    without_dups = de_novos[-families.duplicated()]
    without_recurrences = remove_within_person_recurrences(without_dups)
    
    return without_recurrences
