# script to obtain independent de novo events (ie exclude de novos that recur
# within families).

library(mupit)

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"
TRIOS_PATH = file.path(DATAFREEZE_DIR, "family_relationships.txt")

consequences = c("transcript_ablation", "splice_donor_variant",
    "splice_acceptor_variant", "stop_gained", "frameshift_variant",
    "stop_lost", "initiator_codon_variant", "transcript_amplification",
    "inframe_insertion", "inframe_deletion", "missense_variant",
    "splice_region_variant", "incomplete_terminal_codon_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_exon_variant", "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation",
    "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "feature_elongation", "feature_truncation",
    "intergenic_variant")
severity = data.frame(consequence=consequences, rank=seq(1:length(consequences)),
    stringsAsFactors=FALSE)

#' opens a file defining the DDD individuals, and their families
open_ddd_families <- function(path) {
    families = read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    families$affected = NULL
    families$path_to_vcf = NULL
    families$sex = NULL
    
    return(families)
}

#' get the most severe consequence from a list of VEP consequences
get_most_severe <- function(consequences) {
    
    best_rank = NA
    
    for (consequence in consequences) {
        
        # get consequence and severity rank in the current transcript
        temp_rank = severity$rank[severity$consequence == consequence]
        
        # check if this is the most severe consequence; prefer HGNC transcripts
        if (is.na(best_rank) | temp_rank < best_rank)  {
            best_rank = temp_rank
        }
    }
    
    return(severity$consequence[severity$rank == best_rank])
}

#' find the de novos that are recurrent within a single individual in a
#' single gene. We shall treat these as a single de novo event. Prioritise
#' including the most severe event within a gene, then take the first variant
#' left after that.
remove_within_person_recurrences <- function(de_novos) {
    
    # find the variants which are recurrent within a person in a single gene
    from_start = duplicated(de_novos[, c("person_stable_id", "symbol")])
    from_end = duplicated(de_novos[, c("person_stable_id", "symbol")], fromLast=TRUE)
    
    person_dups = from_start | from_end
    in_person_dups = de_novos[person_dups, ]
    
    # split the dataset, so we can process gene by gene
    genes = split(in_person_dups, in_person_dups$person_stable_id, in_person_dups$symbol)
    
    # construct a blank dataframe to add rows to
    dups = in_person_dups[0, ]
    dups$most_severe = logical(0)
    dups$exclude = logical(0)
    
    # pick a variant for each person, ie the first of the most severe consequence
    # This code is awkwardly written, perhaps there is a simpler way to express
    # this.
    for (gene in genes) {
        gene$most_severe = gene$consequence == get_most_severe(gene$consequence)
        gene$exclude = TRUE
        gene$exclude[gene$most_severe] = FALSE
        gene$exclude[which(gene$most_severe)[-1]] = TRUE
        dups = rbind(dups, gene)
    }
    
    # remove the selected de novos
    person_dups[person_dups] = dups$exclude
    de_novos = de_novos[!person_dups, ]
    
    return(de_novos)
}

#' remove de novos that are shared between multiple probands of a family or that
#' are recurrent in a single gene in a single person
get_independent_de_novos <- function(de_novos) {
    families = open_ddd_families(TRIOS_PATH)
    
    # merge the family IDs with the de novo events, sort the de novos first, so
    # that we can correctly match to the duplicate status
    if ("person_stable_id" %in% names(de_novos)) {
        de_novos = de_novos[order(de_novos$person_stable_id), ]
        dups = merge(de_novos, families, by.x="person_stable_id", by.y="individual_id", all.x=TRUE)
    } else {
        de_novos = de_novos[order(de_novos$person_id), ]
        dups = merge(de_novos, families, by.x="person_id", by.y="individual_id", all.x=TRUE)
    }
    
    # get the family ID, and de novo coordinates, which are sufficient to
    # identify duplicates within families
    if ("pos" %in% names(dups)) {
        dups = dups[, c("family_id", "chrom", "pos", "ref", "alt")]
    } else {
        dups = dups[, c("family_id", "chrom", "start_pos", "ref_allele", "alt_allele")]
    }
    
    # restrict ourselves to the non-duplicates (this retains the first de novo
    # for each family)
    de_novos = de_novos[!duplicated(dups), ]
    
    de_novos = remove_within_person_recurrences(de_novos)
    
    return(de_novos)
}
