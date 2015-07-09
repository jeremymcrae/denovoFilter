
#' run some preliminary filtering of de novos
#'
#' We want to filter out de novos with high MAF, where they are not present in
#' the child VCF, or are present in the parental VCFs
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return data frame of de novos, but where we have excluded sites with high
#'         allele frequency and where the variant appears in one or more parents.
preliminary_filtering <- function(de_novos) {
    de_novos$max_af = fix_maf(de_novos)
    
    # keep sites in child vcf, and not in parental vcfs
    de_novos = de_novos[de_novos$in_child_vcf == 1 & de_novos$in_father_vcf == 0 & de_novos$in_mother_vcf == 0, ]
    
    # remove sites with high population frequencies
    de_novos = de_novos[de_novos$max_af < 0.01 & !is.na(de_novos$max_af), ]
    
    # remove de_novos in samples with >> too many DNMs, focus on too many DNMs at
    # high quality
    # NEED TO UPDATE WITH CAROLINE'S NEW SAMPLE FILE LIST, OR USE PRE-FILTERED SET OF DNMS
    sample.fails = c("276227", "258876", "273778", "258921", "272110", "260337",
        "264083", "264084", "269602", "265624")
    de_novos = de_novos[!de_novos$decipher_id %in% sample.fails, ]
    
    # Annotate de_novos hitting coding exons or splice sites
    coding_splicing = c("coding_sequence_variant", "frameshift_variant",
        "inframe_deletion", "inframe_insertion", "initiator_codon_variant",
        "missense_variant", "protein_altering_variant", "splice_acceptor_variant",
        "splice_donor_variant", "splice_region_variant", "start_lost",
        "stop_gained", "stop_lost", "synonymous_variant")
        
    de_novos$coding = de_novos$consequence %in% coding_splicing
    
    return(de_novos)
}

#' fix a few MAF values that are blank, or have lists of MAF values
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return vector of max allele frequency values for each site
fix_maf <- function(de_novos) {
    
    # annotate with numeric max allele freq
    max_af = de_novos$max_af
    
    # fix the blank and null max AF values
    max_af[max_af == "" | max_af == "."] = 0
    
    # one de novo has a comma separated list of MAF values (both above 0.1)
    max_af[grepl(",", max_af)] = NA
    max_af = as.numeric(max_af)
    
    return(max_af)
}
