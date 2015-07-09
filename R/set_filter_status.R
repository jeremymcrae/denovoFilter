
#' set flags for filtering, fail samples with strand bias < threshold, or any 2 of
#'  (i) both parents have ALTs
#'  (ii) site-specific parental alts < threshold,
#'  (iii) gene-specific parental alts < threshold, if > 1 sites called in gene
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return vector of true/false for whether each variant passes the filters
get_filter_status <- function(de_novos) {
    
    overall_pass = rep(TRUE, nrow(de_novos))
    
    # fail SNVs with excessive strand bias
    overall_pass[de_novos$SB_pval < P_CUTOFF & de_novos$var_type == "DENOVO-SNP"] = FALSE
    
    # find if each de novo has passed each of three different filtering strategies
    # fail sites with gene-specific parental alts, only if >1 sites called per gene
    gene_fail = de_novos$PA_pval_gene < P_CUTOFF & de_novos$count.genes > 1
    site_fail = de_novos$PA_pval_site < P_CUTOFF
    excess_alts = de_novos$min_parent_alt > 0
    
    # exclude sites that fail two of three classes
    sites = data.frame(gene_fail, site_fail, excess_alts)
    overall_pass[apply(sites, 1, sum) >= 2] = FALSE
    
    return(overall_pass)
}

#' subset down to specific columns
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return data frame with only the pertinent columns included
subset_de_novos <- function(de_novos) {
    
    de_novos = subset(de_novos, select=c("person_stable_id", "gender",
        "mother_stable_id", "father_stable_id", "chrom", "pos", "ref", "alt",
        "symbol", "var_type", "consequence", "ensg", "enst", "max_af", "pp_dnm",
        "child_alt_prp", "coding", "overall_pass",
        "child_ref_F", "child_ref_R", "child_alt_F", "child_alt_R",
        "mother_ref_F", "mother_ref_R", "mother_alt_F", "mother_alt_R",
        "father_ref_F", "father_ref_R", "father_alt_F", "father_alt_R"))
    
    return(de_novos)
}
