
#' checks for strand bias per de novo site using the ref and alt counts
#'
#' @param site dataframe of de novo variants, all at a single site, or within a
#'             single gene.
#' @export
#'
#' @return p-value for whether the forward or reverse are biased in the
#'        proportion of ref and alt alleles.
site_strand_bias <- function(site) {
    x = c(site[["REF.F"]], site[["REF.R"]], site[["ALT.F"]], site[["ALT.R"]])
    x = matrix(x, nrow=2)
    return(fisher.test(x)$p.value)
}

#' tests each site for deviation from expected behaviour
#'
#' @param de_novos dataframe of de novo variants
#'
#' @return p-value for whether the forward or reverse are biased in the
#'        proportion of ref and alt alleles.
test_sites <- function(de_novos) {
    alleles = subset(de_novos, select=c("key",
        "child.REF.F", "child.REF.R", "child.ALT.F", "child.ALT.R",
        "mother.REF.F", "mother.REF.R", "mother.ALT.F", "mother.ALT.R",
        "father.REF.F", "father.REF.R", "father.ALT.F", "father.ALT.R"))
    
    cat("splitting by site\n")
    variants = split(alleles, alleles$key)
    
    # count the ref and alt alleles for each de novo site
    cat("counting alleles\n")
    counts = lapply(variants, get_allele_counts)
    results = data.frame(matrix(unlist(counts), ncol=length(counts[[1]]), byrow=TRUE))
    names(results) = names(counts[[1]])
    results$key = names(counts)
    
    # check for overabundance of parental alt alleles using binomial test
    cat("testing parental alt overabundance\n")
    parent_counts = data.frame(results$parent.ALT, (results$parent.ALT + results$parent.REF))
    PA_pval_site = apply(parent_counts, 1, binom.test, p=ERROR_RATE, alternative="greater")
    results$PA_pval_site = as.numeric(sapply(PA_pval_site, "[", "p.value"))
    
    # check for strand bias by fishers exact test on the allele counts
    cat("testing strand bias\n")
    allele_counts = lapply(split(results, seq_along(results[, 1])), as.list)
    results$SB_pval = sapply(allele_counts, site_strand_bias)
    
    de_novos = merge(de_novos, results, by="key", all.x=TRUE)
    
    return(de_novos)
}

#' checks if the variants in a gene have more parental ALTs than expected
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return p-value for whether the forward or reverse are biased in the
#'        proportion of ref and alt alleles within each gene.
test_genes <-function(de_novos) {
    
    # exclude de novo SNVs that fail the strand bias filter, otherwise these
    # skew the parental alts within genes
    sites = de_novos[!(de_novos$SB_pval < P_CUTOFF & de_novos$var_type == "DENOVO-SNP"), ]
    
    # loop to calculate gene-specific PA values after SB filtering
    sites = subset(sites, select=c("key", "symbol",
        "mother.REF.F", "mother.REF.R", "mother.ALT.F", "mother.ALT.R",
        "father.REF.F", "father.REF.R", "father.ALT.F", "father.ALT.R"))
    
    # count the number of parental alleles within genes, and restructure data
    # for testing
    cat("counting alleles per gene\n")
    genes = split(sites, sites$symbol)
    counts = lapply(genes, get_allele_counts, gene=TRUE)
    results = data.frame(matrix(unlist(counts), ncol=length(counts[[1]]), byrow=TRUE))
    names(results) = names(counts[[1]])
    results$symbol = names(counts)
    
    # check for overabundance of parental alt alleles using binomial test
    cat("testing parental alt overabundance\n")
    parent_counts = data.frame(results$gene.ALT, (results$gene.ALT + results$gene.REF))
    PA_pval = apply(parent_counts, 1, binom.test, p=ERROR_RATE, alternative="greater")
    results$PA_pval_gene = as.numeric(sapply(PA_pval, "[", "p.value"))
    
    de_novos = merge(de_novos, results, by="symbol", all.x=TRUE)
    
    return(de_novos)
}
