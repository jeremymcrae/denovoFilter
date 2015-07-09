

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

#' adds gene symbols to variants lacking them, using a mupit function
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return data frame of de novos, but with additional annotations for many
#'         variants previously lacking a HGNC symbol.
fix_missing_gene_symbols <- function(de_novos) {
    
    # get the variants with no gene annotation, ensure chrom, start and stop
    # positions columns exist
    missing_genes = de_novos[de_novos$symbol == "", ]
    missing_genes$start_pos = as.character(missing_genes$pos)
    missing_genes$end_pos = as.character(as.numeric(missing_genes$start_pos) +
        (nchar(missing_genes$ref) - 1))
    
    # find the HGNC symbols (if any) for the variants
    hgnc_symbols = apply(missing_genes, 1, mupit::get_gene_id_for_variant, verbose=TRUE)
    
    # add the HGNC symbols to the rows that need it
    de_novos$symbol[de_novos$symbol == ""] = hgnc_symbols
    
    # 360 out of 17000 de novos still lack HGNC symbols. Their consequences are:
    #
    #   consequence                count
    #   ========================   =====
    #   downstream_gene_variant      17
    #   intergenic_variant          259
    #   regulatory_region_variant    63
    #   upstream_gene_variant        27
    #
    # In spot checks, these are sufficiently distant from genes that we can't
    # add them to the analysis of their nearest gene. We shall analyse these
    # per site by giving them mock gene symbols.
    missing_genes = de_novos[de_novos$symbol == "", ]
    de_novos$symbol[de_novos$symbol == ""] = paste("fake_symbol.",
        missing_genes$chrom, "_", missing_genes$pos, sep="")
    
    return(de_novos)
}

#' extract ALT and REF counts of forward and reverse reads for the members of
#' each trio
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return data frame of de novos, but with an extra columns for the read
#'         depths of forward and reverse reads for each member of the trio.
extract_alt_and_ref_counts <- function(de_novos) {
    
    members = c("child", "father", "mother")
    
    for (member in members) {
        dp4 = strsplit(de_novos[[paste("dp4_", member, sep="")]], ",")
        
        # split the forward and reverse counts by ref and alt for each member
        de_novos[[paste(member, ".REF.F", sep="")]] = as.numeric(sapply(dp4, "[", 1))
        de_novos[[paste(member, ".REF.R", sep="")]] = as.numeric(sapply(dp4, "[", 2))
        de_novos[[paste(member, ".ALT.F", sep="")]] = as.numeric(sapply(dp4, "[", 3))
        de_novos[[paste(member, ".ALT.R", sep="")]] = as.numeric(sapply(dp4, "[", 4))
    }
    
    de_novos$count.child.alt = de_novos$child.ALT.F + de_novos$child.ALT.R
    
    # get the minimum alternate allele count from the parents
    alts = data.frame(de_novos$mother.ALT.F + de_novos$mother.ALT.R,
        de_novos$father.ALT.F + de_novos$father.ALT.R)
    de_novos$min.parent.ALT = apply(alts, 1, min)
    
    return(de_novos)
}

#' count of number of sites called per gene
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return data frame of de novos, but with an extra column indicating the
#        number of candidate sites in the gene.
count_gene_recurrence <- function(de_novos) {
    
    de_novos$key = paste(de_novos$chrom, de_novos$pos, de_novos$alt, sep="_")
    
    # count of number times that gene is called
    table.genes = table(de_novos$symbol)
    de_novos$count.genes = table.genes[match(de_novos$symbol, row.names(table.genes))]
    
    return(de_novos)
}

#' counts REF and ALT alleles in reads for candidate de novo events
#'
#' @param vars dataframe of de novo variants, all at a single site, or within a
#'             single gene.
#' @param gene whether the set of variants of all within a single gene (rather
#'        than at a site)
#' @export
#'
#' @return list of forward reference, forward alternate.
get_allele_counts <- function(vars, gene=FALSE) {
    
    # count parental alt and ref
    parent.ALT = sum(vars[, c("mother.ALT.F", "father.ALT.F", "mother.ALT.R", "father.ALT.R")])
    parent.REF = sum(vars[, c("mother.REF.F", "father.REF.F", "mother.REF.R", "father.REF.R")])
    values = list()
    
    if (!gene) {
        values$REF.F = sum(vars[, c("child.REF.F", "mother.REF.F", "father.REF.F")])
        values$REF.R = sum(vars[, c("child.REF.R", "mother.REF.R", "father.REF.R")])
        values$ALT.F = sum(vars[, c("child.ALT.F", "mother.ALT.F", "father.ALT.F")])
        values$ALT.R = sum(vars[, c("child.ALT.R", "mother.ALT.R", "father.ALT.R")])
        values$parent.ALT = parent.ALT
        values$parent.REF = parent.REF
    } else {
        values$gene.ALT = parent.ALT
        values$gene.REF = parent.REF
    }
    
    return(values)
}

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
    excess_alts = de_novos$min.parent.ALT > 0
    
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
        "child.REF.F", "child.REF.R", "child.ALT.F", "child.ALT.R",
        "mother.REF.F", "mother.REF.R", "mother.ALT.F", "mother.ALT.R",
        "father.REF.F", "father.REF.R", "father.ALT.F", "father.ALT.R"))
    
    return(de_novos)
}
