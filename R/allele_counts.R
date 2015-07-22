
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
        de_novos[[paste(member, "_ref_F", sep="")]] = as.numeric(sapply(dp4, "[", 1))
        de_novos[[paste(member, "_ref_R", sep="")]] = as.numeric(sapply(dp4, "[", 2))
        de_novos[[paste(member, "_alt_F", sep="")]] = as.numeric(sapply(dp4, "[", 3))
        de_novos[[paste(member, "_alt_R", sep="")]] = as.numeric(sapply(dp4, "[", 4))
    }
    
    # get the minimum alternate allele count from the parents
    alts = data.frame(de_novos$mother_alt_F + de_novos$mother_alt_R,
        de_novos$father_alt_F + de_novos$father_alt_R)
    de_novos$min_parent_alt = apply(alts, 1, min)
    
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
    
    #  sum the counts from different de novo calls
    counts = colSums(vars[, 2:ncol(vars)])
    
    # count parental alt and ref
    parent_alt = sum(counts[c("mother_alt_F", "father_alt_F", "mother_alt_R", "father_alt_R")])
    parent_ref = sum(counts[c("mother_ref_F", "father_ref_F", "mother_ref_R", "father_ref_R")])
    values = list()
    
    if (!gene) {
        values$ref_F = sum(counts[c("child_ref_F", "mother_ref_F", "father_ref_F")])
        values$ref_R = sum(counts[c("child_ref_R", "mother_ref_R", "father_ref_R")])
        values$alt_F = sum(counts[c("child_alt_F", "mother_alt_F", "father_alt_F")])
        values$alt_R = sum(counts[c("child_alt_R", "mother_alt_R", "father_alt_R")])
        values$parent_alt = parent_alt
        values$parent_ref = parent_ref
    } else {
        values$gene_alt = parent_alt
        values$gene_ref = parent_ref
    }
    
    return(values)
}
