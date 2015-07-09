
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
