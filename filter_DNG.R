# filter DNG callset initially on basis of

# - in_child_vcf, not in parent vcf
# - max_af<0.01
# - remove dodgy samples with too many DNM calls

# annotate DNG callset

# - with coding/splicing
# - with presence in DDG2P
# - with strand-specitic counts of REF and ALT reads
# - with site-specific strand bias (SB) and parental alt (PA) frequency p values against null
# - with gene-specific parental alt (PA) frequency, after removal of SB filtered sites

ddg2p_path = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2014-11-04/DDG2P/DDG2P_freeze_with_gencode19_genomic_coordinates_20141118_fixed.txt"
de_novos_path = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/denovo_gear_trios_extracted_passed_variants_18.12.14.tsv"

# estimated error rate at 0.0012 from DNM calls in parents in DDG2P genes, set 
# slightly higher to be conservative
error.rate = 0.002 

# threshold for removing sites with high strand bias, or parental alt frequency
p_cutoff = 1e-3 

#' run some preliminary filtering of de novos
#' 
#' We want to filter out de novos with high MAF, where they are not present in
#' the child VCF, or are present in the parental VCFs
preliminary_filtering <- function(dnms) {
    # annotate with numeric max allele freq
    # fix the blank and null max AF values
    dnms$max_af[dnms$max_af == "" | dnms$max_af == "."] = 0
    
    # one de novo has a comma separated list of MAF values (both above 0.1)
    dnms$max_af[grepl(",", dnms$max_af)] = NA
    dnms$max_af = as.numeric(dnms$max_af)
    
    # keep sites in child vcf not in parental vcfs (leaves 123919)
    dnms = dnms[dnms$in_child_vcf == 1 & dnms$in_father_vcf == 0 & dnms$in_mother_vcf == 0, ]
    
    # filter on max_af (leaves 23490)
    dnms = dnms[dnms$max_af < 0.01 & !is.na(dnms$max_af), ]
    
    # remove dnms in samples with >> too many DNMs, focus on too many DNMs at 
    # high quality (leaves 17794)
    # NEED TO UPDATE WITH CAROLINE'S NEW SAMPLE FILE LIST, OR USE PRE-FILTERED SET OF DNMS
    sample.fails = c("276227", "258876", "273778", "258921", "272110", "260337",
     "264083", "264084", "269602", "265624")
    dnms = dnms[!dnms$decipher_id %in% sample.fails, ]
    
    # Annotate dnms hitting coding exons or splice sites
    coding_splicing = c("frameshift_variant", "inframe_deletion", 
        "inframe_insertion", "initiator_codon_variant", "missense_variant",  
        "splice_acceptor_variant", "splice_donor_variant", "stop_gained", 
        "stop_lost", "synonymous_variant")
    dnms$coding = dnms$consequence %in% coding_splicing
    
    # filtering down to coding (leaves 9057)
    # dnms = dnms[dnms$coding, ]
    
    # filter to high quality de novos (leaves 5536)
    # dnms = dnms[dnms$pp_dnm > 0.9, ]
    
    return(dnms)
}


#' annotate with presence in DNG monoallelic/XL genes
#' NEED TO UPDATE FOR MOST RECENT DDG2P
annotate_with_ddg2p <- function(dnms) {
    ddg2p = read.table(ddg2p_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # find the dominantly inherited genes
    allelic = c("Monoallelic", "Hemizygous", "X-linked dominant", "Both")
    mono.xl.ddg2p = ddg2p[ddg2p$Allelic_requirement %in% allelic, ]
    mono.xl.dd.genes = unique(mono.xl.ddg2p$GeneSymbol)
    
    dnms$dnms.mono.xl.ddg2p = dnms$symbol %in% mono.xl.dd.genes
    
    return(dnms)
}

#' extract ALT and REF counts of forward and reverse reads for the members of 
#' each trio
extract_alt_and_ref_counts <- function(dnms) {
    
    dnms$dp4.child = strsplit(dnms$dp4_child, split=",")
    
    dnms$child.REF.F = as.numeric(sapply(dnms$dp4.child, "[", 1))
    dnms$child.REF.R = as.numeric(sapply(dnms$dp4.child, "[", 2))
    dnms$child.ALT.F = as.numeric(sapply(dnms$dp4.child, "[", 3))
    dnms$child.ALT.R = as.numeric(sapply(dnms$dp4.child, "[", 4))
    
    dnms$dp4.father = strsplit(as.character(dnms$dp4_father), split=",")
    
    dnms$father.REF.F = as.numeric(sapply(dnms$dp4.father, "[", 1))
    dnms$father.REF.R = as.numeric(sapply(dnms$dp4.father, "[", 2))
    dnms$father.ALT.F = as.numeric(sapply(dnms$dp4.father, "[", 3))
    dnms$father.ALT.R = as.numeric(sapply(dnms$dp4.father, "[", 4))
    
    dnms$dp4.mother = strsplit(as.character(dnms$dp4_mother), split=",")
    
    dnms$mother.REF.F = as.numeric(sapply(dnms$dp4.mother, "[", 1))
    dnms$mother.REF.R = as.numeric(sapply(dnms$dp4.mother, "[", 2))
    dnms$mother.ALT.F = as.numeric(sapply(dnms$dp4.mother, "[", 3))
    dnms$mother.ALT.R = as.numeric(sapply(dnms$dp4.mother, "[", 4))
    
    dnms$count.child.alt = dnms$child.ALT.F + dnms$child.ALT.R
    
    # remove the DP4 columns (since they are now lists)
    dnms$dp4.child = NULL
    dnms$dp4.mother = NULL
    dnms$dp4.father = NULL
    
    return(dnms)
}

#'' count of number of times that site is called
count_site_and_gene_recurrence <- function(dnms) {
    
    dnms$key = paste(dnms$chrom, dnms$pos, dnms$alt, sep="_")
    table.key = table(dnms$key)
    dnms$count.sites = table.key[match(dnms$key, row.names(table.key))]
    
    # count of number times that gene is called
    table.genes = table(dnms$symbol)
    dnms$count.genes = table.genes[match(dnms$symbol, row.names(table.genes))]
    
    return(dnms)
}

#' counts REF and ALT alleles in reads for de novo events at a single site
get_allele_counts <- function(site_vars) {
    
    REF.F = sum(site_vars[, c("child.REF.F", "mother.REF.F", "father.REF.F")])
    REF.R = sum(site_vars[, c("child.REF.R", "mother.REF.R", "father.REF.R")])
    ALT.F = sum(site_vars[, c("child.ALT.F", "mother.ALT.F", "father.ALT.F")])
    ALT.R = sum(site_vars[, c("child.ALT.R", "mother.ALT.R", "father.ALT.R")])
    
    # count parental alt and ref
    parent.ALT = sum(site_vars[, c("mother.ALT.F", "father.ALT.F", "mother.ALT.R", "father.ALT.R")])
    parent.REF = sum(site_vars[, c("mother.REF.F", "father.REF.F", "mother.REF.R", "father.REF.R")])
    
    return(list(REF.F=REF.F, REF.R=REF.R, ALT.F=ALT.F, ALT.R=ALT.R, 
        parent.ALT=parent.ALT, parent.REF=parent.REF))
}

#' counts REF and ALT alleles in reads for de novo events within a single gene
get_parental_counts <- function(gene_vars) {
    
    parent.ALT = sum(gene_vars[, c("mother.ALT.F", "father.ALT.F", "mother.ALT.R", "father.ALT.R")])
    parent.REF = sum(gene_vars[, c("mother.REF.F", "father.REF.F", "mother.REF.R", "father.REF.R")])
    
    return(list(gene.ALT=parent.ALT, gene.REF=parent.REF))
}

#' checks for strand bias per de novo site using the ref and alt counts
site_strand_bias <- function(site) {
    x = c(site[["REF.F"]], site[["REF.R"]], site[["ALT.F"]], site[["ALT.R"]])
    x = matrix(x, nrow=2)
    return(fisher.test(x)$p.value)
}

#' tests each site for deviation from expected behaviour
test_sites <- function(dnms) {
    alleles = subset(dnms, select=c("key", 
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
    PA_pval_site = apply(parent_counts, 1, binom.test, p=error.rate, alternative="greater")
    results$PA_pval_site = as.numeric(sapply(PA_pval_site, "[", "p.value"))
    
    # check for strand bias by fishers exact test on the allele counts
    cat("testing strand bias\n")
    allele_counts = lapply(split(results, seq_along(results[, 1])), as.list)
    results$SB_pval = sapply(allele_counts, site_strand_bias)
    
    return(results)
}

#' tests each gene for deviation from expected behaviour
test_genes <-function(dnms) {
    # exclude de novo SNVs that fail the strand bias filter, otherwise these 
    # skew the parental alts within genes
    alleles = dnms[!(dnms$SB_pval < p_cutoff & dnms$var_type == "DENOVO-SNP"), ]
    
    # loop to calculate gene-specific PA values after SB filtering
    alleles = subset(alleles, select=c("key", "symbol", 
        "mother.REF.F", "mother.REF.R", "mother.ALT.F", "mother.ALT.R",
        "father.REF.F", "father.REF.R", "father.ALT.F", "father.ALT.R"))
    
    # count the number of parental alleles within genes, and restructure data
    # for testing
    cat("counting alleles per gene\n")
    genes = split(alleles, alleles$symbol)
    counts = lapply(genes, get_parental_counts)
    results = data.frame(matrix(unlist(counts), ncol=length(counts[[1]]), byrow=TRUE))
    names(results) = names(counts[[1]])
    results$symbol = names(counts)
    
    # check for overabundance of parental alt alleles using binomial test
    cat("testing parental alt overabundance\n")
    parent_counts = data.frame(results$gene.ALT, (results$gene.ALT + results$gene.REF))
    PA_pval = apply(parent_counts, 1, binom.test, p=error.rate, alternative="greater")
    results$PA_pval_gene = as.numeric(sapply(PA_pval, "[", "p.value"))
    
    return(results)
}

#' set flags for filtering, fail samples with strand bias < threshold, or any 2 of 
#'  (i) both parents have ALTs 
#'  (ii) site-specific parental alts < threshold, 
#'  (iii) gene-specific parental alts < threshold, if > 1 sites called in gene
set_filter_flag <- function(de_novos) {
    
    de_novos$overall.pass = "PASS"
    
    # fail SNVs with excessive strand bias
    de_novos$overall.pass[de_novos$SB_pval < p_cutoff & de_novos$var_type == "DENOVO-SNP"] = "FAIL"
    
    # fail sites with gene-specific parental alts, only if >1 sites called per gene
    # NEED TO ADAPT TO AVOID SITES WITHOUT GENE SYMBOL ANNOTATED BEING REGARDED AS BEING IN THE SAME GENE
    gene_fail = de_novos$PA_pval_gene < p_cutoff & de_novos$count.genes > 1
    site_fail = de_novos$PA_pval_site < p_cutoff
    excess_alts = de_novos$min.parent.ALT > 0
    
    # fail sites with parental alts, any two of three classes
    parental_alts = data.frame(gene_fail, site_fail, excess_alts)
    parental_alts_fail = apply(parental_alts, 1, sum) >= 2
    
    # fail sites with excess parental alt reads
    de_novos$overall.pass[parental_alts_fail] = "FAIL"
    
    return(de_novos)
}


main <- function() {
    de_novos = read.table(de_novos_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    de_novos = preliminary_filtering(de_novos)
    de_novos = annotate_with_ddg2p(de_novos)
    de_novos = extract_alt_and_ref_counts(de_novos)
    de_novos = count_site_and_gene_recurrence(de_novos)
    
    site_results = test_sites(de_novos)
    de_novos = merge(de_novos, site_results, by="key", all.x=TRUE)
    
    # get the minimum alternate allele count from the parents
    alts = data.frame(de_novos$mother.ALT.F + de_novos$mother.ALT.R, de_novos$father.ALT.F + de_novos$father.ALT.R)
    de_novos$min.parent.ALT = apply(alts, 1, min)
    
    gene_results = test_genes(de_novos)
    de_novos = merge(de_novos, gene_results, by="symbol", all.x=TRUE)
    
    de_novos = set_filter_flag(de_novos)
    
    a = de_novos[de_novos$overall.pass == "PASS" & de_novos$coding == TRUE, ]
    
    write.table(de_novos, file="~/de_novos.all.2.txt", quote=F, row.names=F, sep="\t")
}


main()

