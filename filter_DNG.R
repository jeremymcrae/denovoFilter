# filter DNG callset initially on basis of

library(mupit)

# - in_child_vcf, not in parent vcf
# - max_af < 0.01
# - remove dodgy samples with too many DNM calls

# annotate DNG callset

# - with coding/splicing
# - with presence in DDG2P
# - with strand-specitic counts of REF and ALT reads
# - with site-specific strand bias (SB) and parental alt (PA) frequency p values against null
# - with gene-specific parental alt (PA) frequency, after removal of SB filtered sites
DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"
DDG2P_PATH = file.path(DATAFREEZE_DIR, "DDG2P_freeze_with_gencode19_genomic_coordinates_20141118_fixed.txt")
DE_NOVOS_PATH = file.path(DATAFREEZE_DIR, "denovo_gear_trios_extracted_passed_variants_17.04.15.tsv")

# estimated error rate at 0.0012 from DNM calls in parents in DDG2P genes, set
# slightly higher to be conservative
ERROR_RATE = 0.002

# threshold for removing sites with high strand bias, or parental alt frequency
P_CUTOFF = 1e-3

#' fix a few MAF values that are clank, or have lists of MAF values
fix_maf <- function(de_novos) {
    # annotate with numeric max allele freq
    # fix the blank and null max AF values
    de_novos$max_af[de_novos$max_af == "" | de_novos$max_af == "."] = 0
    
    # one de novo has a comma separated list of MAF values (both above 0.1)
    de_novos$max_af[grepl(",", de_novos$max_af)] = NA
    de_novos$max_af = as.numeric(de_novos$max_af)
    
    return(de_novos)
}

#' run some preliminary filtering of de novos
#'
#' We want to filter out de novos with high MAF, where they are not present in
#' the child VCF, or are present in the parental VCFs
preliminary_filtering <- function(de_novos) {
    de_novos = fix_maf(de_novos)
    
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
    coding_splicing = c("frameshift_variant", "inframe_deletion",
        "inframe_insertion", "initiator_codon_variant", "missense_variant",
        "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
        "stop_gained", "stop_lost", "synonymous_variant")
        
    de_novos$coding = de_novos$consequence %in% coding_splicing
    
    return(de_novos)
}

#' adds gene symbols to variants lacking them, using a mupit function
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


#' annotate with presence in DNG monoallelic/XL genes
annotate_with_ddg2p <- function(de_novos) {
    
    ddg2p = read.table(DDG2P_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # find the dominantly inherited genes
    allelic = c("Monoallelic", "Hemizygous", "X-linked dominant", "Both")
    dominant_genes = ddg2p[ddg2p$Allelic_requirement %in% allelic, ]
    dominant_symbols = unique(dominant_genes$ddg2p_gene_name)
    
    de_novos$in_dominant_ddg2p = de_novos$symbol %in% dominant_symbols
    
    return(de_novos)
}

#' extract ALT and REF counts of forward and reverse reads for the members of
#' each trio
extract_alt_and_ref_counts <- function(de_novos) {
    
    de_novos$dp4.child = strsplit(de_novos$dp4_child, split=",")
    
    de_novos$child.REF.F = as.numeric(sapply(de_novos$dp4.child, "[", 1))
    de_novos$child.REF.R = as.numeric(sapply(de_novos$dp4.child, "[", 2))
    de_novos$child.ALT.F = as.numeric(sapply(de_novos$dp4.child, "[", 3))
    de_novos$child.ALT.R = as.numeric(sapply(de_novos$dp4.child, "[", 4))
    
    de_novos$dp4.father = strsplit(as.character(de_novos$dp4_father), split=",")
    
    de_novos$father.REF.F = as.numeric(sapply(de_novos$dp4.father, "[", 1))
    de_novos$father.REF.R = as.numeric(sapply(de_novos$dp4.father, "[", 2))
    de_novos$father.ALT.F = as.numeric(sapply(de_novos$dp4.father, "[", 3))
    de_novos$father.ALT.R = as.numeric(sapply(de_novos$dp4.father, "[", 4))
    
    de_novos$dp4.mother = strsplit(as.character(de_novos$dp4_mother), split=",")
    
    de_novos$mother.REF.F = as.numeric(sapply(de_novos$dp4.mother, "[", 1))
    de_novos$mother.REF.R = as.numeric(sapply(de_novos$dp4.mother, "[", 2))
    de_novos$mother.ALT.F = as.numeric(sapply(de_novos$dp4.mother, "[", 3))
    de_novos$mother.ALT.R = as.numeric(sapply(de_novos$dp4.mother, "[", 4))
    
    de_novos$count.child.alt = de_novos$child.ALT.F + de_novos$child.ALT.R
    
    # remove the DP4 columns (since they are now lists)
    de_novos$dp4.child = NULL
    de_novos$dp4.mother = NULL
    de_novos$dp4.father = NULL
    
    return(de_novos)
}

#' count of number of times that site is called
count_site_and_gene_recurrence <- function(de_novos) {
    
    de_novos$key = paste(de_novos$chrom, de_novos$pos, de_novos$alt, sep="_")
    table.key = table(de_novos$key)
    de_novos$count.sites = table.key[match(de_novos$key, row.names(table.key))]
    
    # count of number times that gene is called
    table.genes = table(de_novos$symbol)
    de_novos$count.genes = table.genes[match(de_novos$symbol, row.names(table.genes))]
    
    return(de_novos)
}

#' counts REF and ALT alleles in reads for de novo events at a single site
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
site_strand_bias <- function(site) {
    x = c(site[["REF.F"]], site[["REF.R"]], site[["ALT.F"]], site[["ALT.R"]])
    x = matrix(x, nrow=2)
    return(fisher.test(x)$p.value)
}

#' tests each site for deviation from expected behaviour
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
    
    return(results)
}

#' tests each gene for deviation from expected behaviour
test_genes <-function(de_novos) {
    
    stopifnot("PA_pval_site" %in% names(de_novos))
    
    # exclude de novo SNVs that fail the strand bias filter, otherwise these
    # skew the parental alts within genes
    alleles = de_novos[!(de_novos$SB_pval < P_CUTOFF & de_novos$var_type == "DENOVO-SNP"), ]
    
    # loop to calculate gene-specific PA values after SB filtering
    alleles = subset(alleles, select=c("key", "symbol",
        "mother.REF.F", "mother.REF.R", "mother.ALT.F", "mother.ALT.R",
        "father.REF.F", "father.REF.R", "father.ALT.F", "father.ALT.R"))
    
    # count the number of parental alleles within genes, and restructure data
    # for testing
    cat("counting alleles per gene\n")
    genes = split(alleles, alleles$symbol)
    counts = lapply(genes, get_allele_counts, gene=TRUE)
    results = data.frame(matrix(unlist(counts), ncol=length(counts[[1]]), byrow=TRUE))
    names(results) = names(counts[[1]])
    results$symbol = names(counts)
    
    # check for overabundance of parental alt alleles using binomial test
    cat("testing parental alt overabundance\n")
    parent_counts = data.frame(results$gene.ALT, (results$gene.ALT + results$gene.REF))
    PA_pval = apply(parent_counts, 1, binom.test, p=ERROR_RATE, alternative="greater")
    results$PA_pval_gene = as.numeric(sapply(PA_pval, "[", "p.value"))
    
    return(results)
}

#' set flags for filtering, fail samples with strand bias < threshold, or any 2 of
#'  (i) both parents have ALTs
#'  (ii) site-specific parental alts < threshold,
#'  (iii) gene-specific parental alts < threshold, if > 1 sites called in gene
set_filter_flag <- function(de_novos, keep_all=FALSE) {
    
    de_novos$overall.pass = "PASS"
    
    # fail SNVs with excessive strand bias
    de_novos$overall.pass[de_novos$SB_pval < P_CUTOFF & de_novos$var_type == "DENOVO-SNP"] = "FAIL"
    
    # find if each de novo has passed each of three different filtering strategies
    # fail sites with gene-specific parental alts, only if >1 sites called per gene
    gene_fail = de_novos$PA_pval_gene < P_CUTOFF & de_novos$count.genes > 1
    site_fail = de_novos$PA_pval_site < P_CUTOFF
    
    # fail sites
    excess_alts = de_novos$min.parent.ALT > 0
    
    # fail sites with parental alts, any two of three classes
    sites = data.frame(gene_fail, site_fail, excess_alts)
    sites_fail = apply(sites, 1, sum) >= 2
    
    # fail sites with excess parental alt reads
    de_novos$overall.pass[sites_fail] = "FAIL"
    
    if (!keep_all) {
        # subset down to specific columns
        de_novos = subset(de_novos, select=c("person_stable_id", "gender",
            "mother_stable_id", "father_stable_id", "chrom", "pos", "ref", "alt",
            "symbol", "var_type", "consequence", "ensg", "enst", "max_af", "pp_dnm",
            "child_alt_prp", "in_dominant_ddg2p", "coding", "overall.pass",
            "child.REF.F", "child.REF.R", "child.ALT.F", "child.ALT.R",
            "mother.REF.F", "mother.REF.R", "mother.ALT.F", "mother.ALT.R",
            "father.REF.F", "father.REF.R", "father.ALT.F", "father.ALT.R"))
    }
    
    return(de_novos)
}


main <- function() {
    de_novos = read.table(DE_NOVOS_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    de_novos = preliminary_filtering(de_novos)
    de_novos = fix_missing_gene_symbols(de_novos)
    de_novos = annotate_with_ddg2p(de_novos)
    de_novos = extract_alt_and_ref_counts(de_novos)
    de_novos = count_site_and_gene_recurrence(de_novos)
    
    site_results = test_sites(de_novos)
    de_novos = merge(de_novos, site_results, by="key", all.x=TRUE)
    
    # get the minimum alternate allele count from the parents
    alts = data.frame(de_novos$mother.ALT.F + de_novos$mother.ALT.R,
        de_novos$father.ALT.F + de_novos$father.ALT.R)
    de_novos$min.parent.ALT = apply(alts, 1, min)
    
    # test whether any genes have an excess of parental alts
    gene_results = test_genes(de_novos)
    de_novos = merge(de_novos, gene_results, by="symbol", all.x=TRUE)
    
    de_novos = set_filter_flag(de_novos)
    passed = de_novos[de_novos$overall.pass == "PASS" & de_novos$coding, ]
    passed = get_independent_de_novos(passed)
    
    write.table(passed, file="~/de_novos.ddd_4k.ddd_only.txt",
        quote=FALSE, row.names=FALSE, sep="\t")
}

# R equivalent of "if main"
if(getOption("run.main", default=TRUE)) {
    main()
}
