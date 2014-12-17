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


# set PASS/FAIL status for research analyses of gene-specific enrichment


dnms<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Calling_Pipeline/DNG_0.5.4/y1_y2000_y2_y3_4344_trios_extracted_passed_variants_full_21.11.14.tsv", header=T)

error.rate<-0.002 # estimated at 0.0012 from DNM calls in parents in DDg2P genes, set slightly higher to be conservative

SB.filt<-1e-3 # threshold for removing sites with too high strand bias, or Parental Alt Frequency


# annotate with numeric max allele freq

max_freq<-as.numeric(as.character(dnms$max_af))

max_freq[is.na(max_freq)]<-0

dnms<-cbind(dnms, max_freq)



# keep sites in child vcf not in parental vcfs

dnms.vcf<-dnms[which((dnms$in_child_vcf-dnms$in_father_vcf-dnms$in_mother_vcf)==1),]

# 123919 remaining


# filter on max_af

dnms.rare<-dnms.vcf[which(dnms.vcf$max_freq<0.01),]

# 23490 remaining



# remove dnms in samples with >> too many DNMs, focus on too many DNMs at high quality
# NEED TO UPDATE WITH CAROLINE'S NEW SAMPLE FILE LIST, OR USE PRE-FILTERED SET OF DNMS


sample.fails<-c("276227", "258876", "273778", "258921", "272110", "260337", "264083", "264084", "269602", "265624")



dnms.rare.qc<-dnms.rare[-which((dnms.rare$decipher_id %in% sample.fails)=="TRUE"),]

# 17794 remaining



# Annotate dnms hitting coding exons or splice sites

coding_splicing<-c("frameshift_variant", "inframe_deletion", "inframe_insertion", "initiator_codon_variant", "missense_variant",  "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost", "synonymous_variant")

dnms.coding<-cbind(dnms.rare.qc, dnms.rare.qc$consequence %in% coding_splicing)

#dnms.rare.qc[which((dnms.rare.qc$consequence %in% coding_splicing)=="TRUE"),]
# 9057 remaining


#dnms.coding.high<-dnms.coding[which(dnms.coding$pp_dnm>0.9),]
# 5536 remaining



# annotate with presence in DNG monoallelic/XL genes
# NEED TO UPDATE FOR MOST RECENT DDG2P

ddg2p<-read.delim("~/Desktop/DDD/DDG2P/ddg2p_with_gene_positions_and_ids_20140718.txt", header=T, sep="|")

allelic<-c("Monoallelic", "Hemizygous", "X-linked dominant", "Both")

allelic.index<-ddg2p$Allelic_requirement %in% allelic 

mono.xl.ddg2p<-ddg2p[which(allelic.index=="TRUE"),]

mono.xl.dd.genes<-unique(mono.xl.ddg2p$GeneSymbol)

dnms.mono.xl.ddg2p<-dnms.coding$symbol %in% mono.xl.dd.genes

dnms.coding<-cbind(dnms.coding, dnms.mono.xl.ddg2p)


# extract ALT and REF counts of F and R reads

dng<-dnms.coding

dp4.child<-t(as.data.frame(strsplit(as.character(dng$dp4_child), split=",")))

dp4.child.REF.F<-as.numeric(dp4.child[,1])
dp4.child.REF.R<-as.numeric(dp4.child[,2])
dp4.child.ALT.F<-as.numeric(dp4.child[,3])
dp4.child.ALT.R<-as.numeric(dp4.child[,4])

dp4.father <-t(as.data.frame(strsplit(as.character(dng$dp4_father), split=",")))

dp4.father.REF.F<-as.numeric(dp4.father[,1])
dp4.father.REF.R<-as.numeric(dp4.father[,2])
dp4.father.ALT.F<-as.numeric(dp4.father[,3])
dp4.father.ALT.R<-as.numeric(dp4.father[,4])

dp4.mother<-t(as.data.frame(strsplit(as.character(dng$dp4_mother), split=",")))

dp4.mother.REF.F<-as.numeric(dp4.mother[,1])
dp4.mother.REF.R<-as.numeric(dp4.mother[,2])
dp4.mother.ALT.F<-as.numeric(dp4.mother[,3])
dp4.mother.ALT.R<-as.numeric(dp4.mother[,4])

count.child.alt<-dp4.child.ALT.F + dp4.child.ALT.R

dnms.coding<-cbind(dnms.coding, dp4.child.REF.F, dp4.child.REF.R, dp4.child.ALT.F, dp4.child.ALT.R, dp4.father.REF.F, dp4.father.REF.R, dp4.father.ALT.F, dp4.father.ALT.R, dp4.mother.REF.F, dp4.mother.REF.R, dp4.mother.ALT.F, dp4.mother.ALT.R, count.child.alt)


# add count of number of times that site is called

key<-paste(dnms.coding$chrom, dnms.coding$pos, dnms.coding$alt, sep="_")

table.key<-table(key)

count.sites<-table.key[match(key, row.names(table.key))]

dnms.coding<-cbind(dnms.coding, key, count.sites)


# add count of number times that gene is called

table.genes<-table(dnms.coding$symbol)

count.genes<-table.genes[match(dnms.coding$symbol, row.names(table.genes))]

dnms.coding<-cbind(dnms.coding, count.genes)


# vectors to store information from loop, code 2 = not assessed

dng<-dnms.coding

num.sites<-length(dng[,1])


SB_pval<-rep(2,length(dng[,1]))
PA_pval<-rep(2,length(dng[,1]))
PA_pval_gene<-rep(2,length(dng[,1]))

count.REF.F<-rep(2,length(dng[,1]))
count.REF.R<-rep(2,length(dng[,1]))
count.ALT.F<-rep(2,length(dng[,1]))
count.ALT.R<-rep(2,length(dng[,1]))
count.parent.ALT<-rep(2,length(dng[,1]))
count.parent.REF<-rep(2,length(dng[,1]))
min.parent.ALT<-rep(2,length(dng[,1]))



# loop to calculate site-specific PA and SB values

#for (i in seq(1,10)) { #for testing

for (i in seq(1,num.sites)) {
	
	# annotate sites where both parents have ALT reads
	
	min.parent.ALT[i]<-min(dng$dp4.mother.ALT.F[i]+dng$dp4.mother.ALT.R[i], dng$dp4.father.ALT.F[i]+dng$dp4.father.ALT.R[i])
	
	
	# calculcate site-specific SB p value
	
	dng.site<-dng[which(dng$key==dng$key[i]),]
		
	count.REF.F[i]<-sum(dng.site$dp4.child.REF.F)+sum(dng.site$dp4.mother.REF.F)+sum(dng.site$dp4.father.REF.F)
	count.REF.R[i]<-sum(dng.site$dp4.child.REF.R)+sum(dng.site$dp4.mother.REF.R)+sum(dng.site$dp4.father.REF.R)
	count.ALT.F[i]<-sum(dng.site$dp4.child.ALT.F)+sum(dng.site$dp4.mother.ALT.F)+sum(dng.site$dp4.father.ALT.F)
	count.ALT.R[i]<-sum(dng.site$dp4.child.ALT.R)+sum(dng.site$dp4.mother.ALT.R)+sum(dng.site$dp4.father.ALT.R)
					
	SB_pval[i]<-fisher.test(matrix(c(count.REF.F[i], count.REF.R[i], count.ALT.F[i], count.ALT.R[i]), nrow=2))$p.value
	 
	
	# calculate site-specific PA p value
				
	count.parent.ALT[i]<-sum(dng.site$dp4.mother.ALT.F)+sum(dng.site$dp4.father.ALT.F)+sum(dng.site$dp4.mother.ALT.R)+sum(dng.site$dp4.father.ALT.R)
			
	count.parent.REF[i]<-sum(dng.site$dp4.mother.REF.F)+sum(dng.site$dp4.father.REF.F)+sum(dng.site$dp4.mother.REF.R)+sum(dng.site$dp4.father.REF.R)
	
	PA_pval[i]<-binom.test(as.numeric(count.parent.ALT[i]), as.numeric(count.parent.ALT[i])+as.numeric(count.parent.REF[i]), error.rate, alternative="g")$p.value
		
#	print(i)
		
}


# annotate sites with site-specific SB and PA p values

dnms.coding.pvals<-cbind(dnms.coding,count.REF.F, count.REF.R, count.ALT.F, count.ALT.R, count.parent.ALT, count.parent.REF, min.parent.ALT, SB_pval, PA_pval)



# loop to calculate gene-specific PA values after SB filtering

dng<-dnms.coding.pvals

#for (i in seq(1,10)) { #for testing

for (i in seq(1,num.sites)) {
	
		# calculate PA p value
			
		dng.gene<-dng[which(dng$symbol==dng$symbol[i]),] # collate all DNMs called in the same gene
		
		remove.index<-which(dnms.coding$SB_pval<SB.filt & dnms.coding$var_type=="DENOVO-SNP") # index of sites to remove
			
		if(length(remove.index>0)) dng.gene<-dng.gene[-which(dnms.coding$SB_pval<SB.filt & dnms.coding$var_type=="DENOVO-SNP"),] # remove SNVs failing the SB filter, if any are present
		
		if(length(dng.gene[,1])>0) {
			
			count.parent.ALT[i]<-sum(dng.gene$dp4.mother.ALT.F)+sum(dng.gene$dp4.father.ALT.F)+sum(dng.gene$dp4.mother.ALT.R)+sum(dng.gene$dp4.father.ALT.R) # count number of parental REF alleles across sites
			
			count.parent.REF[i]<-sum(dng.gene$dp4.mother.REF.F)+sum(dng.gene$dp4.father.REF.F)+sum(dng.gene$dp4.mother.REF.R)+sum(dng.gene$dp4.father.REF.R) # count number of parental ALT alleles across sites
	
			PA_pval_gene[i]<-binom.test(as.numeric(count.parent.ALT[i]), as.numeric(count.parent.ALT[i])+as.numeric(count.parent.REF[i]), error.rate, alternative="g")$p.value # calculate p value assuming an error rate of seeing the observed number of ALT reads

			}

#	print(i)
		
}


dnms.coding.all<-cbind(dnms.coding.pvals, PA_pval_gene)


# set flags for filtering, fail samples with SB < threshold, or any 2 of (i) both parents have ALTs (ii) site-specific PA < threshold, (iii) gene-specific PA < threshold, if >1 sites called in gene

overall.pass<-rep("PASS",length(dnms.coding.all[,1])) # store overall PASS status
PA.gene.pass<-rep("PASS",length(dnms.coding.all[,1])) # store intermediate PA gene PASS status
PA.pass<-rep("PASS",length(dnms.coding.all[,1])) # store intermediate PA PASS status

# fail SNVs with SB
overall.pass[which(dnms.coding.all$SB_pval<SB.filt & dnms.coding.all$var_type=="DENOVO-SNP")]<-"FAIL"

# fail sites with gene-specific PA, only if >1 sites called per gene
# NEED TO ADAPT TO AVOID SITES WITHOUT GENE SYMBOL ANNOTATED BEING REGARDED AS BEING IN THE SAME GENE
PA.gene.pass[which(dnms.coding.all$PA_pval_gene<SB.filt & dnms.coding.all$count.genes>1)]<-"FAIL"

# fail sites with PA, any two of three classes
PA.pass[which(PA.gene.pass=="FAIL" & dnms.coding.all$PA_pval<SB.filt)]<-"FAIL"
PA.pass[which(PA.gene.pass=="FAIL" & dnms.coding.all$min.parent.ALT>0)]<-"FAIL"
PA.pass[which(dnms.coding.all$PA_pval<SB.filt & dnms.coding.all$min.parent.ALT>0)]<-"FAIL"

# fail sites with PA
overall.pass[which(PA.pass=="FAIL")]<-"FAIL"

#table(overall.pass)

#table(overall.pass, dnms.coding.all[,41])


dnms.coding.all<-cbind(dnms.coding.all, overall.pass)


# write out

write.table(dnms.coding.all, file="~/Desktop/dnms.coding.all.txt", quote=F, row.names=F, sep="\t")


