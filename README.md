### denovoFilter

sites such that gene-specific DNM enrichment analyses are not appreciably
contaminated with dodgy genes/sites. Further specificity will be achieved
through downstream experimental validation of sites in putative ‘hit’ genes.

After filtering there are about 1.7 candidate DNMs per trio.

The filtering is based on a combination of site-specific strand bias (SB) and
3 parental ALT (PA) metrics. A site fails if it fails either the strand bias
filter, or any two of the three parental ALT metrics. These rules have been set
so as to allow us to pass:
 * good sites in dodgy genes,
 * sites that are mosaic in one parent.
 * high depth sites that have a very small number of ALT reads in both parents
   by chance (typically 1 ALT read in each parent).

Reject site:
  If SB p value < 1e-3 (SNVs only)

OR

If any 2 of the following conditions are met:
 1. ALT reads present in both parents
 2. site-specific PA p value <1e-3
 3. gene-wide PA p value < 1e-3 (only if >1 sites called per gene after SB
    filter applied, otherwise it is redundant with the site specific PA p value)
    
### Install
Clone the repository and install with:

```sh
git clone https://github.com/jeremymcrae/denovoFilter
cd denovoFilter
python setup.py install --user
```

### Filtering *de novos*
Define the paths:
 * candidate de novos listed in tab-separated file containing 23 columns for person_stable_id, gender, mother_stable_id, father_stable_id, chrom, pos, ref, alt, symbol, dp4_child, dp4_mother, dp4_father, pp_dnm, max_af, child_alt_prp, consequence, in_child_vcf, in_mother_vcf, in_father_vcf,  decipher_id, var_type, ensg, enst
