 The rationale is to prioritise sensitivity while removing sufficient dodgy 
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