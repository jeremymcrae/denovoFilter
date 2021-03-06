### denovoFilter

Filter sites such that gene-specific de novo enrichment analyses are not
contaminated with dodgy genes/sites. Further specificity will be achieved
through downstream experimental validation of sites in putative 'hit' genes.

The filtering is based on a combination of site-specific strand bias and
three parental alt metrics. A site fails if it fails either the strand bias
filter, or any two of the three parental alt metrics. These rules have been set
so as to allow us to pass:
 * good sites in dodgy genes
 * sites that are mosaic in one parent
 * high depth sites that have a very small number of alt reads in both parents
   by chance (typically 1 alt read in each parent)

Reject site:
  If strand-bias p-value < 1e-3 (SNVs only)

OR

If any two of the following conditions are met:
 1. alt reads present in both parents
 2. site-specific parental alt p-value <1e-3
 3. gene-wide parental alt p-value < 1e-3 (only if >1 sites called per gene after
    strand bias filter applied, otherwise it is redundant with the site specific
    parental alt p-value)
    
### Install
Clone the repository and install with:

```sh
git clone https://github.com/jeremymcrae/denovoFilter.git
cd denovoFilter
python setup.py install --user
```

### Filtering the candidate *de novos*
A script is included (`scripts/filter_de_novos.py`) which can be used as:
```sh
python scripts/filter_de_novos.py \
  --de-novos DE_NOVO_PATH \
  --families FAMILIES_PATH \
  --output OUTPUT_PATH
```

You can also use other optional flags:
 * `--de-novos-indels MISSING_INDELS_PATH` to filter a set of candidate indels
   that were missed by denovogear.
 * `--sample-fails SAMPLE_FAILS_PATH` to exclude a set of probands who we know
   have problems with de novo calling.
 * `--sample-fails-indels SAMPLE_FAILS_INDELS_PATH` to exclude a set of probands
   from the missing indels candidates, due to those problems showing problems
   in the missing indels dataset.
 * `--annotate-only` to add an extra column ('pass') which has True/False values
   for whether the variants pass rather than filtering to a smaller subset. By
   default the script will exclude site which fail the filtering.
 * `--include-noncoding` to include noncoding sites in the filtered output.
 * `--include-recurrent` to skip screening out sites which occur multiple times
   in a family, or multiple sites in a gene in a single individual.

### Input files
#### Definitions for the required columns in the candidate *de novos* file
| name             | example       | definition                            |
| -----------      | ------------- | -------------                         |
| person_stable_id | "DDDP100001"  | ID of the proband                     |
| chrom            | "1"           | chromosome that the candidate is on   |
| pos              | 1000001       | nucleotide position of the candidate along the chromosome |
| ref              | "A"           | base-code for the reference allele    |
| alt              | "G"           | base-code for the alternate allele    |
| symbol           | "ATRX"        | HGNC symbol of the gene that the variant lies in (or blank for no gene) |
| dp4_child        | 10,12,9,10    | Comma-separated read-depths in the child for the reference allele in the forward orientation, the reference allele in the reverse orientation, the alternate allele in the forward orientation, and the alternate allele in the reverse orientation. Orientation is with respect to the reference genome. |
| dp4_mother       | 15,17,0,0     | Read-depths in the mother as for dp4_child |
| dp4_father       | 13,14,0,0     | Read-depths in the father as for dp4_child |
| pp_dnm           | 0.952         | Posterior probability of de novo mutation from denovogear |
| max_af           | 0.0001        | Highest allele frequency for alternate allele from reference healthy control population |
| consequence      | "missense_variant" | VEP-annotated consequence for the variant |
| in_child_vcf     | 1             | Whether the variant is present in the child's VCF (1=true, 0=false) |
| in_mother_vcf    | 0             | Whether the variant is present in the mother's VCF (1=true, 0=false) |
| in_father_vcf    | 0             | Whether the variant is present in the father's VCF (1=true, 0=false) |

#### Definitions for the required columns in the family relationships file
| name          | example       | definition       |
| -----------   | ------------- | -----            |
| family_id     | "FAM100001"   | ID of the family |
| individual_id | "DDDP100001"  | ID of the proband, matching the person_stable_id of the de novos table |
| sex           | "M"           | sex of the proband (M=male, F=female) |  

#### Definitions for the required columns in the missing indels candidates file
These sites have been filtered for sites that are present in the child, but not
in the parents, and not examined by denovogear.

| name             | example       | definition                            |
| -----------      | ------------- | -------------                         |
| person_stable_id | "DDDP100001"  | ID of the proband                     |
| chrom            | "1"           | chromosome that the candidate is on   |
| pos              | 1000001       | nucleotide position of the candidate along the chromosome |
| ref              | "A"           | base-code for the reference allele    |
| alt              | "AG"          | base-code for the alternate allele    |
| symbol           | "ATRX"        | HGNC symbol of the gene that the variant lies in (or blank for no gene) |
| dp4_child        | 10,12,9,10    | Comma-separated read-depths in the child for the reference allele in the forward orientation, the reference allele in the reverse orientation, the alternate allele in the forward orientation, and the alternate allele in the reverse orientation. Orientation is with respect to the reference genome. |
| dp4_mother       | 15,17,0,0     | Read-depths in the mother as for dp4_child |
| dp4_father       | 13,14,0,0     | Read-depths in the father as for dp4_child |
| max_af           | 0.0001        | Highest allele frequency for alternate allele from reference healthy control population |
| consequence      | "missense_variant" | VEP-annotated consequence for the variant |
