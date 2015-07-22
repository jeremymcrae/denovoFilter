""" script to identify candidate de novos that pass certain criteria
"""

from __future__ import absolute_import

import argparse

import pandas

from denovoFilter.exclude_segdups import exclude_segdups
from denovoFilter.preliminary_filtering import preliminary_filtering
from denovoFilter.missing_symbols import fix_missing_gene_symbols
from denovoFilter.allele_counts import extract_alt_and_ref_counts, \
    count_gene_recurrence
from denovoFilter.site_deviations import test_sites, test_genes
from denovoFilter.set_filter_status import get_filter_status, subset_de_novos
from denovoFilter.remove_redundant_de_novos import get_independent_de_novos

DE_NOVOS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/denovo_gear_trios_extracted_passed_variants_11.05.15.tsv"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt"
OUTPUT_PATH = "de_novos.ddd_4k.ddd_only.txt"

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description="Filters candidate de novo \
        variants for sites with good characteristics.")
    parser.add_argument("--de-novos", default=DE_NOVOS_PATH, \
        help="Path to file listing candidate de novos.")
    parser.add_argument("--families", default=FAMILIES_PATH, \
        help="Path to file listing family relationships (PED file).")
    parser.add_argument("--output", default=OUTPUT_PATH, \
        help="Path to file for filtered de novos.")
    
    args = parser.parse_args()
    
    return args

def main():
    args = get_options()
    
    de_novos = pandas.read_table(args.de_novos, na_filter=False)
    de_novos = preliminary_filtering(de_novos)
    de_novos = exclude_segdups(de_novos)
    de_novos = fix_missing_gene_symbols(de_novos)
    
    de_novos = extract_alt_and_ref_counts(de_novos)
    de_novos = count_gene_recurrence(de_novos)
    
    de_novos = test_sites(de_novos)
    de_novos = test_genes(de_novos)
    
    pass_status = get_filter_status(de_novos)
    passed = de_novos[(pass_status & de_novos["coding"])]
    
    passed = subset_de_novos(passed)
    passed = get_independent_de_novos(passed, args.families)
    
    passed.to_csv(args.output, sep= "\t", index=False)

if __name__ == '__main__':
    main()
