"""
Copyright (c) 2015 Wellcome Trust Sanger Institute

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from __future__ import absolute_import

import argparse
from datetime import date

import pandas

from denovoFilter.exclude_segdups import exclude_segdups
from denovoFilter.preliminary_filtering import preliminary_filtering
from denovoFilter.missing_symbols import fix_missing_gene_symbols
from denovoFilter.allele_counts import extract_alt_and_ref_counts, \
    count_gene_recurrence
from denovoFilter.site_deviations import test_sites, test_genes
from denovoFilter.set_filter_status import get_filter_status, subset_de_novos
from denovoFilter.remove_redundant_de_novos import get_independent_de_novos
from denovoFilter.missing_indels import filter_missing_indels, load_missing_indels
from denovoFilter.change_last_base_sites import change_conserved_last_base_consequence

DE_NOVOS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/denovo_gear_trios_extracted_passed_variants_11.05.15.tsv"
MISSED_INDELS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novo_data/missed_denovo_indels_datafreeze_2015-04-13.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt"
FAIL_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novo_data/de_novo_sample_fails.txt"
INDEL_FAILS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novo_data/de_novo_sample_fails_missed_indels.txt"
LAST_BASE_PATH = "/lustre/scratch113/projects/ddd/users/jm33/last_base_sites_G.json"
OUTPUT_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.{}.txt".format(str(date.today()))

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description="Filters candidate de novo \
        variants for sites with good characteristics.")
    parser.add_argument("--de-novos", default=DE_NOVOS_PATH, \
        help="Path to file listing candidate de novos.")
    parser.add_argument("--de-novos-indels", default=MISSED_INDELS_PATH, \
        help="Path to file listing candidate de novos indels (not found in \
            the standard de novo filtering).")
    parser.add_argument("--families", default=FAMILIES_PATH, \
        help="Path to file listing family relationships (PED file).")
    parser.add_argument("--sample-fails", default=FAIL_PATH, \
        help="Path to file listing family relationships (PED file).")
    parser.add_argument("--sample-fails-indels", default=INDEL_FAILS_PATH, \
        help="Path to file listing family relationships (PED file).")
    parser.add_argument("--last-base-sites", default=LAST_BASE_PATH, \
        help="Path to file of all conserved last base of exon sites in genome.")
    parser.add_argument("--output", default=OUTPUT_PATH, \
        help="Path to file for filtered de novos.")
    
    args = parser.parse_args()
    
    return args

def include_missing_indels(de_novos, indels_path, fails_path, families_path):
    """ load and filter the missing candidate indels
    """
    
    # load the missing indels datasets and filter for good quality sites
    missing_indels = load_missing_indels(indels_path, families_path)
    sample_fails = [ x.strip() for x in open(fails_path) ]
    missing_indels = filter_missing_indels(missing_indels, sample_fails)
    missing_indels = subset_de_novos(missing_indels)
    
    de_novos = de_novos.append(missing_indels)
    
    return de_novos

def main():
    args = get_options()
    
    # load the datasets
    de_novos = pandas.read_table(args.de_novos, na_filter=False)
    families = pandas.read_table(args.families)
    de_novos = de_novos.merge(families, how="left", left_on="person_stable_id", right_on="individual_id")
    
    sample_fails = [ x.strip() for x in open(args.sample_fails) ]
    
    # run some initial screening, and annotate sites with ref and alt depths
    de_novos = preliminary_filtering(de_novos, sample_fails)
    de_novos = exclude_segdups(de_novos)
    de_novos = fix_missing_gene_symbols(de_novos)
    de_novos = extract_alt_and_ref_counts(de_novos)
    de_novos = count_gene_recurrence(de_novos)
    
    # check if sites deviate from expected strand bias and parental alt depths
    de_novos = test_sites(de_novos)
    de_novos = test_genes(de_novos)
    
    # get the variants that pass the filtering criteria
    pass_status = get_filter_status(de_novos)
    passed = de_novos[(pass_status & de_novos["coding"])]
    
    passed = subset_de_novos(passed)
    if args.de_novos_indels is not None:
        passed = include_missing_indels(passed, args.de_novos_indels, args.sample_fails_indels, args.families)
    
    passed = change_conserved_last_base_consequence(passed, args.last_base_sites)
    
    passed = get_independent_de_novos(passed, args.families)
    passed.to_csv(args.output, sep= "\t", index=False)

if __name__ == '__main__':
    main()
