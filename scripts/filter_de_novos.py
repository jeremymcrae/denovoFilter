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
FAIL_PATH = "/nfs/users/nfs_j/jm33/apps/denovoFilter/data/sample_fails.txt"
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
    parser.add_argument("--sample-fails", default=FAIL_PATH, \
        help="Path to file listing family relationships (PED file).")
    parser.add_argument("--output", default=OUTPUT_PATH, \
        help="Path to file for filtered de novos.")
    
    args = parser.parse_args()
    
    return args

def main():
    args = get_options()
    
    # load the datasets
    de_novos = pandas.read_table(args.de_novos, na_filter=False)
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
    passed = get_independent_de_novos(passed, args.families)
    
    passed.to_csv(args.output, sep= "\t", index=False)

if __name__ == '__main__':
    main()
