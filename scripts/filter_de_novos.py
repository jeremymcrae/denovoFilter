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

from denovoFilter.preliminary_filtering import preliminary_filtering, check_coding
from denovoFilter.exclude_segdups import check_segdups
from denovoFilter.missing_symbols import fix_missing_gene_symbols
from denovoFilter.filter_denovogear_sites import filter_denovogear_sites
from denovoFilter.standardise import standardise_columns
from denovoFilter.remove_redundant_de_novos import get_independent_de_novos
from denovoFilter.missing_indels import filter_missing_indels, load_missing_indels
from denovoFilter.change_last_base_sites import change_conserved_last_base_consequence

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description="Filters candidate de novo "
        "variants for sites with good characteristics.")
    parser.add_argument("--de-novos",
        help="Path to file listing candidate de novos.")
    parser.add_argument("--de-novos-indels",
        help="Path to file listing candidate de novos indels (not found in"
            "the standard de novo filtering).")
    parser.add_argument("--families", required=True,
        help="Path to file listing family relationships (PED file).")
    parser.add_argument("--sample-fails",
        help="Path to file listing problematic samples for the denovogear calls.")
    parser.add_argument("--sample-fails-indels",
        help="Path to file listing problematic samples for the indel calls.")
    parser.add_argument("--last-base-sites",
        help="Path to file of all conserved last base of exon sites in genome,"
            " default is to not change anything if this option is not used.")
    
    parser.add_argument("--fix-missing-genes", action='store_true', default=False,
        help="Whether to attempt re-annotation of gene symbols for variants"
            "lacking gene symbols.")
    parser.add_argument("--include-noncoding", action='store_true', default=False,
        help="Whether to include noncoding variants in the output. By default "
            "this will filter to sites within coding regions")
    parser.add_argument("--annotate-only", action='store_true', default=False,
        help="Use if instead of filtering, you want to annotate all candidates"
            "for whether they pass or not.")
    
    parser.add_argument("--output", help="Path to file for filtered de novos.")
    
    args = parser.parse_args()
    
    return args

def check_denovogear_sites(de_novos_path, fails_path, fix_missing_genes=True,
        annotate_only=False):
    ''' load and filter sites identified by denovogear
    '''
    # load the datasets
    de_novos = pandas.read_table(de_novos_path, na_filter=False)
    sample_fails = [ x.strip() for x in open(fails_path) ]
    
    # run some initial screening
    statuses = preliminary_filtering(de_novos, sample_fails)
    segdup_statuses = check_segdups(de_novos)
    
    if fix_missing_genes:
        de_novos['symbol'] = fix_missing_gene_symbols(de_novos)
    
    # get the variants that pass the filtering criteria
    pass_status = filter_denovogear_sites(de_novos, statuses & segdup_statuses)
    
    if not annotate_only:
        de_novos = de_novos[pass_status]
    else:
        de_novos['pass'] = pass_status
    
    de_novos = standardise_columns(de_novos)
    
    return de_novos

def check_missing_indels(indels_path, fails_path, annotate_only=False):
    """ load and filter the missing candidate indels
    """
    
    # load the missing indels datasets and filter for good quality sites
    missing_indels = load_missing_indels(indels_path)
    sample_fails = [ x.strip() for x in open(fails_path) ]
    
    # run some initial screening
    status = preliminary_filtering(missing_indels, sample_fails, maf_cutoff=0)
    segdup_status = check_segdups(missing_indels)
    
    pass_status = filter_missing_indels(missing_indels) & status & segdup_status
    
    if not annotate_only:
        missing_indels = missing_indels[pass_status]
    else:
        missing_indels['pass'] = pass_status
    
    missing_indels = standardise_columns(missing_indels)
    
    return missing_indels

def main():
    args = get_options()
    
    # set a blank dataframe
    de_novos = pandas.DataFrame(columns=["person_stable_id", "chrom", "pos",
        "ref", "alt", "symbol", "consequence", "max_af", "pp_dnm"])
    
    if args.de_novos is not None:
        de_novos = check_denovogear_sites(args.de_novos, args.sample_fails,
            args.fix_missing_genes, args.annotate_only)
    
    if args.de_novos_indels is not None:
        indels = check_missing_indels(args.de_novos_indels,
            args.sample_fails_indels, args.annotate_only)
        de_novos = de_novos.append(indels, ignore_index=True)
    
    if not args.include_noncoding and not args.annotate_only:
        # make sure we are only looking at sites within coding regions
        de_novos = de_novos[check_coding(de_novos)]
    
    if args.last_base_sites is not None:
        de_novos = change_conserved_last_base_consequence(de_novos, args.last_base_sites)
    
    # include the proband sex, since this is necessary to check chrX candidates
    # for being likely pathogenic.
    families = pandas.read_table(args.families, sep='\t')
    sex = dict(zip(families['individual_id'], families['sex']))
    de_novos['sex'] = de_novos['person_stable_id'].map(sex)
    
    de_novos = get_independent_de_novos(de_novos, args.families, args.annotate_only)
    de_novos.to_csv(args.output, sep= "\t", index=False, na_value='NA')

if __name__ == '__main__':
    main()
