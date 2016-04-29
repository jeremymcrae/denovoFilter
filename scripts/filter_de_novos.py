"""
Copyright (c) 2016 Genome Research Ltd.

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
import sys

import pandas

from denovoFilter.screen_candidates import screen_candidates
from denovoFilter.filter_denovogear_sites import filter_denovogear_sites
from denovoFilter.missing_indels import filter_missing_indels
from denovoFilter.change_last_base_sites import change_conserved_last_base_consequence
from denovoFilter.check_independence import check_independence

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
    parser.add_argument("--include-recurrent", action='store_true', default=False,
        help="Use if you want to retain recurrent de novos within a family, or "
            "within an individual (per gene).")
    
    parser.add_argument("--output", default=sys.stdout,
        help="Path to file for filtered de novos. Defaults to standard out.")
    
    args = parser.parse_args()
    
    return args

def main():
    args = get_options()
    
    # set a blank dataframe
    de_novos = pandas.DataFrame(columns=["person_stable_id", "chrom", "pos",
        "ref", "alt", "symbol", "consequence", "max_af", "pp_dnm"])
    
    if args.de_novos is not None:
        de_novos = screen_candidates(args.de_novos, args.sample_fails,
            filter_denovogear_sites, args.fix_missing_genes, args.annotate_only)
    
    if args.de_novos_indels is not None:
        indels = screen_candidates(args.de_novos_indels, args.sample_fails_indels,
            filter_missing_indels, 0.0, args.fix_missing_genes, args.annotate_only)
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
    
    if not args.include_recurrent:
        family_ids = dict(zip(families['individual_id'], families['family_id']))
        independent = check_independence(de_novos, family_ids, args.annotate_only)
        
        if args.annotate_only:
            de_novos['pass'] = de_novos['pass'] & independent
        else:
            de_novos = de_novos[independent]
    
    de_novos.to_csv(args.output, sep= "\t", index=False, na_value='NA')

if __name__ == '__main__':
    main()
