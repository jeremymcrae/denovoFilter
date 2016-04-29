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

from denovoFilter.load_candidates import load_candidates
from denovoFilter.preliminary_filtering import preliminary_filtering, check_coding
from denovoFilter.exclude_segdups import check_segdups
from denovoFilter.missing_symbols import fix_missing_gene_symbols
from denovoFilter.standardise import standardise_columns

def screen_candidates(de_novos_path, fails_path, filter_function, maf=0.01,
        fix_missing_genes=True, annotate_only=False):
    """ load and optionally filter candidate de novo mutations.
    
    Args:
        de_novos_path: path to table of unfiltered canddiate DNMs
        fails_path: path to file listing samples which failed QC, and therefore
            all of their candidates need to be excluded.
        filter_function: function for filtering the candidates, either
            filter_denovogear_sites(), or filter_missing_indels().
        maf: MAF threshold for filtering. This is 0.01 for denovogear sites,
            and 0 for the missing indels.
        fix_missing_genes: whether to annotate HGNC symbols for candidates
            missing these.
        annotate_only: whether to include a column indicating pass status, rather
            than excluding all candidates which fail the filtering.
    
    Returns:
        pandas DataFrame of candidate de novo mutations.
    """
    
    # load the datasets
    de_novos = load_candidates(de_novos_path)
    sample_fails = []
    if fails_path is not None:
        sample_fails = [ x.strip() for x in open(fails_path) ]
    
    # run some initial screening
    status = preliminary_filtering(de_novos, sample_fails, maf_cutoff=maf)
    segdup = check_segdups(de_novos)
    
    if fix_missing_genes:
        de_novos['symbol'] = fix_missing_gene_symbols(de_novos)
    
    pass_status = filter_function(de_novos, status & segdup) & status & segdup
    
    if annotate_only:
        de_novos['pass'] = pass_status
    else:
        de_novos = de_novos[pass_status]
    
    return standardise_columns(de_novos)
