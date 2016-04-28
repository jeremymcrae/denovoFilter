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

import pandas
import numpy

from denovoFilter.most_severe import get_most_severe

def get_person_recurrences(de_novos):
    """ remove de novos recurrent in a gene within individuals.
    
    Find the de novos that are recurrent within a single individual in a
    single gene. We shall treat these as a single de novo event. Prioritise
    including the most severe event within a gene, then take the first variant
    left after that.
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        pandas Series for whether each canddiate is a duplicate of not
    """
    
    # find the variants which are recurrent within a person in a single gene
    from_start = de_novos.duplicated(["person_stable_id", "symbol"])
    from_end = de_novos.duplicated(["person_stable_id", "symbol"], take_last=True)
    
    person_dups = from_start | from_end
    in_person_dups = de_novos[person_dups]
    
    # split the dataset, so we can process gene by gene
    genes = in_person_dups.groupby(["person_stable_id", "symbol"])
    
    # pick a variant for each person, the first of the most severe consequence
    retain = pandas.Series([], dtype=numpy.bool_)
    for key, gene in genes:
        consequence = get_most_severe(gene["consequence"])
        first = gene[gene["consequence"] == consequence].index[0]
        
        gene_retain = pandas.Series([True] * len(gene), index=gene.index)
        gene_retain[first] = False
        retain = retain.append(gene_retain)
    
    # set the selected de novos
    person_dups.loc[retain.index] = retain
    
    return person_dups

def get_independent_de_novos(de_novos, families_path, annotate_only=False):
    """ exclude de novos that would otherwise lead to double counting.
    
    Remove de novos that are shared between multiple probands of a family or
    that are recurrent in a single gene in a single person.
    
    Args:
        de_novos: dataframe of de novo variants
        families_path: path to pedigree file for cohort.
        annotate_only: whether to annotate the canddiate
    
    Returns:
        dataframe with duplicated sites removed
    """
    
    families = pandas.read_table(families_path, na_filter=False)
    families = dict(zip(families['individual_id'], families['family_id']))
    de_novos['family_id'] = de_novos['person_stable_id'].map(families)
    
    # restrict ourselves to the non-duplicates (this retains the first de novo
    # for each family)
    family_dups = de_novos[["family_id", "chrom", "pos", "ref", "alt"]].duplicated()
    person_dups = get_person_recurrences(de_novos)
    
    # find the candidates which are not duplicates
    non_dup = ~family_dups & ~person_dups
    
    if not annotate_only:
        de_novos = de_novos[non_dup]
    else:
        de_novos['pass'] = de_novos['pass'] & non_dup
    
    de_novos = de_novos.drop('family_id', axis=1)
    
    return de_novos
