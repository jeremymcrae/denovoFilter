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

from denovoFilter.most_severe import get_most_severe

def tidy_consequence_and_gene_symbols(table):
    ''' clean up the consequence and HGNC symbols if required.
    '''
    
    cq = table['consequence'].copy()
    hgnc = table['symbol'].copy()
    
    # for comma-separated entries,select the first element of the entry, since
    # these are all just duplicates.
    hgnc = hgnc.str.split(',').apply(lambda x: x[0])
    cq = cq.str.split(',').apply(lambda x: x[0])
    
    if not any(cq.str.contains('\|')):
        return table
    
    cq = cq.str.split('|')
    hgnc = hgnc.str.split('|')
    
    # a few variants lack any VEP prediction, just give them an intergenic
    # annotation for now.
    cq = [ x if x != ['-'] else ['intergenic_variant'] for x in cq ]
    
    most_severe = [ get_most_severe(x) for x in cq ]
    
    table['symbol'] = pick_gene_symbols(hgnc, cq, most_severe)
    table['consequence'] = most_severe
    
    return table

def pick_gene_symbols(hgnc, cq, most_severe):
    ''' get the gene symbols for the most severe consequences
    
    Args:
        hgnc: pandas Series of HGNC symbols, each entry has been split by '|'
            separators, so we can have multiple HGNC symbols per candidate.
        cq: list of consequence lists for each candidate.
        most_severe: list of most severe consequence annotations, one for each
            candidate.
    
    Returns:
        list of HGNC symbols, one symbol per candidate.
    '''
    
    # for each variants, identify the positions of the consequences that match
    # the most severe consequence. This is so we can also identify the HGNC
    # symbol matching the most severe consequence.
    positions = []
    for (lis, conseq) in zip(cq, most_severe):
        indices = [ i for i, x in enumerate(lis) if x == conseq ]
        positions.append(indices)
    
    symbols = []
    for pos, genes in zip(positions, hgnc):
        try:
            allowed = [ genes[x] for x in pos ]
        except IndexError:
            allowed = ['.']
        
        # missing symbols are '-' or '.'. Standardise missing symbols to '.'.
        allowed = [ x if x != '-' else '.''' for x in allowed ]
        
        allowed = sorted(set(allowed))
        if len(allowed) > 1 and '.' in allowed:
            allowed.remove('.')
        
        symbols.append(allowed[0])
    
    return symbols
