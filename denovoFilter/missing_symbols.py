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

import time
import sys
import json

try:
    import urllib.request as request
except ImportError:
    import urllib2 as request

import pandas

PREV_TIME = time.time()
IS_PYTHON3 = sys.version[0] == "3"

def fix_missing_gene_symbols(de_novos):
    """ adds gene symbols to variants lacking them.
    
    Args:
        de_novos: dataframe of de novo variants
    
    Returns:
        data frame of de novos, but with additional annotations for many
        variants previously lacking a HGNC symbol.
    """
    
    # get the variants with no gene annotation, ensure chrom, start and stop
    # positions columns exist
    missing = de_novos[de_novos["symbol"] == ""].copy()
    missing["start"] = missing["pos"]
    missing["end"] = [ x["start"] + len(x["ref"]) - 1 for i, x in missing.iterrows() ]
    
    # find the HGNC symbols (if any) for the variants
    symbols = [ get_gene_id(x["chrom"], x["start"], x["end"], verbose=True) for i, x in missing.iterrows() ]
    de_novos["symbol"][de_novos["symbol"] == ""] = symbols
    
    # 360 out of 17000 de novos still lack HGNC symbols. Their consequences are:
    #
    #   consequence                count
    #   ========================   =====
    #   downstream_gene_variant      17
    #   intergenic_variant          259
    #   regulatory_region_variant    63
    #   upstream_gene_variant        27
    #
    # In spot checks, these are sufficiently distant from genes that we can't
    # add them to the analysis of their nearest gene. We shall analyse these
    # per site by giving them mock gene symbols.
    missing = de_novos[de_novos["symbol"] == ""]
    
    pandas.set_option("mode.chained_assignment", None)
    de_novos["symbol"][de_novos["symbol"] == ""] = \
        [ "fake_symbol.{}_{}".format(x["chrom"], x["pos"]) for i, x in missing.iterrows() ]
    
    return de_novos

def open_url(url, headers):
    """ open url with python libraries
    
    Args:
        url:
        headers:
    
    Returns:
        tuple of http response, http response code, and http headers
    """
    
    req = request.Request(url, headers=headers)
    
    try:
        handler = request.urlopen(req)
    except Exception as e:
        handler = e
    
    status_code = handler.getcode()
    response = handler.read()
    if IS_PYTHON3:
        response = response.decode("utf-8")
    
    # parse the headers into a key, value dictionary
    headers = {}
    for key, value in zip(handler.headers.keys(), handler.headers.values()):
        headers[key.lower()] = value
    
    return response, status_code, headers

def rate_limit_requests():
    """ limit ensembl requests to one per 0.067 s
    """
    
    global PREV_TIME
    rate_limit = 0.0667
    
    diff = time.time() - PREV_TIME
    if diff < rate_limit:
        time.sleep(rate_limit - diff)
    
    PREV_TIME = time.time()

def get_gene_id(chrom, start_pos, end_pos, build="grch37", verbose=False, attempts=0):
    """find the hgnc symbol overlapping a variant position
    
    Args:
        variant: data frame or list for a variant, containing columns named
            "chrom", "start_pos", and "end_pos" for a single variant
        build: genome build to find consequences on
        verbose: flag indicating whether to print variants as they are checked
    
    Returns:
        a character string containing the HGNC symbol.
    """
    
    attempts += 1
    if attempts > 5:
        raise ValueError("too many attempts, figure out why its failing")
    
    rate_limit_requests()
    
    # define parts of the URL
    ext = "overlap/region/human/{0}:{1}-{2}?feature=gene".format(chrom, start_pos, end_pos)
    server_dict = {"grch37": "grch37.", "grch38": ""}
    base_url = "http://{}rest.ensembl.org".format(server_dict[build])
    url = "{0}/{1}".format(base_url, ext)
    headers = {"Content-Type" : "application/json"}
    
    if verbose:
        print("chr{0}:{1}   {2}".format(chrom, start_pos, ext))
    
    response, status_code, requested_headers = open_url(url, headers)
    
    if status_code == 429:
        if "retry-after" in requested_headers:
            time.sleep(float(requested_headers["retry-after"]))
        elif "x-ratelimit-reset" in requested_headers:
            time.sleep(int(requested_headers["x-ratelimit-reset"]))
        return get_gene_id_for_variant(chrom, start_pos, end_pos, build, verbose, attempts)
    elif status_code in [503, 504]:
        time.sleep(30)
        return get_gene_id_for_variant(chrom, start_pos, end_pos, build, verbose, attempts)
    elif status_code != 200:
        raise ValueError("Invalid Ensembl response: {0} for {1}.\nSubmitted URL was: {2}{3}\nheaders: {4}\nresponse: {5}".format(status_code, sequence_id, \
                base_url, ext, requested_headers, response))
    
    json_text = json.loads(response)
    if len(json_text) > 0:
        return json_text[0]["external_name"]
    
    return ""