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

import os

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

import pandas
import seaborn

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

user_dir = os.path.expanduser("~")

ddd_1k_validations_path = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/DNG_Validation_1133trios_20140130.tsv"
ddd_4k_validations_path = os.path.join(user_dir, "de_novos.ddd_4k.validation_results.2015-09-02.xlsx")
# de_novos_path = os.path.join(user_dir, "de_novos_for_validation.txt")
de_novos_path = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-09-02.txt"
top_genes_path = "/lustre/scratch113/projects/ddd/users/jm33/enriched_genes.ddd_4k.2015-02-02.tsv"
outpath = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-09-22.txt"

def load_de_novo_calls(path):
    """
    """
    
    de_novos = pandas.read_table(path)
    
    # rename the columns
    recode = {"person_stable_id": "person_id", "ref": "ref_allele",
        "alt": "alt_allele", "pos": "start_pos", "symbol": "hgnc"}
    de_novos = de_novos.rename(columns=recode)
    
    de_novos["chrom"] = de_novos["chrom"].astype("str")
    de_novos["end_pos"] = de_novos["start_pos"] + de_novos["ref_allele"].str.len() - 1
    
    de_novos = de_novos[['person_id', 'sex', 'chrom', 'start_pos', 'end_pos',
        'ref_allele', 'alt_allele', 'hgnc', 'var_type', 'consequence', 'max_af',
         'pp_dnm']]
    
    return de_novos

def load_ddd_1k_validations(path):
    """ load the data for the DDD 1K validation efforts
    """
    
    validations = pandas.read_table(path, sep="\t")
    
    # rename the columns I need
    recode = {"person_stable_id": "person_id", "chr": "chrom",
        "ref": "ref_allele", "alt": "alt_allele", "pos": "start_pos",
        "gene_name": "hgnc", "validation_result": "status"}
    validations = validations.rename(columns=recode)
    
    validations["end_pos"] = validations["start_pos"] + validations["ref_allele"].str.len() - 1
    
    # recode the validation status to something easier to interpret
    recode = {"DNM": "de_novo", "FP": "false_positive", "U": "uncertain",
        "INH": "inherited", "P/U": "uncertain", "R": "uncertain"}
    validations["status"] = validations["status"].map(recode)
    
    # restrict the validations table to only a subset of columns. I don't want
    # to know the read depth, genotypes, conserved status etc.
    validations = validations[["person_id", "chrom", "start_pos", "end_pos", \
        "ref_allele", "alt_allele", "status"]]
    
    return validations

def load_ddd_4k_validations(path):
    """ load the data for the DDD 4K validation efforts
    """
    
    validations = pandas.read_excel(path, sheetname="for_Jeremy")
    # make sure the chromosome columns are string types
    validations["chrom"] = validations["CHR"].astype("str")
    validations["person_id"] = validations["ID"]
    validations["start_pos"] = validations["POS"]
    validations["ref_allele"] = validations["REF"]
    validations["alt_allele"] = validations["ALT"]
    validations["end_pos"] = validations["start_pos"] + validations["ref_allele"].str.len() - 1
    
    # recode the validation status
    status = validations["final manual score - edit format to match Best_model"]
    amended_calls = validations["final review JM/Di 02.09.15"]
    status[~amended_calls.isnull()] = amended_calls[~amended_calls.isnull()]
    validations["status"] = status
    
    recode = {"did not validate (pcr/seq fail in child)" : "uncertain",
        "did not validate (pcr/seq fail in child and mum)": "uncertain",
        "dnm": "de_novo", "DNM": "de_novo", "DNM ": "de_novo",
        "DNM mosaic?": "de_novo",
        "DNM mosaic in mother": "de_novo",
        "DNM mosaic in father": "de_novo",
        "DNM (potentially mosaic in father)": "de_novo",
        "fp": "false_positive", "FP": "false_positive",
        "inherited": "inherited", "Inherited": "inherited",
        "inherited/mosaic?": "inherited",
        "Inherited mum (child pcr fail)": "inherited",
        "mosaic?": "inherited",
        "present in proband (P/U)": "uncertain",
        "present in proband (P/U) ": "uncertain",
        "review in meeting - probably de novo": "de_novo", "unclear": "uncertain"}
    validations["status"] = validations["status"].map(recode)
    
    validations = validations[["person_id", "chrom", "start_pos", "end_pos", \
        "ref_allele", "alt_allele", "status"]]
    
    return validations

def find_matching_site(row, de_novos):
    """ this identifies the correct position for variants
    
    Some variants have changed their coordinates during the validation efforts.
    These are all indels, so are due to realignment calling the variant at
    slightly different locations. We can identify the initial position, by
    finding the variant that is closest to the validation coordinates.
    
    Args:
        row: pandas Dataframe of a single row for a candidate, with columns for
            person_id, chrom, start_pos, ref_allele
        de_novos: pandas DataFrame of all candidate sites
    
    Returns:
        correct start position for the variant
    """
    
    # identify the candidate sites for the proband on the same chromosome as the
    # given variant, the figure out how far they are from the given variant. We
    # expect the correct posirion to be within the distance of the ref allele.
    rows = de_novos[(de_novos.person_id == row.person_id) & (de_novos.chrom == row.chrom)]
    rows.delta = abs(rows.start_pos - row.start_pos)
    matches = rows.delta < rows.ref_allele.str.len()
    
    if sum(matches) > 1:
        return None
    elif sum(matches) == 0:
        return None
    
    return int(rows[matches].start_pos)

def fix_incorrect_positions(validations, de_novos):
    """ fix the indels with incorrect positions, by comparing them to the sites
    submitted for validations
    """
    
    # merge the de novo dataset, which includes a HGNC symbol for each candidate
    validations = validations.merge(de_novos, how="left", \
        on=["person_id", "chrom", "start_pos", "end_pos", "ref_allele", "alt_allele"])
    
    # some of the sites have changed positions during the validation efforts
    missing = validations[validations.var_type.isnull()]
    correct_positions = missing.apply(find_matching_site, axis=1, de_novos=de_novos)
    validations.start_pos[validations.var_type.isnull()] = correct_positions
    
    validations = validations[["person_id", "chrom", "start_pos", "end_pos", "ref_allele", \
        "alt_allele", "status"]]
    
    return validations

def count_validated_per_gene(validations):
    """ count the number of sites which validated in each gene
    """
    
    # define whether each site validated or not
    validated = pandas.Series([False] * len(validations["chrom"]))
    validated[validations["status"].isin(["de_novo", "uncertain"])] = True
    validations["validated"] = validated
    
    counts = validations.pivot_table(values=["chrom"], rows=["hgnc"], cols=["validated"], aggfunc=len)
    
    # count the number of validated and invalidated candidates per gene
    hgnc = list(counts.index)
    counts = pandas.DataFrame({"hgnc": hgnc, \
        "validated": counts["chrom"][True].values, \
        "invalidated": counts["chrom"][False].values})
    
    counts["validated"][counts["validated"].isnull()] = 0
    counts["invalidated"][counts["invalidated"].isnull()] = 0
    counts["total"] = counts[["validated", "invalidated"]].sum(axis=1)
    counts["proportion"] = counts["validated"]/counts["total"]
    
    return counts

def plot_validated_per_gene(counts):
    """
    """
    
    counts["3+ candidates"] = counts["total"] > counts["total"].median()
    
    # plot the proportion of candidates validated per gene
    fig = seaborn.FacetGrid(counts, col="3+ candidates", size=6)
    fig.map(pyplot.hist, "proportion")
    fig.set_xlabels("validated per gene")
    fig.set_ylabels("Frequency")
    fig.savefig("proportion_validated.pdf", format="pdf")

def main():
    de_novos = load_de_novo_calls(de_novos_path)
    top_genes = pandas.read_table(top_genes_path)
    
    ddd_1k_results = load_ddd_1k_validations(ddd_1k_validations_path)
    ddd_4k_results = load_ddd_4k_validations(ddd_4k_validations_path)
    ddd_4k_results = fix_incorrect_positions(ddd_4k_results, de_novos)
    validations = ddd_4k_results.append(ddd_1k_results)
    
    # merge the de novo dataset, which includes a HGNC symbol for each candidate
    validations = validations.merge(de_novos, how="left", \
        on=["person_id", "chrom", "start_pos", "end_pos", "ref_allele", "alt_allele"])
    
    validations = validations[["person_id", "chrom", "start_pos", "end_pos", \
        "ref_allele", "alt_allele", "hgnc", "consequence", "status"]]
    validations.to_csv(outpath, sep="\t", index=False)
    
    # counts = count_validated_per_gene(validations)
    # plot_validated_per_gene(counts)

if __name__ == '__main__':
    main()