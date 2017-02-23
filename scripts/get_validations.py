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

import argparse

from datetime import date

import pandas

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description="Generate a table of all the"
        "validations results.")
    parser.add_argument("--ddd-1k-validations", \
        default="/nfs/ddd0/Data/datafreeze/1133trios_20131218/DNG_Validation_1133trios_20140130.tsv", \
        help="Path to file listing candidate de novos.")
    parser.add_argument("--updated-validations", \
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novo_data/de_novos.validation_results.2017-02-22.xlsx", \
        help="Path to file listing candidate de novos indels (not found in \
            the standard de novo filtering).")
    parser.add_argument("--low-pp-dnm", \
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novo_data/de_novos.validation_results.low_pp_dnm.2017-02-22.xlsx", \
        help="Path to file listing family relationships (PED file).")
    parser.add_argument("--de-novos", \
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_8k.2016-11-28.txt", \
        help="Path to file listing family relationships (PED file).")
    parser.add_argument("--output", \
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.{}.txt".format(str(date.today())), \
        help="Path to file for filtered de novos.")
    
    return parser.parse_args()

def load_de_novo_calls(path):
    """ load a dataset of filtered de novo calls
    
    Args:
        path: path to de novo candidates dataset
    
    Returns:
        pandas dataframe of candidate de novos, restricted to to most pertinent
        columns.
    """
    
    de_novos = pandas.read_table(path)
    
    # rename the columns
    recode = {"person_stable_id": "person_id", "ref": "ref_allele",
        "alt": "alt_allele", "pos": "start_pos", "symbol": "hgnc"}
    de_novos = de_novos.rename(columns=recode)
    
    de_novos["chrom"] = de_novos["chrom"].astype("str")
    de_novos["end_pos"] = de_novos["start_pos"] + de_novos["ref_allele"].str.len() - 1
    
    return de_novos[['person_id', 'sex', 'chrom', 'start_pos', 'end_pos',
        'ref_allele', 'alt_allele', 'hgnc', 'consequence', 'max_af', 'pp_dnm']]

def load_ddd_1k_validations(path):
    """ load the data for the DDD 1K validation efforts
    
    Args:
        path: path to dataset for DDD 1K validation results
    
    Returns:
        pandas dataframe of candidates, restricted to specific columns
    """
    
    data = pandas.read_table(path, sep="\t")
    
    # rename the columns I need
    recode = {"person_stable_id": "person_id", "chr": "chrom",
        "ref": "ref_allele", "alt": "alt_allele", "pos": "start_pos",
        "gene_name": "hgnc", "validation_result": "status"}
    data = data.rename(columns=recode)
    
    data["end_pos"] = data["start_pos"] + data["ref_allele"].str.len() - 1
    
    # recode the validation status to something easier to interpret
    recode = {"DNM": "de_novo", "FP": "false_positive", "U": "uncertain",
        "INH": "inherited", "P/U": "uncertain", "R": "uncertain"}
    data["status"] = data["status"].map(recode)
    
    # restrict the validations table to only a subset of columns. I don't want
    # to know the read depth, genotypes, conserved status etc.
    return data[["person_id", "chrom", "start_pos", "end_pos",
        "ref_allele", "alt_allele", "status"]]

def load_updated_validations(path):
    """ load the data for the DDD 4K validation efforts
    
    Args:
        path: path to dataset for DDD 4K validation results
    
    Returns:
        pandas dataframe of candidates, restricted to specific columns
    """
    
    data = pandas.read_excel(path, sheetname="Sheet1")
    # make sure the chromosome columns are string types
    data["chrom"] = data["chrom"].astype("str")
    data["end_pos"] = data["start_pos"] + data["ref_allele"].str.len() - 1
    
    return data[["person_id", "chrom", "start_pos", "end_pos", \
        "ref_allele", "alt_allele", "status"]]

def load_low_pp_dnm_validations(path, de_novos):
    """ load the data for the DDD 4K low pp_dnm validation efforts
    
    Args:
        path: path to dataset for DDD 4K validation results
    
    Returns:
        pandas dataframe of candidates, restricted to specific columns
    """
    
    data = pandas.read_excel(path, sheetname="Sheet1")
    # make sure the chromosome columns are string types
    data["chrom"] = data["chrom"].astype("str")
    data["start_pos"] = data["pos"]
    
    data = data.merge(de_novos, how="inner",
        on=["person_id", "chrom", "start_pos"])
    
    data["end_pos"] = data["start_pos"] + data["ref_allele"].str.len() - 1
    
    return data[["person_id", "chrom", "start_pos", "end_pos", \
        "ref_allele", "alt_allele", "status"]]

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
        return -9999999
    elif sum(matches) == 0:
        return -9999999
    
    return int(rows[matches].start_pos)

def fix_incorrect_positions(data, de_novos):
    """ fix the indels with incorrect positions, by comparing them to the sites
    submitted for validations
    
    Args:
        data: pandas dataframe of validation results, where some have
            been swapped to incorrect coordinates.
        de_novos: pandas dataframe of all de novo candidate, so we can match the
            original call coordinates.
    
    Returns:
        pandas DataFrame of validations. where the sites with incorrect
        positions have been fixed.
    """
    
    # merge the de novo dataset, which includes a HGNC symbol for each candidate
    data = data.merge(de_novos, how="left",
        on=["person_id", "chrom", "start_pos", "end_pos", "ref_allele", "alt_allele"])
    
    # some of the sites have changed positions during the validation efforts
    missing = data[data.consequence.isnull()]
    correct_positions = missing.apply(find_matching_site, axis=1, de_novos=de_novos)
    data.start_pos[data.consequence.isnull()][correct_positions > 0] = correct_positions[correct_positions > 0]
    
    return data[["person_id", "chrom", "start_pos", "end_pos", "ref_allele", \
        "alt_allele", "status"]]

def count_validated_per_gene(data):
    """ count the number of sites which validated in each gene
    """
    
    # define whether each site validated or not
    validated = pandas.Series([False] * len(data["chrom"]))
    validated[data["status"].isin(["de_novo", "uncertain"])] = True
    data["validated"] = validated
    
    counts = data.pivot_table(values=["chrom"], rows=["hgnc"], cols=["validated"], aggfunc=len)
    
    # count the number of validated and invalidated candidates per gene
    hgnc = list(counts.index)
    counts = pandas.DataFrame({"hgnc": hgnc,
        "validated": counts["chrom"][True].values,
        "invalidated": counts["chrom"][False].values})
    
    counts["validated"][counts["validated"].isnull()] = 0
    counts["invalidated"][counts["invalidated"].isnull()] = 0
    counts["total"] = counts[["validated", "invalidated"]].sum(axis=1)
    counts["proportion"] = counts["validated"]/counts["total"]
    
    return counts

def remove_duplicates(data):
    """ removes duplicate validations
    
    Args:
        data: dataframe of validation data
    
    Returns:
        dataframe of validation data with duplicate rows removed, and where
        duplicates exist, we check if at least one of the pair has been validated.
    """
    
    columns = ["person_id", "chrom", "start_pos"]
    first = data.duplicated(take_last=False, cols=columns)
    second = data.duplicated(take_last=True, cols=columns)
    dups = data[first | second]
    without_dups = data[~(first | second)]
    
    # some of the duplicates have different validation status codes, such as one
    # being annotated a "uncertain, while the other is annotated as "de_novo".
    # We want to capture if at least on of the pair is "de_novo".
    fixed = pandas.DataFrame(columns=dups.columns)
    for key, x in dups.groupby(["person_id", "chrom", "start_pos"]):
        row = dict(x.iloc[0])
        if "de_novo" in list(x["status"]):
            row["status"] = "de_novo"
        fixed = fixed.append(row, ignore_index=True)
    
    return without_dups.append(fixed)

def main():
    args = get_options()
    de_novos = load_de_novo_calls(args.de_novos)
    
    ddd_1k_results = load_ddd_1k_validations(args.ddd_1k_validations)
    updated_results = load_updated_validations(args.updated_validations)
    updated_results = fix_incorrect_positions(updated_results, de_novos)
    validations = updated_results.append(ddd_1k_results)
    
    low_pp_dnm = load_low_pp_dnm_validations(args.low_pp_dnm, de_novos)
    validations = validations.append(low_pp_dnm)
    validations = remove_duplicates(validations)
    
    # merge the de novo dataset, which includes a HGNC symbol for each candidate
    validations = validations.merge(de_novos, how="left", \
        on=["person_id", "chrom", "start_pos", "end_pos", "ref_allele", "alt_allele"])
    
    validations = validations[["person_id", "chrom", "start_pos", "end_pos", \
        "ref_allele", "alt_allele", "hgnc", "consequence", "status"]]
    validations.to_csv(args.output, sep="\t", index=False)

if __name__ == '__main__':
    main()
