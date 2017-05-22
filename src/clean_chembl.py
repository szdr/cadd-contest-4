import math
import argparse
import re
import numpy as np
import pandas as pd


def ic50_um_to_pic50(ic50_um):
    return - (np.log10(ic50_um) - 6)


def ic50_to_relevance_score(relation, standard_value, standard_unit,
                            ic50_um_threshold=100):
    """
    relevance_score: large -> good compound
    relation:
    * "=": pIC50 (ic50_um < ic50_um_threshold) else 0
    * ">": 0 (ic50_um >= ic50_um_threshold) else None
    * "<": pIC50 (ic50_um < ic50_um_threshold) else None
    * "NaN": 0 (the compound may not be inhibitor)
    """
    if standard_unit == "nM":
        ic50_um = standard_value * 1e-3
    elif isinstance(standard_unit, float) and math.isnan(standard_unit):
        return 0
    else:
        return None

    if relation == "=":
        if ic50_um < ic50_um_threshold:
            return ic50_um_to_pic50(ic50_um)
        else:
            return 0
    elif relation == ">":
        if ic50_um >= ic50_um_threshold:
            return 0
        else:
            return None
    elif relation == "<":
        if ic50_um < ic50_um_threshold:
            return ic50_um_to_pic50(ic50_um)
        else:
            return None
    else:
        raise ValueError("{} {} {}".format(
            relation, standard_value, standard_unit))


def inhibition_to_relevance_score(relation, standard_value, standard_unit):
    """
    relevance_score: large -> good compound
    relation:
    * "=": max(standard_value, 0)
    * ">": max(standard_value, 0)
    * "<": None
    * "NaN": 0 (the compound may not be inhibitor)
    """
    if standard_unit != "%":
        return None
    elif isinstance(standard_unit, float) and math.isnan(standard_unit):
        return 0

    if relation == "=":
        return max(standard_value, 0)
    elif relation == ">":
        return max(standard_value, 0)
    elif relation == "<":
        return None
    else:
        raise ValueError("{} {} {}".format(
            relation, standard_value, standard_unit))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("chembl_tsv")
    parser.add_argument("outfile")
    args = parser.parse_args()

    df = pd.read_csv(args.chembl_tsv, delimiter="\t")
    df = df[["CMPD_CHEMBLID", "DOC_CHEMBLID", "CANONICAL_SMILES",
             "STANDARD_TYPE", "RELATION", "STANDARD_VALUE", "STANDARD_UNITS"]]

    df_ic50 = df[df.STANDARD_TYPE == "IC50"]
    df_ic50_uniq = df_ic50.drop_duplicates(
        subset=["CMPD_CHEMBLID", "DOC_CHEMBLID"])
    df_ic50_uniq["relevance"] = df_ic50_uniq.apply(
        lambda row: ic50_to_relevance_score(
            row.RELATION, row.STANDARD_VALUE, row.STANDARD_UNITS),
        axis=1)
    df_ic50_uniq.dropna(subset=["relevance"], inplace=True)

    df_inhibition = df[df.STANDARD_TYPE == "Inhibition"]
    df_inhibition_uniq = df_inhibition.drop_duplicates(
        subset=["CMPD_CHEMBLID", "DOC_CHEMBLID"])
    df_inhibition_uniq["relevance"] = df_inhibition_uniq.apply(
        lambda row: inhibition_to_relevance_score(
            row.RELATION, row.STANDARD_VALUE, row.STANDARD_UNITS),
        axis=1)
    df_inhibition_uniq.dropna(subset=["relevance"], inplace=True)

    """
    qid: (SIRT_N, {IC50|INHIBITION}, DOC_ID)
    """
    doc_set = set(list(df_ic50_uniq.DOC_CHEMBLID) +
                  list(df_inhibition_uniq.DOC_CHEMBLID))
    doc_ind_map = {}
    for ind, doc in enumerate(doc_set):
        doc_ind_map[doc] = ind + 1

    sirt_number = re.search("SIRT[0-9]", args.chembl_tsv).group(0)[-1]
    df_ic50_uniq["qid"] = df_ic50_uniq.apply(
        lambda row: sirt_number + "1" + str(doc_ind_map[row.DOC_CHEMBLID]),
        axis=1)
    df_inhibition_uniq["qid"] = df_inhibition_uniq.apply(
        lambda row: sirt_number + "2" + str(doc_ind_map[row.DOC_CHEMBLID]),
        axis=1)

    df_output = pd.concat([df_ic50_uniq, df_inhibition_uniq])

    """
    For some qid, all relevances are 0.0
    They should be removed !
    """
    max_relevances_by_qid = df_output.groupby(["qid"])["relevance"].max()
    drop_qids = set(
        [qid for qid, val in
         max_relevances_by_qid[max_relevances_by_qid < 1e-12].iteritems()])
    df_output = df_output[~df_output.qid.isin(drop_qids)]

    df_output.to_csv(path_or_buf=args.outfile, index=False)


if __name__ == '__main__':
    main()
