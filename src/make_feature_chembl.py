import argparse
import pandas as pd
from sklearn.datasets import dump_svmlight_file
from util import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_chembl_csv")
    parser.add_argument("output_svmlight")
    args = parser.parse_args()

    df = pd.read_csv(args.input_chembl_csv)
    desc_names = get_desc_names()
    feature_data_dict = {row[1]["CMPD_CHEMBLID"]: smiles2feature(
        row[1]["CANONICAL_SMILES"], desc_names) for row in df.iterrows()}

    df_feature = pd.DataFrame.from_dict(feature_data_dict, orient="index")
    df_feature["CMPD_CHEMBLID"] = df_feature.index

    df_merged = pd.merge(df, df_feature, on="CMPD_CHEMBLID")
    # for XGBoost
    df_merged_sorted_by_qid = df_merged.sort_values(by="qid")
    feature_names = get_feature_names()
    # feature_namesの順序は一意
    X = df_merged_sorted_by_qid[feature_names].as_matrix()
    ys = df_merged_sorted_by_qid[["relevance"]].as_matrix().flatten()
    qid = df_merged_sorted_by_qid[["qid"]].as_matrix()
    print("X shape: ", X.shape)
    dump_svmlight_file(X, ys, args.output_svmlight, query_id=qid)


if __name__ == '__main__':
    main()
