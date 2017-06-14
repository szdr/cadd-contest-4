"""
ZINCからランダムサンプリングした化合物の特徴量生成に使う
ランダムサンプリング結果のファイルはqid,SMILESのcsv
"""

import argparse
import csv
import numpy as np
from sklearn import datasets
from util import get_desc_names, smiles2feature_as_array


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_smi_csv")
    parser.add_argument("output_svmlight")
    parser.add_argument("--dummy_relevance", type=int, default=0)
    args = parser.parse_args()

    X = []
    ys = []
    qid_lst = []
    desc_names = get_desc_names()
    with open(args.input_smi_csv) as input_fp:
        reader = csv.reader(input_fp)
        for i, row in enumerate(reader):
            if i % 1000 == 0:
                print("{} reading...".format(i))
            qid = int(row[0])
            smi = row[1]
            try:
                feature = smiles2feature_as_array(smi, desc_names)
                X.append(feature)
                ys.append(args.dummy_relevance)
                qid_lst.append(qid)
            except:
                print("{} mol is error".format(i))
                pass
    X = np.array(X)
    ys = np.array(ys)
    qid_lst = np.array(qid_lst)
    print("X shape: ", X.shape)
    datasets.dump_svmlight_file(X, ys, args.output_svmlight, query_id=qid_lst)


if __name__ == '__main__':
    main()
