"""
各qidに紐づく化合物数のsample_ratio倍数の化合物をZINCからランダムサンプリングする
[output] -> csv
qid, smiles
"""

import csv
import argparse
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input_svmlight")
    parser.add_argument("output")
    parser.add_argument("--sample_ratio", type=int, default=10)  # 決め打ち
    parser.add_argument("--zinc_smiles_file", default="../data/6_p0.smi")
    parser.add_argument("--skip_zinc_read", type=int,
                        default=10)  # 全て読み込むとメモリ重い
    args = parser.parse_args()

    # svmlightを読み込んで各qidに紐づく化合物数を取得する
    n_compound_by_qid = {}
    with open(args.input_svmlight) as fp:
        reader = csv.reader(fp, delimiter=" ")
        for row in reader:
            qid = row[1].split(":")[1]
            n_compound_by_qid[qid] = n_compound_by_qid.get(qid, 0) + 1
    print("Done: read svmlight")

    # zincのsmilesをリストとして読み込む, ランダムサンプリングはindexで行う
    zinc_smiles_list = []
    with open(args.zinc_smiles_file) as fp:
        # zincのファイル形式: SMILES ZINC_ID
        reader = csv.reader(fp, delimiter=" ")
        for i, row in enumerate(reader):
            if i % args.skip_zinc_read == 0:
                zinc_smiles_list.append(row[0])
    N_zinc_smiles_list = len(zinc_smiles_list)
    print("Done: read zinc file")

    # 各qidについてランダムサンプリングし、ファイルに書き込む
    with open(args.output, "w") as out_fp:
        writer = csv.writer(out_fp)
        for qid, n_compound in n_compound_by_qid.items():
            print("qid: {} sampling...".format(qid))
            n_random_sample = n_compound * args.sample_ratio
            random_index = np.random.choice(
                N_zinc_smiles_list, size=n_random_sample, replace=False)
            for ind in random_index:
                writer.writerow([qid, zinc_smiles_list[ind]])
    print("Done: random sampling")
