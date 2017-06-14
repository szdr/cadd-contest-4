import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from rdkit.Chem import Descriptors


def get_desc_names():
    ret = []
    for desc_name, _ in Descriptors.descList:
        if "Charge" in desc_name:
            continue
        ret.append(desc_name)
    return ret


def get_feature_names(nbits=2048):
    ecfp4_names = ["ECFP4_{}".format(i + 1) for i in range(nbits)]
    desc_names = get_desc_names()
    feature_names = ecfp4_names + desc_names
    return feature_names


def smiles2feature(smiles, desc_names):
    """
    featureは辞書で返す
    """
    mol = Chem.MolFromSmiles(smiles)

    # ECFP4
    ecfp4 = Chem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    ecfp4_arr = np.zeros((1, ))
    DataStructs.ConvertToNumpyArray(ecfp4, ecfp4_arr)
    ret = {}
    for i, bit in enumerate(ecfp4_arr):
        ret["ECFP4_{}".format(i + 1)] = int(bit)

    # Descriptor
    desc_name_set = set(get_desc_names())
    for desc_name, desc_func in Descriptors.descList:
        if desc_name not in desc_name_set:
            continue
        ret[desc_name] = desc_func(mol)
    return ret


def smiles2feature_as_array(smiles, desc_names):
    """
    featureはnumpy arrayで返す
    ECFP4, descriptorsの順番
    """
    mol = Chem.MolFromSmiles(smiles)

    # ECFP4
    ecfp4 = Chem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    ecfp4_arr = np.zeros((1, ))
    DataStructs.ConvertToNumpyArray(ecfp4, ecfp4_arr)

    # Descriptor
    # WARNING: Descriptors.descListとdesc_namesの順番は同じ
    descriptor_arr = np.zeros(len(desc_names))
    desc_name_set = set(desc_names)
    ind = 0
    for desc_name, desc_func in Descriptors.descList:
        if desc_name not in desc_name_set:
            continue
        descriptor_arr[ind] = desc_func(mol)
        ind += 1

    return np.r_[ecfp4_arr, descriptor_arr]
