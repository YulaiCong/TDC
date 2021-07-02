# import torch
import numpy as np
import pandas as pd

pd.set_option('display.width', None)

from tdc.metadata import dataset_names
from tdc.utils import retrieve_label_name_list


data_path = 'E:/Data/TDC/'


def print_stats(split, Drug='Drug', Target=None):

    print(split['train'].columns)

    print('======= Train Data 1 =======')
    print(split['train'].head(1))

    print('======= Test Data 1 =======')
    print(split['test'].head(1))

    print('Nsample Train-Valid-Test: %d-%d-%d' % (split['train'].shape[0],
                                                  split['valid'].shape[0],
                                                  split['test'].shape[0]))
    if isinstance(split['train'][Drug][0], str):
        print('MaxLen Train-Valid-Test: %d-%d-%d' % (split['train'][Drug].str.len().max(),
                                                     split['valid'][Drug].str.len().max(),
                                                     split['test'][Drug].str.len().max()))
    elif isinstance(split['train'][Drug][0], np.ndarray):
        print('ArrayShape Train-Valid-Test: %d-%d-%d' % (split['train'][Drug][0].shape[0],
                                                         split['valid'][Drug][0].shape[0],
                                                         split['test'][Drug][0].shape[0]))

    if Target is not None:
        print('MaxLen' + Target + ' Train-Valid-Test: %d-%d-%d' % (split['train'][Target].str.len().max(),
                                                     split['valid'][Target].str.len().max(),
                                                     split['test'][Target].str.len().max()))


def print_stats_yields(split):
    print('Nsample Train-Valid-Test: %d-%d-%d' % (split['train'].shape[0],
                                                  split['valid'].shape[0],
                                                  split['test'].shape[0]))

    print(split['train'].columns)

    print('======= Train Data 1 =======')
    print(split['train'].head(1))

    print('======= Test Data 1 =======')
    print(split['test'].head(1))

    print('============================')
    train_data = pd.DataFrame(list(split['train'].Reaction))
    print('Reaction:', train_data.columns)
    train_catalyst_maxlen = train_data['catalyst'].str.len().max()
    train_reactant_maxlen = train_data['reactant'].str.len().max()
    train_product_maxlen = train_data['product'].str.len().max()
    del train_data
    valid_data = pd.DataFrame(list(split['valid'].Reaction))
    valid_catalyst_maxlen = valid_data['catalyst'].str.len().max()
    valid_reactant_maxlen = valid_data['reactant'].str.len().max()
    valid_product_maxlen = valid_data['product'].str.len().max()
    del valid_data
    test_data = pd.DataFrame(list(split['test'].Reaction))
    test_catalyst_maxlen = test_data['catalyst'].str.len().max()
    test_reactant_maxlen = test_data['reactant'].str.len().max()
    test_product_maxlen = test_data['product'].str.len().max()
    del test_data

    print('catalystMaxLen Train-Valid-Test: %d-%d-%d' % (train_catalyst_maxlen,
                                                         valid_catalyst_maxlen,
                                                         test_catalyst_maxlen))
    print('reactantMaxLen Train-Valid-Test: %d-%d-%d' % (train_reactant_maxlen,
                                                         valid_reactant_maxlen,
                                                         test_reactant_maxlen))
    print('productMaxLen Train-Valid-Test: %d-%d-%d' % (train_product_maxlen,
                                                        valid_product_maxlen,
                                                        test_product_maxlen))

# ===========================================  Single-instance Prediction Datasets  ==================================================

# -------------------ADME-----------------------
if False:
    from tdc.single_pred import ADME

    for name in dataset_names['ADME']:
        data = ADME(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- ADME ' + name + ' --------------------------------')
        print_stats(split)

# ---------------------Tox---------------------
if False:

    from tdc.single_pred import Tox

    for name in dataset_names['Toxicity']:

        if name == 'tox21':
            label_list = retrieve_label_name_list('Tox21')
            data = Tox(name=name, path=data_path, label_name=label_list[0])
        else:
            data = Tox(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- Tox ' + name + ' --------------------------------')
        print_stats(split)

# ---------------------HTS---------------------
if False:
    from tdc.single_pred import HTS

    for name in dataset_names['HTS']:
        data = HTS(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- HTS ' + name + ' --------------------------------')
        print_stats(split)

# --------------------QM----------------------
if False:
    from tdc.single_pred import QM

    for name in dataset_names['QM']:
        label_list = retrieve_label_name_list(name)
        data = QM(name=name, path=data_path, label_name=label_list[0])
        split = data.get_split()

        print('-------------------------------- QM ' + name + ' --------------------------------')
        print_stats(split)

# ------------------------Yields------------------
if False:
    from tdc.single_pred import Yields

    for name in dataset_names['Yields']:
        data = Yields(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- Yields ' + name + ' --------------------------------')
        print_stats_yields(split)


# ------------------------Paratope ------------------
if False:
    from tdc.single_pred import Paratope

    for name in dataset_names['Paratope']:
        data = Paratope(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- Paratope ' + name + ' --------------------------------')
        print_stats(split, Drug='Antibody')


# ------------------------Epitope ------------------
if False:
    from tdc.single_pred import Epitope

    for name in dataset_names['Epitope']:
        data = Epitope(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- Epitope ' + name + ' --------------------------------')
        print_stats(split, Drug='Antigen')


# ------------------------Develop ------------------
if False:
    from tdc.single_pred import Develop

    for name in dataset_names['Develop']:

        label_list = retrieve_label_name_list(name) if name == 'tap' else [None, ]
        data = Develop(name=name, path=data_path, label_name=label_list[0])
        split = data.get_split()

        print('-------------------------------- Develop ' + name + ' --------------------------------')
        print_stats(split, Drug='Antibody')


# ------------------------Develop ------------------
if False:
    from tdc.single_pred import CRISPROutcome

    for name in dataset_names['CRISPROutcome']:

        label_list = retrieve_label_name_list(name)
        data = CRISPROutcome(name=name, path=data_path, label_name=label_list[0])
        split = data.get_split()

        print('-------------------------------- CRISPROutcome ' + name + ' --------------------------------')
        print_stats(split, Drug='GuideSeq')



# ===========================================  Multi-instance Prediction Datasets  ==================================================


if True:
    from tdc.multi_pred import DTI

    for name in dataset_names['DTI']:

        data = DTI(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- DTI ' + name + ' --------------------------------')
        print_stats(split, Drug='Drug', Target='Target')




a = 1
