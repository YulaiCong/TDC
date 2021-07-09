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

# ------------------------CRISPROutcome ------------------
if False:
    from tdc.single_pred import CRISPROutcome

    for name in dataset_names['CRISPROutcome']:
        label_list = retrieve_label_name_list(name)
        data = CRISPROutcome(name=name, path=data_path, label_name=label_list[0])
        split = data.get_split()

        print('-------------------------------- CRISPROutcome ' + name + ' --------------------------------')
        print_stats(split, Drug='GuideSeq')

# ===========================================  Multi-instance Prediction Datasets  ==================================================


# ------------------------DTI ------------------
if False:
    from tdc.multi_pred import DTI

    for name in dataset_names['DTI']:
        data = DTI(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- DTI ' + name + ' --------------------------------')
        print_stats(split, Drug='Drug', Target='Target')

# ------------------------DDI ------------------
if False:
    from tdc.multi_pred import DDI

    for name in dataset_names['DDI']:

        data = DDI(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- DDI ' + name + ' --------------------------------')
        print_stats(split, Drug='Drug1', Target='Drug2')
        print('Number of Classes: ', split['valid']['Y'].max())

        if name == 'drugbank':
            from tdc.utils import get_label_map

            print(get_label_map(name='DrugBank', task='DDI', path=data_path))

# ------------------------PPI ------------------
if False:
    from tdc.multi_pred import PPI

    for name in dataset_names['PPI']:
        data = PPI(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- PPI ' + name + ' --------------------------------')
        print_stats(split, Drug='Protein1', Target='Protein1')
        print('Number of Classes: ', split['valid']['Y'].max())

        # Note: (1) For genes that associate with multiple protein sequences,
        # we separate by * symbol.
        # (2) The dataset contains only positive pairs.
        # All of the unobserved pairs are real negative PPIs, tested experimentally.
        # To get the negative samples, you can call:
        data = data.neg_sample(frac=1)

# ------------------------GDA ------------------
if False:
    from tdc.multi_pred import GDA

    for name in dataset_names['GDA']:
        data = GDA(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- GDA ' + name + ' --------------------------------')
        print_stats(split, Drug='Gene', Target='Disease')
        # print('Number of Classes: ', split['valid']['Y'].max())

# ------------------------DrugRes ------------------
if False:
    from tdc.multi_pred import DrugRes

    for name in dataset_names['DrugRes']:
        data = DrugRes(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- DrugRes ' + name + ' --------------------------------')
        print_stats(split, Drug='Drug', Target='Cell Line')
        # print('Number of Classes: ', split['valid']['Y'].max())

# ------------------------DrugSyn ------------------
if False:
    from tdc.multi_pred import DrugSyn

    for name in dataset_names['DrugSyn']:
        data = DrugSyn(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- DrugSyn ' + name + ' --------------------------------')
        print_stats(split, Drug='Drug1', Target='Drug2')
        # print('Number of Classes: ', split['valid']['Y'].max())

# ------------------------PeptideMHC ------------------
if False:
    from tdc.multi_pred import PeptideMHC

    for name in dataset_names['PeptideMHC']:
        data = PeptideMHC(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- PeptideMHC ' + name + ' --------------------------------')
        print_stats(split, Drug='Peptide', Target='MHC')
        # print('Number of Classes: ', split['valid']['Y'].max())

# ------------------------AntibodyAff ------------------
if False:
    from tdc.multi_pred import AntibodyAff

    for name in dataset_names['AntibodyAff']:
        data = AntibodyAff(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- AntibodyAff ' + name + ' --------------------------------')
        print_stats(split, Drug='Antibody', Target='Antigen')
        # print('Number of Classes: ', split['valid']['Y'].max())

# ------------------------MTI ------------------
if False:
    from tdc.multi_pred import MTI

    for name in dataset_names['MTI']:
        data = MTI(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- MTI ' + name + ' --------------------------------')
        print_stats(split, Drug='miRNA', Target='Target')
        # print('Number of Classes: ', split['valid']['Y'].max())

        # Note: The dataset contains only positive pairs.
        # To get the negative samples, you can call:
        data = data.neg_sample(frac=1)

# ------------------------Catalyst ------------------
if False:
    from tdc.multi_pred import Catalyst

    for name in dataset_names['Catalyst']:
        data = Catalyst(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- Catalyst ' + name + ' --------------------------------')
        print_stats(split, Drug='Reactant', Target='Product')
        print('Number of Classes: ', split['valid']['Y'].max())

        # Note: To know what type of catalyst the label index corresponds to, use:
        from tdc.utils import get_label_map

        get_label_map(name='USPTO_Catalyst', task='Catalyst', path=data_path)

# ===========================================  Generation Datasets  ==================================================


# ------------------------MolGen ------------------
if False:
    from tdc.generation import MolGen

    for name in dataset_names['MolGen']:
        data = MolGen(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- MolGen ' + name + ' --------------------------------')
        print_stats(split, Drug='smiles')

# ------------------------RetroSyn ------------------
if False:
    from tdc.generation import RetroSyn

    for name in dataset_names['RetroSyn']:
        data = RetroSyn(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- RetroSyn ' + name + ' --------------------------------')
        print_stats(split, Drug='input', Target='output')

# ------------------------Reaction ------------------
if True:
    from tdc.generation import Reaction

    for name in dataset_names['Reaction']:
        data = Reaction(name=name, path=data_path)
        split = data.get_split()

        print('-------------------------------- Reaction ' + name + ' --------------------------------')
        print_stats(split, Drug='input', Target='output')

a = 1
