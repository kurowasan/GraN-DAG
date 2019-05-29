import os
import glob
import numpy as np

from collections import namedtuple

DATA = 'data{0}.npy'
DAG = 'DAG{0}.npy'
CPDAG = 'CPDAG{0}.npy'

DAGData = namedtuple('DAGData', 'dag cpdag train valid')

def load_from_id(folder, file_id, valid_split=0.2):
    dag_filename = os.path.join(folder, DAG.format(file_id))
    cpdag_filename = os.path.join(folder, CPDAG.format(file_id))
    data_filename = os.path.join(folder, DATA.format(file_id))

    if os.path.isfile(dag_filename):
        with open(dag_filename, 'rb') as f:
            dag = np.load(f)
    else:
        raise IOError('The DAG file `{0}` was not found.'.format(dag_filename))

    if os.path.isfile(cpdag_filename):
        with open(cpdag_filename, 'rb') as f:
            cpdag = np.load(f)
    else:
        raise IOError('The CPDAG file `{0}` was not found.'.format(cpdag_filename))

    if os.path.isfile(data_filename):
        with open(data_filename, 'rb') as f:
            data = np.load(f)
    else:
        raise IOError('The data file `{0}` was not found.'.format(data_filename))

    num_samples = data.shape[0]
    num_valid = int(valid_split * num_samples)
    data_valid, data_train = data[:num_valid], data[num_valid:]

    return DAGData(dag=dag, cpdag=cpdag, train=data_train, valid=data_valid)


def load_from_folder(folder, valid_split=0.2):
    if not os.path.isdir(folder):
        raise IOError('The folder `{0}` was not found or is not a '
            'folder.'.format(folder))

    dag_filenames = glob.glob(os.path.join(folder, DAG.format('*')))
    num_files = len(dag_filenames)

    try:
        data = [load_from_id(folder, idx, valid_split=valid_split) for idx in range(1, num_files + 1)]
    except IOError:
        raise IOError()

    return data
