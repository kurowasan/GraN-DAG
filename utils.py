import os
import pickle

import torch
import numpy as np


def dump(obj, exp_path, name, txt=False):
    if not txt:
        with open(os.path.join(exp_path, name + ".pkl"), "wb") as f:
            pickle.dump(obj, f)
    else:
        with open(os.path.join(exp_path, name + ".txt"), "w") as f:
            f.write(str(obj))


def load(exp_path, name):
    with open(os.path.join(exp_path, name), "rb") as f:
        obj = pickle.load(f)
    return obj

