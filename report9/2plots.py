import matplotlib.pyplot as plt
import numpy as np

def raw_to_spec(path):
    with open(path, "r") as raw:
        return np.array([float(n) for n in raw.read().splitlines()])
    
cs_spec = raw_to_spec("dat/cs137.txt")
cs_bins = np.array([i for i in range(len(cs_spec))])

na_spec = raw_to_spec("dat/na22.txt")
na_bins = np.array([i for i in range(len(na_spec))])
