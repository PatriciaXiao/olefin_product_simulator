import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def rxnmix(molecules, nrrxn):
    dimstart = len(molecules)
    moleculestart = np.random.choice(molecules)
    range_doublebond = range(1, moleculestart)
    doublebond = np.random.choice(range_doublebond)
    moleculecleaved = [doublebond, moleculestart - doublebond]
    rest = np.random.choice(moleculecleaved)
    # Start of the reaction
    index = np.argwhere(molecules == moleculestart)[0]
    moleculesreaction = np.delete(molecules, index)
    # Replace here means select nrrxn samples from moleculesreaction, could select the same element twice
    range_molecule_index = range(len(moleculesreaction))
    indexmoleculeselected = np.random.choice(range_molecule_index, nrrxn, replace=True)
    for i in range(nrrxn):
        tmp_index = indexmoleculeselected[i]
        moleculeselected = moleculesreaction[tmp_index]
        range_doublebond = range(1, moleculeselected)
        doublebond = np.random.choice(range_doublebond)
        resttemp = np.random.choice([doublebond, moleculeselected - doublebond])
        moleculesreaction[tmp_index] = rest + moleculeselected - resttemp
        rest = resttemp
    return (moleculesreaction)

def init_molecules12():
    return np.array([12] * 200)

current_init = init_molecules12
n_steps = 200
fname = "C12_{0}.csv".format(n_steps)

repeats = 100
max_slot = 110
min_slot = 0
n_slots = 110

x_sum = np.zeros(n_slots)

for j in range(repeats):
    molecules = current_init()
    mix = rxnmix(molecules, n_steps) # 200 or 1000 or 2000 or 20000
    densedat, bin_edges = np.histogram(mix, bins=n_slots, density=True, range=(min_slot, max_slot))
    x_sum += densedat

y = x_sum / float(repeats)

d = {'y': y}
df = pd.DataFrame(data=d)
df.to_csv(fname, sep='\t')

