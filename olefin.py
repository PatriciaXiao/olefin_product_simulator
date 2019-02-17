import numpy as np
import matplotlib.pyplot as plt

def rxnmix(molecules, nrrxn):
    dimstart = len(molecules)
    moleculestart = np.random.choice(molecules)
    range_doublebond = range(1, moleculestart)
    doublebond = np.random.choice(range_doublebond)
    moleculecleaved = [doublebond, moleculestart - doublebond]
    rest = np.random.choice(moleculecleaved)
    # Start of the reaction
    moleculesreaction = molecules[:] # copy of the list
    moleculesreaction.remove(moleculestart)
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
repeats = 1 #100
n_steps = 200
n_slots = 110

x_sum = np.zeros(n_slots)

for j in range(repeats):
    molecules = current_init()
    mix = rxnmix(molecules, n_steps) # 200 or 1000 or 2000 or 20000
    densedat, bin_edges = np.histogram(mix, bins=110, density=True)
    x_sum += densedat

y = x_sum / float(repeats)
print y