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

def run_simulator(current_init, n_steps):
    repeats = 100
    max_slot = 110
    min_slot = 0
    n_slots = 110

    x_sum = np.zeros(n_slots)

    for j in range(repeats):
        molecules = current_init()
        mix = rxnmix(molecules, n_steps) 
        densedat, bin_edges = np.histogram(mix, bins=n_slots, density=True, range=(min_slot, max_slot))
        x_sum += densedat

    y = x_sum / float(repeats)
    return y

def save_file(y_dict, fname):
    df = pd.DataFrame(data=y_dict)
    df.to_csv(fname, sep='\t')

def plot_data(fname, total_width=0.8, legend_format=lambda x: x, slice_range=(1, 38), save_fig=None):
    df = pd.read_csv(fname, sep='\t', index_col=0)
    df = df[slice_range[0]: slice_range[1]]
    size = len(df)
    x = np.arange(size)

    y_dict = dict()
    for y in df:
        y_dict[y] = df[y].values
    
    n = len(y_dict.keys())
    width = total_width / float(n)
    x = x - (total_width - width) / 2.0

    i = 0
    for k in y_dict.keys():
        plt.bar(x + width * i, y_dict[k], width=width, label = legend_format(k))
        i += 1

    n_ticks = (slice_range[1] + 2) / 4
    x_positions = [i * 4 + 1 for i in range(n_ticks)]
    x_values = ["C{0}".format(i * 4 + 2) for i in range(n_ticks)]
    plt.xticks(x_positions, x_values)

    plt.legend()

    if save_fig:
        plt.savefig(save_fig)
    plt.show()

def steps_format(key):
    return "{0} steps".format(key)

def init_molecules12():
    return np.array([12] * 200)

def init_molecules_C6_C18(portion = (100, 100)):
    return np.array([6] * portion[0] + [18] * portion[1])

TASK = [1, 2]

# task 1
def __main__():
    if 1 in TASK:
        '''
        print("conducting task 1")
        y = dict()
        n_steps_list = [200, 1000, 2000, 20000]
        for n_steps in n_steps_list:
            print("\t simulating the process of steps {0}".format(n_steps))
            y[n_steps] = run_simulator(init_molecules12, n_steps)
        save_file(y, "task_1.csv")
        '''
        plot_data("task_1.csv", legend_format=steps_format, save_fig="task_1.png")
    if 2 in TASK:
        '''
        print("conducting task 2: 2000 steps, C12 vs C6/C18 1:1")
        n_steps = 2000
        components_dict = {"C12": init_molecules12, "C6/C18 1:1": init_molecules_C6_C18}
        y = dict()
        for component in components_dict.keys():
            print("\t simulating the process of {0}".format(component))
            y[component] = run_simulator(init_molecules12, n_steps)
        save_file(y, "task_2.csv")
        '''
        plot_data("task_2.csv", save_fig="task_2.png")

__main__()


