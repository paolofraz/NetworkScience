import pandas as pd
import dynetx as dn
import matplotlib.pyplot as plt
import numpy as np
import os
import json


def burstiness(g, window_size):

    interactions = []
    for t in range(0, g.temporal_snapshots_ids()[-1] + 1):
        val = g.size(t) if g.number_of_nodes(t=t) > 0 else 0
        interactions.append(val)

    interactions = pd.Series(interactions)

    if window_size > 1:
        rolling = interactions.rolling(window=window_size)
        rolling = rolling.mean()

        val = 0
        for i in range(0, window_size-1):
            val = val + interactions[i]
            rolling[i] = val / window_size

        interactions = rolling

    interactions = interactions.tolist()
    return interactions


def burstiness_window(g,windows):
    data = {}
    for window in windows:
        interactions = burstiness(g, window)
        fano = compute_fano(interactions)
        data[window] = (fano, interactions)

    return data


def compute_fano(x): # Idea: compute fano factor over different windows, base window size on average time in museum
    fano = np.var(x)/np.mean(x)
    return fano


def preprocess(windows):
    data = {}
    for filename in os.listdir('Data/Preprocessed'):
        if 'preprocessed_0' in filename:
            g = dn.read_snapshots(str('Data/Preprocessed/' + filename), nodetype=int, timestamptype=int)
            day = filename[13:18]
            print(day)
            data[day] = burstiness_window(g, windows)

    with open('Data/Preprocessed_bursts.txt', 'w') as outfile:
        json.dump(data, outfile)


def read_data():
    with open('Data/Preprocessed_bursts.txt') as json_file:
        data = json.load(json_file)

    print('retrieved data')
    return data

def read_data_long():
    with open('Data/Preprocessed_bursts_large.txt') as json_file:
        data = json.load(json_file)

    print('retrieved data')
    return data

def plot_window(data, day, min):
    fig, ax = plt.subplots()

    for window, result in data.items():
        window = int(window)
        fano = result[0]
        text = '20 seconds' if window == 1 else 'W = ' +  str(window)
        text = text + ', Burstiness: ' + str(fano)[:4]
        plt.plot(result[1], label=text, alpha=0.8, linewidth=1)

    ax.legend(loc='upper left', frameon=False, title='Window Size')
    plt.xlabel("Time (t)")
    plt.ylabel("Mean of number of interactions in window")
    filemap = 'min' if min else 'max'
    location = 'Data/Results_bursts/' + filemap + '/figure_' + day
    plt.savefig(location)
    plt.show()


def compute_max_min(data, size):
    sorted_data = sorted((value['10'][0], key) for (key, value) in data.items())
    for i in range(0, size):
        day = sorted_data[i][1]
        plot_window(data[day], day, True)

    for i in range(len(sorted_data)-size, len(sorted_data)):
        day = sorted_data[i][1]
        plot_window(data[day], day, False)



data = read_data()
compute_max_min(data, 15) #compute figures for 15 largest and smallest burstiness days

