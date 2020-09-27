from collections import Counter

import dynetx as dn
import matplotlib.pyplot as plt
import matplotlib
import json
matplotlib.rcParams['text.usetex'] = True
import numpy as np
import os
import pickle

def load_graphs(dir = 'Data/Preprocessed/'):
    """
    Construct dictionary containing results of may indexed by the day month (e.g. 05-09)
    """
    all_days = {}
    for filename in os.listdir(dir):
        if 'preprocessed_' in filename:
            path = dir + filename
            all_days[filename[13:]] = dn.read_snapshots(path, nodetype=int, timestamptype=int)
    print("Loaded", len(all_days), "days.")

    return all_days

class Day:
    """
    Class to analyze and simulate SI for any given day.
    It takes as input a preprocessed graph, alpha = transmission probability, prior = probability of a infected visitor.
    """
    def __init__(self, g, alpha = 0.05, prior = 0.1):
        self.g        = g
        self.alpha    = alpha
        self.prior    = prior
        self.status   = {}
        self.time_map = {}
        self.n_inf    = 0

    def set_alpha(self,a):
        self.alpha = a

    def run_si(self, only_inf):
        """
        Run 1 SI simulation. If only_inf=True returns only the total amount of infected at closing time (faster).
        If only_inf=False it stores the current amount of visitors and infected at every timestamp (slower).
        """
        self.status = {}
        count  = 0
        self.time_map = {}
        self.n_inf = 0

        def add_new(u,v,n):
            """
            Add a new node, Infected with probability = prior, otherwise is susceptible.
            """
            if u not in self.status:
                p = np.random.rand()
                self.status[u] = 'I' if p < self.prior else 'S'
                n+=1
            if v not in self.status:
                p = np.random.rand()
                self.status[v] = 'I' if p < self.prior else 'S'
                n+=1
            return n

        def count_infected(status):
            return Counter(status.values())['I']  # faster than looping

        def one_interaction(t, count):
            """
            Simulate one interaction for every S-I edge present at time t.
            :param count: number of nodes up time t
            """
            for u, v, _ in self.g.interactions(t=t):
                count = add_new(u,v,count)
                if self.status[u] != self.status[v]: #not equal corresponds to XOR
                    p = np.random.rand()
                    if p < self.alpha:
                        self.status[u] = 'I'
                        self.status[v] = 'I'

        if only_inf == False:
            for t in range(0,self.g.temporal_snapshots_ids()[-1]+1):
                one_interaction(t,count)
                current_n = self.g.number_of_nodes(t=t)
                infected = count_infected(self.status)
                self.time_map[t] = {'n': count, 'n_current': current_n, 'infected': infected}
        else: # store only final amount of infected
            for t in range(0,self.g.temporal_snapshots_ids()[-1]): # ! the range is different
                one_interaction(t,count)
            one_interaction(self.g.temporal_snapshots_ids()[-1]+1, count)
            self.n_inf = count_infected(self.status)


    def infected_end(self):
        """
        Returns the final amount of infected.
        """
        return self.time_map[self.g.temporal_snapshots_ids()[-1]]['infected'] if len(self.time_map) != 0 else self.n_inf


    def plot_infected(self):
        """
        Returns the time evolution of the total number of infected.
        """
        x = []
        y = []
        for t, val in self.time_map.items():
            x.append(t)
            y.append(val['infected'])

        plt.plot(x, y)
        plt.xlabel(r"Time")
        plt.ylabel(r"\# infected")
        plt.show()

        return x,y

    def simulate(self, alpha, n = 10, final_plot = True): # to be optimized with np.array
        """
        Simulate n SI processes for an array of alpha values.
        :param final_plot: If True, plot the average final amount of infected for every alpha.
        :param n: Number of different runs.
        """
        infected_total = np.empty((alpha.shape[0], n), dtype=np.int16 )
        for j, a in enumerate(alpha):
            self.alpha = a
            x = []
            for i in range(n):
                self.run_si(only_inf = True)
                x.append(self.infected_end())
            infected_total[j, :] = x
            print(j, end=" ")

        infected_total = (infected_total / self.g.number_of_nodes()) * 100 # take percentage
        if final_plot == True:
            plt.errorbar(alpha, infected_total.mean(axis=1), yerr=infected_total.std(axis=1), capsize=2)
            plt.hlines(100 * prior, alpha.min(), alpha.max(), 'r', linestyles="dashed")
            plt.xlabel(r"$\textit{Infection Probability} \: \alpha$")
            plt.ylabel(r"$\langle i \rangle (\alpha) \: [\%]$")
            plt.show()

        return infected_total
def plot_time_evolution(days,alpha):
    """
    Function that plots and saves the time evolution for given days and beta values.
    :param days: Arrays of days labels to loop over. Need object all_days
    :param alpha: Array of beta transmission rates
    :return: Dictionary for each day, with all the simulated SI porcesses and number of infected at every time step
    """
    results = {}
    for day in days:
        print("Day:"+day)
        graph = all_days[day]
        d = Day(graph)
        for j,a in enumerate(alpha):
            d.set_alpha(a)
            d.run_si(only_inf=False)
            if j == 0:
                x = np.array(list(d.time_map.keys()))
                y = np.empty((len(alpha), len(x)), dtype=np.int16)
            for t, val in d.time_map.items():
                y[j,t] = val['infected']

        plt.plot(x, y.T)
        plt.xlabel(r"Time $[t]$")
        plt.ylabel(r"$i(\beta)$")
        plt.legend(np.round(alpha, decimals=4), title=r"$\beta$", fontsize='small')
        plt.savefig("figures/TimeEvolution/TimeEvolution_"+day+".png", dpi=200)
        plt.close()

        results[day] = {'x': x, 'y' : y}

    return results



def final_results():
    for day in days:
        path = 'Data/Preprocessed/preprocessed_' + day
        graph = dn.read_snapshots(path, nodetype=int, timestamptype=int)
        d = Day(graph)
        i_total = np.empty((alpha.shape[0], n))  # alpha x simulations
        i_total = d.simulate(alpha, n, final_plot=False)
        with open('Data/Simulations/final_' + day + '.pkl', 'wb') as output:
            pickle.dump(i_total, output, pickle.HIGHEST_PROTOCOL)
        print('finished', n, 'simulations for day', day)


def pkl_to_json():
    for filename in os.listdir('Data/Simulations'):
        day = filename[6:11]
        if day in days:
            input = "Data/Simulations/final_" + day + ".pkl"
            with open(input, 'rb') as f:
                d = pickle.load(f)
                with open('Data/Simulations/final_' + day + '.txt', 'w') as outfile:
                    json.dump(d.tolist(), outfile)


def read_json(day):
    with open('Data/Simulations/final_' + day + '.txt') as f:
        data = json.load(f)  # alpha x simulations

    return np.array(data)


def read_graph(day):
    path = "Data/Preprocessed/preprocessed_" + day
    graph = dn.read_snapshots(path, nodetype=int, timestamptype=int)

    return graph


def read_instance_day(day):
    path = "Data/Simulations/final_" + day + ".pkl"
    with open(path, 'rb') as f:
        d = pickle.load(f)

    return d


def plot_infected_final(d):
    means = []
    deviations = []
    for simulations in d:
        means.append(np.mean(simulations))
        deviations.append(np.std(simulations))


    plt.errorbar(alpha, means, yerr=deviations, capsize=2)
    #plt.hlines(100 * prior, alpha.min(), alpha.max(), 'r', linestyles="dashed")
    #plt.xlabel(r"$\textit{Infection Probability} \: \alpha$")
    #plt.ylabel(r"$\langle i \rangle_1 (\alpha) \: [\%]$")  # Fix for more days

    plt.show()




np.random.seed(42)
all_days = load_graphs()
alpha = np.linspace(0.0001,0.1,15)
n = 50 # number simulations
prior = 0.1

# # Run simulations over all days
# i_total = np.empty( (m, alpha.shape[0], n) ) # day x alpha x realization
# count = 0
# for i,k,g in enumerate(all_days.items()):
#     print(count, ":", end=" ")
#     d = Day(g)
#     i_total[i] = d.simulate(alpha,n, final_plot=False)
#     count += 1
#     print("")
#
# i_total.mean(axis=2).mean(axis=0)
#
# plt.errorbar(alpha, i_total.mean(axis=2).mean(axis=0), yerr=i_total.std(axis=2).std(axis=0), capsize=2)
# plt.xlabel(r"$\textit{Infection Probability} \: \alpha$")
# plt.ylabel(r"$\langle i \rangle_{total} (\alpha) \: [\%]$")  # Fix for more days
# plt.show()

days = ['05_13','05_20','05_22','06_12','05_10','05_30','06_06']
m = len(days)
alpha = np.linspace(0.0001,0.1,15) # alpha values used for the simulations
n = 100  # number simulations

  # Save this object
with open('Data/Simulations/'+ '_alldays_n_5_' + '.pkl', 'wb') as output:
    pickle.dump(i_total, output, pickle.HIGHEST_PROTOCOL)
#
#  # Load pre-existing data
# with open("Data/Simulations/preprocessed_04_28_n_100_alpha_.pkl", 'rb') as input:
#     x_1 = pickle.load(input)
