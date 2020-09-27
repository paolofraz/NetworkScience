import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from visualization import vis
import datetime
import calendar
import os


class day:
    def __init__(self, graph, month_day):
        self.graph = graph
        self.month_day = month_day

    @property
    def weekday(self):
        day_nozeroes = self.month_day[4] if self.month_day[3] == '0' else self.month_day[3:5]
        date = datetime.datetime(2009, int(self.month_day[1]), int(day_nozeroes))
        return [date.weekday(), calendar.day_name[date.weekday()]]

    @property
    def largest_connected_subgraph(self):
        return max(nx.connected_component_subgraphs(self.graph), key=len)

    @property
    def degrees(self):
        return [self.graph.degree(n) for n in self.graph.nodes()]

    @property
    def weighted_degrees(self):
        return [self.graph.degree(n, weight='weight') for n in self.graph.nodes()]

    @property
    def avg_degree(self):
        return sum(self.degrees)/len(self.degrees)

    @property
    def avg_shortest_path(self):
        return nx.average_shortest_path_length(self.largest_connected_subgraph)

    @property
    def diameter(self):
        return nx.diameter(self.largest_connected_subgraph)

    @property
    def avg_clustering(self):
        return nx.average_clustering(self.graph)

    @property
    def clustering(self):
        return nx.clustering(self.graph)

    def visualize(self):
        vis(self.graph)


def plot_degree(G):
    mx = max(G.degrees)
    bins = np.arange(mx+2)-0.5
    plt.xticks(range(mx+1))
    plt.hist(G.degrees, bins, align='mid', density=True, rwidth=0.5, edgecolor='black')
    plt.show()


# compute boxplots for degree distribution for multiple days
def aggregate_hist(all_days):
    agg_hist = []
    for i in range(max_degree+1):
        agg_hist.append([])

    fig, ax = plt.subplots()
    bins = np.arange(max_degree + 2) - 0.5
    for d in all_days:
        hist, bins = np.histogram(all_days[d].degrees, bins=bins)
        total = sum(hist)
        for i in range(len(hist)):
            agg_hist[i].append(hist[i]/total)

    plt.boxplot(agg_hist)
    labels = [item.get_text() for item in ax.get_xticklabels()]
    for index, item in enumerate(ax.get_xticklabels()):
        labels[index] = index  # for some reason the x-tick labels are 1 too large, this decreases them by one
        item.set_visible((index % 5) == 0)
    ax.set_xticklabels(labels)
    plt.ylabel('P(degree)')
    plt.xlabel('Degree')
    plt.show()


def aggregate_graphs(all_days):
    degrees = []
    nodes = []
    edges = []

    for d in all_days:
        G = all_days[d]
        nodes.append(G.graph.number_of_nodes())
        edges.append(G.graph.number_of_edges())
        degrees.append(G.avg_degree)

    plt.scatter(nodes, edges, c=degrees, cmap='YlOrRd', alpha=1)
    plt.xlabel('Number of Nodes')
    plt.ylabel('Number of Edges')
    plt.grid(True)
    plt.colorbar().set_label("Average Degree")
    plt.show()


def aggregate_distances(all_days):
    plt.grid()
    for d in all_days:
        G = all_days[d]
        shortest = G.avg_shortest_path
        diameter = G.diameter
        plt.scatter(diameter, shortest, c='purple', alpha=0.4)

    plt.xlabel('Network Diameter')
    plt.ylabel('Average Shortest Path')
    plt.show()


def aggregate_clustering(all_days):
    clustering = []
    for d in all_days:
        G = all_days[d]
        clustering.append(G.avg_clustering)

    bins = np.arange(0.30, 0.675, 0.025)
    average = sum(clustering)/len(clustering)
    plt.hist(clustering, bins=bins, rwidth=0.95)
    plt.xlabel(
        'Network Clustering Coëfficient for ' + str(len(clustering)) + ' days. (Average = ' + str(average)[0:4] + ')'
    )
    plt.xticks(bins[::2])
    plt.grid()
    plt.show()


def correlate_weighted_degrees(all_days):
    degrees = np.array([])
    weighted_degrees = np.array([])
    for day in all_days:
        G = all_days[day]
        degrees = np.append(degrees, G.degrees)
        weighted_degrees = np.append(weighted_degrees, G.weighted_degrees)

    plt.plot(degrees, weighted_degrees, alpha=0.1, c='Purple', linestyle='', marker='o', markersize=6.0)
    plt.xlabel('Degree')
    plt.ylabel('Weighted Degree')
    plt.show()
    cor_coef = np.corrcoef(degrees, weighted_degrees)
    print("correlation is: ", cor_coef)


def read_data():
    # construct dictionary containing results of may indexed by the day month (e.g. 05-09)
    data = {}
    max_degree = 0  # maintain max degree in data for more efficient later use
    max_avg_degree = 0  # maintain max_avg degree for more efficient later use
    for filename in os.listdir('Data/INFECTIOUS_daily'):
        if 'cumulative_2009-0' in filename:
            path = 'Data/INFECTIOUS_daily/' + filename
            data[filename[16:21]] = day(nx.read_gml(str(path), label='id'), filename[16:21])
            max_degree = max(max_degree, max(data[filename[16:21]].degrees))
            max_avg_degree = max(max_avg_degree, data[filename[16:21]].avg_degree)

    print('finished collecting data')
    return data, max_avg_degree, max_degree


def compute_results(data):
    aggregate_hist(data)
    aggregate_graphs(data)
    aggregate_distances(data)
    aggregate_clustering(data)
    correlate_weighted_degrees(data)
    print('Finished')


#collect data and call function to compute all results
all_days, max_avg_degree, max_degree = read_data()
compute_results(all_days)




