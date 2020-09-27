import os
import numpy as np
import pandas as pd 
from datetime import datetime, timedelta
import dynetx as dn


def read_data(filename):
    df = pd.read_csv(filename, sep="\t", header=None)
    df.columns = ["final_time", "Source", "Target"]
    df['final_time'] = pd.to_datetime(df['final_time'],unit='s')
    df['start_time'] = df['final_time'] - timedelta(seconds=30)

    # label encode the target and source IDs
    df["src"] = df["Source"].astype('category').cat.codes
    df["trg"] = df["Target"].astype('category').cat.codes

    # convert time format from UNIX time
    df['start_time_UNIX'] = df['start_time'].astype(np.int64) // 10**9
    df['final_time_UNIX'] = df['final_time'].astype(np.int64) // 10**9
    tmin = df['start_time_UNIX'].min()

    # rescale time intervals
    df['start_time_norm'] = (df['start_time_UNIX'] - tmin) // 20

    # time interval for Gephi
    df['Time Interval'] = "<["+df['start_time'].dt.strftime('%H:%M:%S').values+','+df['final_time'].dt.strftime('%H:%M:%S').values+"]>"

    # export
    # df[['Source','Target','start_time','final_time']].to_excel(filename[:-4]+'_preprocessing.xls',header=True,index=False)
    #nodes = pd.unique(df[['Source', 'Target']].values.ravel('K'))

    # get nodes arrival and leaving time in a format suited for Gephi, i.e.
    # node ID / arrival time / leaving time
    s1 = df.groupby(['Source']).min()['start_time']
    s2 = df.groupby(['Target']).min()['start_time']
    s3 = pd.concat([s1, s2], axis=1).min(axis=1)

    s4 = df.groupby(['Source']).max()['final_time']
    s5 = df.groupby(['Target']).max()['final_time']
    s6 = pd.concat([s4, s5], axis=1).max(axis=1)

    s7 = pd.concat([s3,s6], axis=1)
    s7.index.name = 'Id'

    # DyNetx https://dynetx.readthedocs.io/en/latest/reference/
    g = dn.DynGraph(data=None, edge_removal=True)

    for index, row in df.iterrows(): # a for loop may not be the most efficient way...
        #print(row['start_time_UNIX'] - tmin)
        start_time = row['start_time_norm']
        u = row['src']
        v = row['trg']
        g.add_interaction(u = u, v = v, t = start_time)
    dn.classes.function.freeze(g) # freeze the graph to avoid unwanted changes

    return g


def filehandler():
    for filename in os.listdir('Data/INFECTIOUS_raw'):
        if 'listcontacts_2009_0' in filename:
            path = 'Data/INFECTIOUS_raw/' + filename
            graph = read_data(path)  # graph is a dynetx data structure for the temporal graph
            output_filename = 'Data/Preprocessed/preprocessed_' + filename[18:23]
            dn.write_snapshots(graph, output_filename)
            print('finished preprocessing day ', filename[18:23])
    print('finished preprocessing files')


filehandler()
# to read a file:
# g = dn.read_snapshots('Data/Preprocessed/preprocessed_04_28', nodetype=int, timestamptype=int)
