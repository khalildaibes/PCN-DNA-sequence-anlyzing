
from prerequest import *
from matplotlib import pyplot as plt

def get_spike_for_kmer(key):
    appearances_count = spike_dic.get(key)
    if float(appearances_count[0] + appearances_count[1]) * 0.9 <= float(appearances_count[0]):
        return int(1)
    elif float(appearances_count[0] + appearances_count[1]) * 0.9 <= float(appearances_count[1]):
        return int(0)


def get_spike_plus_sum(key):
    appearances_count = spike_dic.get(key)
    return appearances_count[0]


def get_spike_minus_sum(key):
    appearances_count = spike_dic.get(key)
    return appearances_count[1]


def plot_spikes(figure_directory):
    data = pd.DataFrame()
    spike_plus_sum = 0
    spike_minus_sum = 0
    data["X"] = list(spike_dic.keys())
    tmplst = []
    for x in spike_dic.keys():
        tmplst.append(get_spike_for_kmer(x))
        spike_plus_sum += get_spike_plus_sum(x)
        spike_minus_sum += get_spike_minus_sum(x)
    data["Y"] = tmplst
    font1 = {'family': 'serif', 'color': 'blue', 'size': 20}
    font2 = {'family': 'serif', 'color': 'darkred', 'size': 15}
    plt.title("number of spike plus =" + str(spike_plus_sum) + " number of spike minus =" + str(spike_minus_sum), fontdict=font1)
    plt.scatter(x="X", y="Y", data=data)
    plt.xlabel("AA sub sequence 3 mers", fontdict=font2)
    plt.ylabel("'spike+'=1; 'spike-'=0", fontdict=font2)
    plt.savefig("_spike_fig.png")
    plt.show()

