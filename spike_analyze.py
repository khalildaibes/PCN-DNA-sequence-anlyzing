import os
import sys

from matplotlib import pyplot as plt
import matplotlib
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile, join
import glob
import re
import pickle
from tqdm import tqdm
import csv

import matplotlib

from venn.pyvenn import venn


global init
init = False


def get_overlap(spikes_path, threshold, spikes_dict):
    total_spikes_avg = 0
    temp = spikes_path.split("ID")
    subject_id = temp[1][1:2]
    with open(spikes_path) as spikes_file:
        for line in spikes_file:
            key_value = None
            if "all_spikes.txt" in spikes_path:
                (tri_part, spike_plus_part, spike_minus_part, spike_n_part) = line.split(",")
            else:
                (tri_part, spike_plus_part, spike_minus_part,) = line.split(",")
            (junk, triple_k_mer_sequence) = tri_part.split(":")
            (junk, spike_plus_count) = spike_plus_part.split(":")
            (junk, spike_minus_count) = spike_minus_part.split(":")
            spike_minus_count = spike_minus_count.strip('\n')
            triple_k_mer_sequence_sum = int(spike_plus_count) + int(spike_minus_count)
            if str(triple_k_mer_sequence) in spikes_dict.keys():
                if str(triple_k_mer_sequence) == "SIK":
                    print("here sik")
                if spikes_dict.get(str(triple_k_mer_sequence))[3] == "None":
                    key_value[3] = []

                key_value = spikes_dict.get(str(triple_k_mer_sequence))
                templst = spikes_dict[str(triple_k_mer_sequence)][3]
                templst.append(subject_id) if subject_id not in key_value[3] else key_value[3]
                spikes_dict[str(triple_k_mer_sequence)][0] = key_value[0] + int(triple_k_mer_sequence_sum)
                spikes_dict[str(triple_k_mer_sequence)][1] = key_value[1] + int(spike_plus_count)
                spikes_dict[str(triple_k_mer_sequence)][2] = key_value[2] + int(spike_minus_count)
                spikes_dict[str(triple_k_mer_sequence)][3] = templst
                total_spikes_avg += int(triple_k_mer_sequence_sum)
            else:
                subjects_list = []
                subjects_list.append(subject_id)
                if subject_id == "None":
                    print("here")
                spikes_dict[str(triple_k_mer_sequence)] = [
                    int(triple_k_mer_sequence_sum),
                    int(spike_plus_count),
                    int(spike_minus_count),
                    subjects_list
                ]
                total_spikes_avg += int(triple_k_mer_sequence_sum)


# all_spikes_file = open(all_spikes_path, "r")
# print(all_spikes_file.read())
def save_pie_figure(spikes_path, threshold):
    spikes_dict = {}
    total_spikes_avg = 0
    temp = spikes_path.split("ID")
    subject_id = temp[1][1:2]
    with open(spikes_path) as spikes_file:
        for line in spikes_file:
            try:

                # print(line.split(","))
                (tri_part, spike_plus_part, spike_minus_part) = line.split(",")
                (junk, triple_k_mer_sequence) = tri_part.split(":")
                (junk, spike_plus_count) = spike_plus_part.split(":")
                (junk, spike_minus_count) = spike_minus_part.split(":")
                spike_minus_count = spike_minus_count.strip('\n')
                triple_k_mer_sequence_sum = int(spike_plus_count) + int(spike_minus_count)
                if str(triple_k_mer_sequence) in spikes_dict.keys():
                    key_value = spikes_dict.get(str(triple_k_mer_sequence))
                    spikes_dict[str(triple_k_mer_sequence)] = [
                        key_value[0] + int(triple_k_mer_sequence_sum),
                        key_value[1] + int(spike_plus_count),
                        key_value[2] + int(spike_minus_count),
                        key_value[3].append(subject_id)
                    ]
                    total_spikes_avg += int(triple_k_mer_sequence_sum)
                else:
                    subjects_list = []
                    subjects_list.append(subject_id)
                    spikes_dict[str(triple_k_mer_sequence)] = [
                        int(triple_k_mer_sequence_sum),
                        int(spike_plus_count),
                        int(spike_minus_count),
                        subjects_list
                    ]
                    total_spikes_avg += int(triple_k_mer_sequence_sum)
            except BaseException as ex:
                print("an error", ex, "save_pie_figure", spikes_path)

            # print(triple_k_mer_sequence + " has appeared => " + str(triple_k_mer_sequence_sum) + "times")
            # print(triple_k_mer_sequence + " has appeared in spike + => " + str(spike_plus_count) + "times")
            # print(triple_k_mer_sequence + " has appeared in spike + => " + str(spike_minus_count) + "times")
            # print(triple_k_mer_sequence + " has appeared in spike + => " + str(spike_N_count) + "times")
        try:
            labels = []
            sizes = []

            for seq, spikes_dict_val in spikes_dict.items():
                if spikes_dict_val[0] < total_spikes_avg / len(spikes_dict):
                    pass
                else:
                    label_title = seq + " total " + str(spikes_dict_val[0]) + " spike +" + str(
                        spikes_dict_val[1]) + " spike -" + str(
                        spikes_dict_val[2]) + str(spikes_dict_val[3])
                    labels.append(label_title)
                    sizes.append(spikes_dict_val[0])
            pie_plot = plt.pie(sizes, labels=labels)
            plt.axis('equal')
            # plt.show()
            with open(spikes_path + "_pie_fig" + '.pkl', 'wb')as fid:
                pickle.dump(pie_plot, fid)

            # # Two  lines to make our compiler able to draw:
            plt.savefig(spikes_path + "_pie_fig.png")
            x = np.random.normal(sizes)
            plt.title('spikes appearance distribution ')
            plt.xlabel('appearance times')
            plt.ylabel('how many tri mers')
            plt.hist(sizes)
            with open(spikes_path + "_normal_fig" + '.pkl', 'wb')as fid:
                pickle.dump(x, fid)

            # plt.show()

        except BaseException as ex:
            print("an error", ex, "save_pie_figure", spikes_path)


def save_all_spikes_pie_figure(spikes_path, threshold):
    spikes_dict = {}
    with open(spikes_path) as spikes_file:
        for line in spikes_file:
            (tri_part, spike_plus_part, spike_minus_part, spike_N_part,) = x = line.split(",")
            (junk, triple_k_mer_sequence) = tri_part.split(":")
            (junk, spike_plus_count) = spike_plus_part.split(":")
            (junk, spike_minus_count) = spike_minus_part.split(":")
            (junk, spike_N_count) = spike_N_part.split(":")
            spike_N_count = spike_N_count.strip('\n')
            triple_k_mer_sequence_sum = int(spike_plus_count) + int(spike_minus_count)
            spikes_dict[str(triple_k_mer_sequence)] = [
                int(triple_k_mer_sequence_sum), int(spike_plus_count), int(spike_minus_count), int(spike_N_count)
            ]

            # print(triple_k_mer_sequence + " has appeared => " + str(triple_k_mer_sequence_sum) + "times")
            # print(triple_k_mer_sequence + " has appeared in spike + => " + str(spike_plus_count) + "times")
            # print(triple_k_mer_sequence + " has appeared in spike + => " + str(spike_minus_count) + "times")
            # print(triple_k_mer_sequence + " has appeared in spike + => " + str(spike_N_count) + "times")

    labels = []
    sizes = []

    for seq, spikes_dict_val in spikes_dict.items():
        if spikes_dict_val[0] < threshold:
            pass
        else:
            labels.append(
                seq + " total " + str(spikes_dict_val[0]) + " spike +" + str(spikes_dict_val[1]) + " spike -" + str(
                    spikes_dict_val[2]) + " spike N" + str(spikes_dict_val[3]))
            sizes.append(spikes_dict_val[0])
    pie_plot = plt.pie(sizes, labels=labels)
    # plt.legend()
    plt.title('all spikes pie')
    plt.axis('equal')
    with open(spikes_path + "_pie_fig" + '.pkl', 'wb')as fid:
        pickle.dump(pie_plot, fid)
    # plt.show()
    x = np.random.normal(sizes)
    plt.title('all spikes appearance distribution ')
    plt.xlabel('appearance times')
    plt.ylabel('how many tri mers')
    n, bins, patches = plt.hist(x, 50, density=True, facecolor='g', alpha=0.75)
    with open(spikes_path + "_normal_fig" + '.pkl', 'wb')as fid:
        pickle.dump(x, fid)
    # plt.show()

    # # Two  lines to make our compiler able to draw:
    plt.savefig(spikes_path + "_pie_fig.png")


mypath = "../"
onlyfiles = [f for f in listdir(mypath) if (join(mypath, f))]
txtfiles = []
# for file in glob.glob("*.txt"):
#     txtfiles.append(file)
#     pwd = '2022-08-26-01-04-28-DB-covid_vaccine_new_ID-3_CDR-0_P-All The Sequence_S0_E0_Rall_Results/Amino_acid'
number_of_subjects = 4
directory_list = list()
for root, dirs, files in os.walk("../", topdown=False):
    # for name in files:
    #     print(os.path.join(root, name))
    for name in dirs:
        if name == "Amino_acid":
            directory_list.append(os.path.join(root, name))
            print(os.path.join(root, name))
        # if name == "nucleotides":
        #     directory_list.append(os.path.join(root, name))
        # if name == "pajek":
        #     directory_list.append(os.path.join(root, name))
spike_plus_dictionary = {}
spike_minus_dictionary = {}
spike_n_dictionary = {}
all_spikes_dictionary = {}
print(directory_list)
for name in directory_list:
    pwd = str(name)
    # print(pwd)
    for root, dirs, files in os.walk(name, topdown=False):
        for name in files:

            if name.endswith("all_spikes.txt"):
                print(name)
                all_spikes_path = pwd + '/' + 'all_spikes.txt'
                # print("spike all")
                get_overlap(all_spikes_path, 3000, all_spikes_dictionary)
                save_all_spikes_pie_figure(all_spikes_path, 3000)
            if name.endswith("spike+"):
                # print("spike plus")
                plus_spikes_path = pwd + '/' + 'spike+'
                save_pie_figure(plus_spikes_path, 100)
                get_overlap(plus_spikes_path, 100, spike_plus_dictionary)
            if name.endswith("spike_-.txt"):
                # print("spike mimus")
                plus_spikes_path = pwd + '/' + 'spike_-.txt'
                save_pie_figure(plus_spikes_path, 100)
                get_overlap(plus_spikes_path, 100, spike_minus_dictionary)
            if name.endswith("spike_neutral.txt"):
                # print("spike neutral")
                plus_spikes_path = pwd + '/' + 'spike_neutral.txt'
                save_pie_figure(plus_spikes_path, 1000)
                get_overlap(plus_spikes_path, 1000, spike_n_dictionary)
                with open(plus_spikes_path + "_normal_fig" + '.pkl', 'rb') as fid:
                    ax = pickle.load(fid)
                # plt.show()

# subject_trimers_lists = {3: [], 4: [], 5: [], 6: [], 7: []}
# with open('all_spikes_overlapping.txt', 'a+') as fid:
#     for item in all_spikes_dictionary:
#         key_value = all_spikes_dictionary.get(item)
#         # if len(key_value[3]) >= 4:  # grater than number of subjects
#         fid.write(str(item) + " here is the victor" + str(key_value[3]) + "\n")
#         for subject in key_value[3]:
#             if item not in subject_trimers_lists[subject]:
#                 subject_trimers_lists[subject].append(item)
# with open('spikes_plus_overlapping.txt', 'a+') as fid:
#     for item in spike_plus_dictionary:
#         key_value = spike_plus_dictionary.get(item)
#         # if len(key_value[3]) >= 4:  # grater than number of subjects
#         fid.write(str(item) + " here is the victor" + str(key_value[3]) + "\n")
# with open('spikes_minus_overlapping.txt', 'a+') as fid:
#     for item in spike_minus_dictionary:
#         key_value = spike_minus_dictionary.get(item)
#         # if len(key_value[3]) >= 4:  # grater than number of subjects
#         fid.write(str(item) + " here is the victor" + str(key_value[3]) + "\n")
# with open('spikes_natural_overlapping.txt', 'a+') as fid:
#     for item in spike_n_dictionary:
#         key_value = spike_n_dictionary.get(item)
#         # if len(key_value[3]) >= 4:  # grater than number of subjects
#         fid.write(str(item) + " here is the victor" + str(key_value[3]) + "\n")
# # venn3(subsets=(a), b, a cross b, c, a cross c, b cross c, a cross b cross c )
# plt.clf()
# acrossb = {}
# temp = []
# for akey in spike_plus_dictionary.keys():
#     if akey in spike_minus_dictionary.keys():
#         temp = [acrossb.get(akey, []), spike_minus_dictionary.get(akey), spike_plus_dictionary.get(akey)]
#         acrossb[akey] = temp
# acrossc = {}
#
# temp = []
# for akey in spike_plus_dictionary.keys():
#     if akey in spike_n_dictionary.keys():
#         temp = [acrossc.get(akey, []), spike_n_dictionary.get(akey), spike_plus_dictionary.get(akey)]
#         acrossc[akey] = temp
# bcrossc = {}
#
# temp = []
# for akey in spike_minus_dictionary.keys():
#     if akey in spike_n_dictionary.keys():
#         temp = [bcrossc.get(akey, []), spike_n_dictionary.get(akey), spike_minus_dictionary.get(akey)]
#         bcrossc[akey] = temp
#
# temp = []
# acrossbcrossc = {}
# for akey in spike_minus_dictionary.keys():
#     if akey in spike_n_dictionary.keys() and akey in spike_plus_dictionary.keys():
#         temp = [bcrossc.get(akey, []), spike_n_dictionary.get(akey), spike_minus_dictionary.get(akey),
#                 spike_n_dictionary.get(akey)]
#         bcrossc[akey] = temp
#
# venn3(subsets=(
#     len(spike_plus_dictionary), len(spike_minus_dictionary), len(acrossb), len(spike_n_dictionary), len(acrossc),
#     len(bcrossc), len(acrossbcrossc))
#     , set_labels=('spike_plus_dictionary', 'spike_minus_dictionary', 'spike_n_dictionary'), alpha=0.5)
# plt.show()
subjects_trimers = {'3': [], '4': [], '5': [], '6': [], '7': []}
for akey in all_spikes_dictionary.keys():
    dict_entery = all_spikes_dictionary.get(akey)
    for subject in dict_entery[3]:
        subjects_trimers[subject].append(akey)

labels = venn.get_labels([subjects_trimers['3'], subjects_trimers['4'], subjects_trimers['5']
                             , subjects_trimers['6'], subjects_trimers['7']], fill=['number', 'logic'])
figure, ax = venn.venn5(labels, names=['subject 3', 'subject 4', 'subject 5', 'subject 6', 'subject 7'])
figure.savefig('venn_overlapping.png')

# pwd = '2022-08-26-01-04-28-DB-covid_vaccine_new_ID-3_CDR-0_P-All The Sequence_S0_E0_Rall_Results/Amino_acid'

# all_spikes_path = pwd + '/' + 'all_spikes'
# save_all_spikes_pie_figure(all_spikes_path, 3000)
# plus_spikes_path = pwd + '/' + 'spike+'
# save_pie_figure(plus_spikes_path, 100)

# plus_spikes_path = pwd + '/' + 'spike_-.txt'
# save_pie_figure(plus_spikes_path, 100)
# plus_spikes_path = pwd + '/' + 'spike_neutral.txt'
# save_pie_figure(plus_spikes_path, 1000)
