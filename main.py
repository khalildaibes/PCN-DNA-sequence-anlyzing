# -*- coding: utf-8 -*-

# NEEDED IMPORT AND LIBRARIES
# A LIBRARY YTO WORK WITH CSV FILES (READ\WRITE)
import csv
# A LIBRARY TO WORK WITH HTML FILES (READ\WRITE)
import numbers
from datetime import datetime
import pandas as pd
from bs4 import BeautifulSoup
# A LIBRARY TO WORK VIEW PROGRESS IN CONSOLE WINDOW
from tqdm import tqdm
# A LIBRARY TO WORK WIT THE NETWORK X AND ANALYZE THE NETWORK ALSO GET DATA ABOUT A NETWORK
import networkx as nx
# A LIBRARY TO WORK WITH DATABASE AND CONNECT TO THE CONNECTION MODULE
import db_connection as db_connection
# A LIBRARY TO WORK VIEW PROGRESS IN CONSOLE WINDOW
import pickle
# A LIBRARY TO WORK WITH THREADS AND PROCESS
from threading import Thread
# A LIBRARY TO WORK WITH THE OPERATION SYSTEM
import os
from matplotlib import pyplot as plt


#####################################################################################################################
# @DESCRIPTION : DECLARING AND CREATOING THE DIORCTORES NEDDED
# @INPUT : DB : DATABASE NAME , ID: SUBJECT ID
# @OUTPUT : ---
#####################################################################################################################
def create_files(DB, ID, cdr_length, selected_analyzed_position, startposition, endposition, rows):
    # CREATING THE DIRECTORY
    global directory
    directory = "{}-DB-{}_ID-{}_CDR-{}_P-{}_S{}_E{}_R{}_Results".format(datetime.today().strftime('%Y-%m-%d-%H-%M-%S'),
                                                                        DB, ID, cdr_length, selected_analyzed_position,
                                                                        startposition, endposition, rows)
    try:
        # CREATE A DIRECTORY
        os.makedirs(directory)
        print("Directory '%s' created successfully" % directory)
        # CATCH EXCEPTION
    except OSError as error:
        print("Directory '%s' can not be created" % directory)

    try:  # DECLARE MUST HAVE FOLDERS AND FILES
        global toFile
        toFile = "{}\{}emp_id{}.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)
        global toFile1
        toFile1 = "{}\{}emp_AA_id{}.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)
        global toFile2
        toFile2 = "{}\{}emp_id{}_seqK.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)
        global u_path
        u_path = "{}\{}_id{}_uniqeKmers.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)
        global Amino_acid
        Amino_acid = "{}\Amino_acid".format(str(os.getcwd() + "\\" + directory))
        isdir = os.path.isdir(Amino_acid)
        if not isdir:
            os.makedirs(Amino_acid)
        global nucleotides
        nucleotides = "{}\\nucleotides".format(str(os.getcwd() + "\\" + directory))
        isdir = os.path.isdir(nucleotides)
        if not isdir:
            os.makedirs(nucleotides)
        global pajek
        pajek = "{}\pajek".format(str(os.getcwd() + "\\" + directory))
        isdir = os.path.isdir(pajek)
        if not isdir:
            os.makedirs(pajek)
    except OSError as error:
        print("Directory '%s' can not be created" % directory)


# A CLASS TO HOLD THE NUCLEOTIDES NODES DATA AND METADATA
class NucleotidesNetNodes:
    def __init__(self, sub_seq=None, sequence=None, subject=None, seq_id=None, position=None):
        self.sub_seq = sub_seq
        self.sequence = sequence
        self.subject = subject
        self.seq_id = seq_id
        self.position = position


# A CLASS TO HOLD THE AMINO ACID NODES DATA AND METADATA
class AminoAcidNetNodes:
    def __init__(self, sequence=None, subject=None):
        self.sequence = sequence
        self.subject = subject


#####################################################################################################################
# @DESCRIPTION : DECLARING THE REQUESTED DATA QUERY BASED ON THE INPUT VALUES
# @INPUT : DB: DATABASE , ID: SUBJECT ID , cdr_length: CDR LENGTH
# @OUTPUT : CMD QUERY
#####################################################################################################################
def select_query(db, id, cdr_length=14):
    cmd = ""
    cdr_length = 14
    # SPECIFIC SUBJECTS AND SPECIFIC CDR LENGTHS
    if id != "all" and cdr_length != "not specific":
        cmd = "select seq.*, coll.*,sammet.* from {}.sequences as seq inner join {}.sequence_collapse as coll on seq.ai=coll.seq_ai inner join {}.sample_metadata as sammet on sammet.sample_id=seq.sample_id WHERE seq.subject_id={} AND seq.functional=1 AND coll.instances_in_subject !=0 AND coll.copy_number_in_subject > 1  AND seq.deletions is null AND  seq.insertions is null AND LENGTH(seq.cdr3_aa)={} AND (sammet.value like'%spike+%' OR sammet.value like'%spike-%') ORDER BY RAND()".format(
            db, db, db, id, cdr_length)
    # ALL SUBJECTS AND SPECIFIC CDR LENGTHS
    elif id == "all" and cdr_length != "not specific":
        cmd = "select seq.*, coll.*,sammet.* from {}.sequences as seq inner join {}.sequence_collapse as coll on seq.ai=coll.seq_ai inner join {}.sample_metadata as sammet on sammet.sample_id=seq.sample_id WHERE seq.functional=1 AND coll.instances_in_subject !=0 AND coll.copy_number_in_subject > 1  AND seq.deletions is null AND  seq.insertions is null AND LENGTH(seq.cdr3_aa)={} ORDER BY RAND() AND(sammet.value like'%spike+%' OR sammet.value like'%spike-%') ORDER BY RAND()".format(
            db, db, db, cdr_length)
    # SPECIFIC SUBJECTS AND ALL CDR LENGTHS
    elif id != "all" and cdr_length == "not specific":
        cmd = "select seq.*, coll.*,sammet.* from {}.sequences as seq inner join {}.sequence_collapse as coll on seq.ai=coll.seq_ai inner join {}.sample_metadata as sammet on sammet.sample_id=seq.sample_id WHERE seq.subject_id={} AND seq.functional=1 AND coll.instances_in_subject !=0 AND coll.copy_number_in_subject > 1  AND seq.deletions is null AND seq.insertions is null AND LENGTH(seq.cdr3_aa)>8 AND  (sammet.value like'%spike+%' OR sammet.value like'%spike-%') ORDER BY RAND()".format(
            db, db, db, id)
    # ALL SUBJECTS AND ALL CDR LENGTHS
    elif id == "all" and cdr_length == "not specific":
        cmd = "select seq.*, coll.*,sammet.* from {}.sequences as seq inner join {}.sequence_collapse as coll on seq.ai=coll.seq_ai inner join {}.sample_metadata as sammet on sammet.sample_id=seq.sample_id WHERE seq.functional=1 AND coll.instances_in_subject !=0 AND coll.copy_number_in_subject > 1  AND seq.deletions is null AND  seq.insertions is null AND LENGTH(seq.cdr3_aa)>8  AND (sammet.value like'%spike+%' OR sammet.value like'%spike-%') ORDER BY RAND()".format(
            db, db, db, id)

    return cmd


#####################################################################################################################
# @DESCRIPTION : GET ROWS BASED ON THE SELECTED ROWS COUNT
# @INPUT : amount: ROWS COUNT  , mycursor: DATABASE EXECUTE OBJECT
# @OUTPUT : A list that holds dicts as rows
#####################################################################################################################
def fetch_based_on_row_count(amount, mycursor):
    if 'all' in str(amount):
        seq = mycursor.fetchall()  # select how many rows
    else:
        seq = mycursor.fetchmany(int(amount))  # select how many rows
    return seq


#####################################################################################################################
# @DESCRIPTION : replace each AA with its most frequented nucleotide way
# @INPUT : DB: DATABASE NAME, ID, cdr_length, selected_analyzed_position, startposition, endposition, rows
# @OUTPUT : A list that holds dicts as rows
#####################################################################################################################
def replace_nucleotide(DB, ID, selected_analyzed_position, startposition, endposition, rows, cdr_length=14):
    print("\nFirst step (translate sequences to amino acids):")
    print("First step strats...\n")
    # CREATING A CONNECTION WITH THE DATABASE
    mydb = db_connection.create_connection()
    from csv import DictReader
    mycursor = mydb.cursor(dictionary=True)
    # DNA codon table
    protein = {"TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
               "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
               "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
               "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",
               "TCT": "s", "CCT": "P", "ACT": "T", "GCT": "A",
               "TCC": "s", "CCC": "P", "ACC": "T", "GCC": "A",
               "TCA": "s", "CCA": "P", "ACA": "T", "GCA": "A",
               "TCG": "s", "CCG": "P", "ACG": "T", "GCG": "A",
               "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
               "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
               "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
               "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",
               "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
               "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
               "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
               "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G",
               "---": "-"
               }
    keys_use = {}
    for key in aa_most_frequent.keys():
        lst = aa_most_frequent[key]
        max_index = 1
        max_value = 0
        for value, index in enumerate(lst):
            if isinstance(value, numbers.Number):
                if value > max_value:
                    max_value = value
                    max_index = index

        keys_use[key] = [lst[0], max_index, max_value]
    # Command for getting the sequences translated:

    command = select_query(DB, ID, cdr_length)
    mycursor.execute(command)
    global Selected_Number_Of_Rows
    kmers_G = nx.Graph(name='kmers_of_position')
    seq = fetch_based_on_row_count(rows, mycursor)
    print("Translating starts ....")
    with open(toFile, 'r+') as seq:
        csv_dict_reader = DictReader(seq)
        # iterate over each line as a ordered dictionary
        for line in tqdm(csv_dict_reader):
            # Check file as empty

            dna = ""
            germ = ""
            protein_sequence = ""
            # fix the germline to match the cdr with N's

            # -------------- To NNNNNNNNNNNNNN #
            dna = line["sequence"]
            seqID = line["seq_id"]
            sampleID = line["sample_id"]
            subjectID = line["subject_id"]
            spike = line["value"]
            # Generate protein sequence
            j = 0
            temp_pos = 0
            start = 0
            end = 0
            triple = 0
            nuclotides_sub_seq = ""
            complete_cdr_length = cdr_length * 3
            startswitcher = {
                'Select Region': 112,
                'All The Sequence': 0,
                'Select Start And End Position': int(startposition),
            }
            endswitcher = {
                'Select Region': 112 + complete_cdr_length,
                'All The Sequence': len(dna) - ((len(dna) % 3) * 3),
                'Select Start And End Position': int(int(endposition) - (int(endposition) % 3) * 3)
            }

            startswitcher = int(startswitcher.get(selected_analyzed_position))
            endswitcher = int(endswitcher.get(selected_analyzed_position))
            tripleAA = ""
            for i in range(startswitcher, endswitcher, 3):

                if dna[i] == "N" and dna[i + 1] == "N" and dna[i + 2] == "N":
                    protein_sequence += "x"
                elif dna[i] == "N" or dna[i + 1] == "N" or dna[i + 2] == "N" and dna[i] != "-" and \
                        dna[i + 1] != "-" and \
                        dna[i + 2] != "-":
                    protein_sequence += "x"
                elif dna[i] == "-" and dna[i + 1] == "-" and dna[i + 2] == "-":
                    protein_sequence += protein[dna[i:i + 3]]
                elif dna[i] == "-" or dna[i + 1] == "-" or dna[i + 2] == "-":
                    i = i + 0
                else:
                    protein_sequence += protein[dna[i:i + 3]]

                    tripleAA += protein[dna[i:i + 3]]
                    nuclotides_use = keys_use.get(protein[dna[i:i + 3]])
                    nuclotides_sub_seq += nuclotides_use[0]

                    triple += 1
                    most_freq_lst = aa_most_frequent.get(protein[dna[i:i + 3]])
                    index = most_freq_lst.index(dna[i:i + 3])
                    most_freq_lst[index + 1] += 1
                    aa_most_frequent[protein[dna[i:i + 3]]] = most_freq_lst

                    if triple == 3:

                        if protein[dna[i:i + 3]] != "-" and protein[dna[i:i + 3]] != "*":
                            a_sequence = line["sequence"]  # protein a node
                            a_subject = line["subject_id"]
                            a_seq_id = line["seq_id"]
                            # protein a node

                            a_node = NucleotidesNetNodes(tripleAA, a_sequence, a_subject, a_seq_id, i)
                            old_list = amino_acids_kmers_nets_G.get(i, None)
                            if old_list:
                                old_list.append(a_node)
                            else:
                                old_list = []
                                old_list.append(a_node)
                            amino_acids_kmers_nets_G[i] = old_list
                        tripleAA = ""

                        triple = 0
                        a_sequence = line["sequence"]  # protein a node
                        a_subject = line["subject_id"]
                        a_position = i
                        a_seq_id = line["seq_id"]
                        # protein a node
                        a_node = NucleotidesNetNodes(nuclotides_sub_seq, a_sequence, a_subject, a_seq_id, a_position)
                        aa_count[dna[i:i + 3]][0] += 1
                        aa_count[dna[i:i + 3]].append(a_node)
                        old_list = replaced_nucoltides_kmers_nets.get(i, [])
                        old_list.append(a_node)
                        replaced_nucoltides_kmers_nets[i] = old_list
                        nuclotides_sub_seq = ""
                j = j + 3

                temp_pos = temp_pos + 9

        print("Translating DONE!")
        save_fasta = open(str(DB) + str(ID) + str(selected_analyzed_position) + str(startposition)
                          + str(endposition) + str(rows) + str(cdr_length) + "_replaced_cdr3_fasta.fasta", "w+")
        save_fasta.write(protein_sequence)
        save_fasta.close()


def startf(DB, ID, selected_analyzed_position, startposition, endposition, rows, cdr_length=14):
    create_files(DB, ID, cdr_length, selected_analyzed_position, startposition, endposition, rows)
    print("\nFirst step (translate sequences to amino acids):")
    print("First step strats...\n")
    sequence_legth = 231
    removed = "trash.csv"

    mydb = db_connection.create_connection()
    mycursor = mydb.cursor(dictionary=True)
    # DNA codon table
    protein = {"TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
               "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
               "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
               "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",
               "TCT": "s", "CCT": "P", "ACT": "T", "GCT": "A",
               "TCC": "s", "CCC": "P", "ACC": "T", "GCC": "A",
               "TCA": "s", "CCA": "P", "ACA": "T", "GCA": "A",
               "TCG": "s", "CCG": "P", "ACG": "T", "GCG": "A",
               "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
               "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
               "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
               "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",
               "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
               "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
               "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
               "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G",
               "---": "-"
               }

    for key in protein.keys():

        if protein.get(key) in aa_most_frequent.keys():

            if key not in aa_most_frequent.get(protein.get(key)):
                old_list = aa_most_frequent[protein.get(key)]
                old_list.append(key)
                old_list.append(0)
                aa_most_frequent[protein.get(key)] = old_list
        else:
            lst = []
            lst.append(key)
            lst.append(0)
            aa_most_frequent[protein.get(key)] = lst

    # Command for getting the sequences translated:

    command = select_query(DB, ID, cdr_length)
    mycursor.execute(command)
    global Selected_Number_Of_Rows
    count = 0
    kmers_G = nx.Graph(name='kmers_of_position')
    seq = fetch_based_on_row_count(rows, mycursor)
    print("Translating starts ....")
    save_fasta = open(str(DB) + str(ID) + str(selected_analyzed_position) + str(startposition)
                      + str(endposition) + str(rows) + str(cdr_length) + "cdr3_fasta.fasta", "w+")
    save_fasta.write(">seq" + str(count))
    with open(toFile, 'w', newline='') as new_file:
        csv_writer = csv.writer(new_file)
        csv_writer.writerow(
            ['seq_id', 'sequence', 'TranslatedSeq', 'ai', 'subject_id', 'clone_id', 'sample_id',
             'cdr3_aa'])
        for line in tqdm(seq):
            dna = ""
            germ = ""
            protein_sequence = ""
            # fix the germline to match the cdr with N's
            germ = line['germline']
            cdr3_seq = line['cdr3_aa']
            spike = line["value"]
            cdr3Length = line['cdr3_num_nts']
            postCDR = line['post_cdr3_length']
            x = int(cdr3Length) + int(postCDR)
            replaced = germ[-x:]
            replaced = replaced.replace('-', 'N', x)
            germ = germ.replace(germ[-x:], replaced)
            # -------------- To NNNNNNNNNNNNNN #
            dna = line['sequence']
            seqID = line['seq_id']
            ai = line['ai']
            cloneID = line['clone_id']
            sampleID = line['sample_id']
            subjectID = line['subject_id']
            # Generate protein sequence
            j = 0
            temp_pos = 0
            start = 0
            end = 0
            triple = 0
            nuclotides_sub_seq = ""
            complete_cdr_length = cdr_length * 3
            startswitcher = {
                'Select Region': 112,
                'All The Sequence': 0,
                'Select Start And End Position': int(startposition),
            }
            endswitcher = {
                'Select Region': 112 + complete_cdr_length,
                'All The Sequence': len(dna) - ((len(dna) % 3) * 3),
                'Select Start And End Position': int(int(endposition) - (int(endposition) % 3) * 3)
            }

            startswitcher = int(startswitcher.get(selected_analyzed_position))
            endswitcher = int(endswitcher.get(selected_analyzed_position))
            tripleAA = ""

            for i in range(startswitcher, endswitcher, 3):

                if dna[i] == "N" and dna[i + 1] == "N" and dna[i + 2] == "N":
                    protein_sequence += "x"
                elif dna[i] == "N" or dna[i + 1] == "N" or dna[i + 2] == "N" and dna[i] != "-" and \
                        dna[i + 1] != "-" and \
                        dna[i + 2] != "-":
                    protein_sequence += "x"
                elif dna[i] == "-" and dna[i + 1] == "-" and dna[i + 2] == "-":
                    protein_sequence += protein[dna[i:i + 3]]
                elif dna[i] == "-" or dna[i + 1] == "-" or dna[i + 2] == "-":
                    i = i + 0
                else:
                    protein_sequence += protein[dna[i:i + 3]]
                    tripleAA += protein[dna[i:i + 3]]
                    nuclotides_sub_seq += dna[i:i + 3]
                    triple += 1

                    most_freq_lst = aa_most_frequent.get(protein[dna[i:i + 3]])
                    index = most_freq_lst.index(dna[i:i + 3])
                    most_freq_lst[index + 1] += 1
                    aa_most_frequent[protein[dna[i:i + 3]]] = most_freq_lst

                    if triple == 3:
                        if tripleAA not in spike_dic.keys():
                            spkie_count = [0, 0]
                            spike_dic[tripleAA] = spkie_count
                        if "spike+" in str(spike):
                            x = spike_dic.get(tripleAA)
                            x[0] += 1
                            spike_dic[tripleAA] = x
                        elif "spike-" in str(spike):
                            x = spike_dic.get(tripleAA)
                            x[1] += 1
                            spike_dic[tripleAA] = x

                        if protein[dna[i:i + 3]] != "-" and protein[dna[i:i + 3]] != "*":
                            a_sequence = line['sequence']  # protein a node
                            a_subject = line['subject_id']
                            a_seq_id = line['seq_id']
                            # protein a node

                            a_node = NucleotidesNetNodes(tripleAA, a_sequence, a_subject, a_seq_id, i)
                            old_list = amino_acids_kmers_nets_G.get(i, None)
                            if old_list:
                                old_list.append(a_node)
                            else:
                                old_list = []
                                old_list.append(a_node)
                            amino_acids_kmers_nets_G[i] = old_list
                        tripleAA = ""

                        triple = 0
                        a_sequence = line['sequence']  # protein a node
                        a_subject = line['subject_id']
                        a_position = i
                        a_seq_id = line['seq_id']
                        # protein a node
                        a_node = NucleotidesNetNodes(nuclotides_sub_seq, a_sequence, a_subject, a_seq_id, a_position)
                        aa_count[dna[i:i + 3]][0] += 1
                        aa_count[dna[i:i + 3]].append(a_node)
                        old_list = nucoltides_kmers_nets.get(i, [])
                        old_list.append(a_node)
                        nucoltides_kmers_nets[i] = old_list
                        nuclotides_sub_seq = ""
                j = j + 3

                temp_pos = temp_pos + 9

            csv_writer.writerow(
                [seqID, dna, protein_sequence, ai, subjectID, cloneID, sampleID, cdr3_seq, spike])
            print("Translating DONE!")
            save_fasta.write(dna)
            count += 1
            save_fasta.write("\n>seq" + str(count) + "\n")
        save_fasta.close()


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


def calculateHammingDistance(protein1, protein2):
    count = 0
    rangeofseq = min(len(protein1), len(protein2))
    for j in range(0, rangeofseq):
        if protein1[j] == protein2[j]:
            pass
        else:
            count = count + 1
    if count > 2:
        return 0

    return count


def calculateAminoAcidsSimilarty(protein1, protein2):
    count = 0
    rangeofseq = min(len(protein1), len(protein2))
    for j in range(0, rangeofseq):
        if protein1[j] == protein2[j]:
            pass
        else:
            count = count + 1
    if count > 1:
        return 0
    return count


def show_Nuclotides_netwrok_stats(networkx_graph, network_name):
    max = 0
    the_key = None
    result = ""
    sequences = ""
    for key in aa_count:
        result += "\n"
        result += (' ' + key + '  ')
        lst = aa_count.get(key)
        result += (str(lst[0]) + '  ' + ' '.join(str(x.sub_seq) for x in lst[1:]))
        if lst[0] > max:
            max = lst[0]
            sequences = aa_count.get(key)
            the_key = key
    result += "\n"
    result += (str(the_key) + '   The max encountered sub sequence is ' + str(max) + '  ' + ' '.join(
        str(x.sub_seq) + ' ' + str(x.position) for x in sequences[1:]))
    most_connected_node_counter = [0, None]
    for node, node_attrs in tqdm(networkx_graph.nodes(data=True)):
        if networkx_graph.degree(node) > most_connected_node_counter[0]:
            most_connected_node_counter[1] = node
            most_connected_node_counter[0] = networkx_graph.degree(node)
    result += "\n"
    result += str(most_connected_node_counter[1]) + ' with this count of connectivity' + str(
        most_connected_node_counter[0])
    with open(str(network_name) + "network_stats.txt", "w") as file:
        file.write(str(result))
    return result


def draw_graph3(networkx_graph, notebook=True, output_filename='empgraph.html', show_buttons=False,
                only_physics_buttons=False):
    from pyvis import network as net
    actual_total_appearances = 0
    real_nodes_count = 0
    real_edges_count = 0
    nuclotide_count_dict = {}
    # make a pyvis network
    pyvis_graph = net.Network(height='500px', width='500px', notebook=notebook)
    pyvis_graph.width = '1000px'
    pyvis_graph = net.Network(height='900px', width='900px', bgcolor='#222222', font_color='white')
    switcher = {
        '1': 'red',
        '2': 'green',
        '3': 'red',
        '4': 'magenta',
        '5': 'cyan',
        '6': 'blue',
        '7': 'black',
        '8': '#eeefff'
    }
    # for each node and its attributes in the networkx graph
    for node, node_attrs in tqdm(networkx_graph.nodes(data=True)):
        color1 = switcher.get(str(node.subject))
        x = nuclotide_count_dict.get(node.sub_seq)
        if x:
            title = "The node Sub sequence is " + str(node.sub_seq)
            if x[0] >= 2:
                x[0] += 1
                nuclotide_count_dict[node.sub_seq] = x
                x.append(node)

                for node1 in x[1:]:
                    title += str(x[0]) + " the subject id is : " + str(
                        node1.subject) + " The sequence id is :" + str(node1.seq_id) + "\n"
                if node.sub_seq not in pyvis_graph.nodes:
                    pyvis_graph.add_node(node.sub_seq, title=title,
                                         label=str(str(node.sub_seq)), color=color1, )
                    real_nodes_count += 1
            else:
                x[0] += 1
                nuclotide_count_dict[node.sub_seq] = x
                x.append(node)
                if x[0] >= 2:
                    for node1 in x[1:]:
                        title += str(x[0]) + " the subject id is : " + str(
                            node1.subject) + " The sequence id is :" + str(node1.seq_id) + "\n"
                    if node.sub_seq not in pyvis_graph.nodes:
                        pyvis_graph.add_node(node.sub_seq, title=title,
                                             label=str(str(node.sub_seq)), color=color1, )
                        real_nodes_count += 1
                        actual_total_appearances += x[0]
        else:
            x = [0]
            x[0] = 1
            x.append(node)
            nuclotide_count_dict[node.sub_seq] = x
    #         print(node,node_attrs)

    # for each edge and its attributes in the networkx graph
    for source, target, edge_attrs in tqdm(networkx_graph.edges(data=True)):
        # if value/width not specified directly, and weight is specified, set 'value' to 'weight'
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            # place at key 'value' the weight of the edge
            edge_attrs['value'] = edge_attrs['weight']
        # add the edge
        x = nuclotide_count_dict.get(source.sub_seq)
        y = nuclotide_count_dict.get(target.sub_seq)
        if x:
            if x[0] >= 2:
                if y:
                    if y[0] >= 2:
                        pyvis_graph.add_edge(source.sub_seq, target.sub_seq, title=str(edge_attrs['weight']),
                                             value=int(edge_attrs['weight']), weight=int(edge_attrs['weight'])
                                             , color='yellow')
                        real_edges_count += 1

    # turn buttons on
    if show_buttons:
        if only_physics_buttons:
            pyvis_graph.show_buttons(filter_=["physics"])
        else:
            pyvis_graph.show_buttons()

    # return and also save
    pyvis_graph.save_graph(output_filename + ".html")
    soup = BeautifulSoup(open(output_filename + ".html"), 'html.parser')
    extra_html = '''
    <div class="record">
        <div class="header">
        <h2>
        '''
    extra_html += "real nodes count =" + str(real_nodes_count) + "\n"
    extra_html += "actual total appearances =" + str(actual_total_appearances) + "\n"
    extra_html += "real edges count =" + str(real_edges_count) + "\n"
    extra_html += "</h2> <div class='title'>"
    extra_html += nx.info(networkx_graph)
    #######################################
    extra_html += "Network density:" + str(nx.density(networkx_graph)) + "\n"
    extra_html += ""
    # plt.plot(list(nuclotide_count_dict.keys()), [nuclotide_count_dict.get(x)[0] for x in nuclotide_count_dict.keys()])
    # plt.savefig("tempfig.png")
    # plt.show()

    extra_html += '''
            </div>
        </div>
    </div>'''

    div = soup.select("#config")
    soup.find_all('div', {"id": "mynetwork"})[-1].insert_after(BeautifulSoup(extra_html, 'html.parser'))
    with open(output_filename + ".html", "w") as file:
        file.write(str(soup))


def export_preproccessed_data():
    pickle.dump(nucoltides_kmers_nets, open('nucoltides_kmers_nets.pkl', 'wb'))
    pickle.dump(aa_count, open('protein_count.pkl', 'wb'))


def readTrueData(name):
    fName = str(name + '.pkl')
    with open(fName, 'rb') as f:  # with statement avoids file leak
        # Match load with file object, and provide encoding for Py2 str
        return pickle.load(f)


def import_preproccessed_data():
    global nucoltides_kmers_nets
    nucoltides_kmers_nets = readTrueData('nucoltides_kmers_nets')
    global aa_count
    aa_count = readTrueData('protein_count')


def draw_kmers_graph3(networkx_graph, notebook=True, output_filename='empgraph', show_buttons=False,
                      only_physics_buttons=False):
    # import
    from pyvis import network as net
    real_nodes_count = 0
    real_edges_count = 0
    # make a pyvis network
    pyvis_graph = net.Network(height='900px', width='900px', bgcolor='#222222', font_color='white')
    switcher = {
        '1': 'red',
        '2': 'green',
        '3': 'red',
        '4': 'magenta',
        '5': 'cyan',
        '6': 'blue',
        '7': 'black',
        '8': '#eeefff'
    }
    nuclotide_count_dict = {}
    # for each node and its attributes in the networkx graph
    nodes_count = len(networkx_graph.nodes(data=True))
    for node, node_attrs in tqdm(networkx_graph.nodes(data=True)):
        color1 = switcher.get(str(node.subject))
        x = nuclotide_count_dict.get(node.sub_seq)
        if x:

            title = "The node Sub sequence is " + str(node.sub_seq)
            if x[0] >= 2:
                x[0] += 1
                nuclotide_count_dict[node.sub_seq] = x
                x.append(node)

                for node1 in x[1:]:
                    title += " the subject id is : " + str(
                        node1.subject) + " The sequence id is :" + str(node1.seq_id) + "\n"
                pyvis_graph.add_node(node.sub_seq, title=title, label=str(str(node.sub_seq)), color=color1, )
            else:
                x[0] += 1
                nuclotide_count_dict[node.sub_seq] = x
                x.append(node)
                if x[0] >= 2:
                    for node1 in x[1:]:
                        title += " the subject id is : " + str(
                            node1.subject) + " The sequence id is :" + str(node1.seq_id) + "\n"
                    if node.sub_seq not in pyvis_graph.nodes:
                        pyvis_graph.add_node(node.sub_seq, title=title,
                                             label=str(str(node.sub_seq)), color=color1, )
                        real_nodes_count += 1

        else:
            x = [0]
            x[0] = 1
            x.append(node)
            nuclotide_count_dict[node.sub_seq] = x

    #         print(node,node_attrs)

    # for each edge and its attributes in the networkx graph
    for source, target, edge_attrs in tqdm(networkx_graph.edges(data=True)):

        # if value/width not specified directly, and weight is specified, set 'value' to 'weight'
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            # place at key 'value' the weight of the edge
            edge_attrs['value'] = edge_attrs['weight']
        # add the edge
        title = ""
        x = nuclotide_count_dict.get(source.sub_seq)
        y = nuclotide_count_dict.get(target.sub_seq)
        if x:
            if x[0] >= 2:
                if y:
                    if y[0] >= 2:
                        pyvis_graph.add_edge(source.sub_seq, target.sub_seq, title=str(3 - int(edge_attrs['weight'])),
                                             value=3 - int(edge_attrs['weight']), weight=3 - int(edge_attrs['weight']))
                        real_edges_count += 1

        # turn buttons on
    if show_buttons:
        if only_physics_buttons:
            pyvis_graph.show_buttons(filter_=["physics"])
        else:
            pyvis_graph.show_buttons()
    # nx.number_weakly_connected_components(G)
    if nx.number_connected_components(networkx_graph):
        nx.number_connected_components(networkx_graph)

    pyvis_graph.save_graph(output_filename + ".html")
    soup = BeautifulSoup(open(output_filename + ".html"), 'html.parser')
    extra_html = '''
       <div class="record">
           <div class="header">
           <h2>
           '''
    extra_html += "real nodes count =" + str(real_nodes_count) + "\n"
    extra_html += "real edges count =" + str(real_edges_count) + "\n"
    extra_html += "</h2> <div class='title'>"
    extra_html += nx.info(networkx_graph)
    #######################################
    extra_html += "Network density:" + str(nx.density(networkx_graph)) + "\n"
    extra_html += "" + plt.plot(nuclotide_count_dict.keys(), nuclotide_count_dict.values)

    extra_html += '''
               </div>
           </div>
       </div>'''
    div = soup.select("#config")
    soup.find_all('div', {"id": "mynetwork"})[-1].insert_after(BeautifulSoup(extra_html, 'html.parser'))
    with open(output_filename + ".html", "w") as file:
        file.write(str(soup))


def compare_nucoltides(node1, nodes_list, nucoltides_kmers_nets_G, key):
    for node2 in nodes_list:
        if node1.sub_seq != node2.sub_seq and node1.sub_seq[0:2] in node2.sub_seq or node1.sub_seq[3:5] in \
                node2.sub_seq or node1.sub_seq[6:8] in node2.sub_seq:
            w = calculateHammingDistance(node1.sub_seq, node2.sub_seq)
            if w:
                nucoltides_kmers_nets_G.add_edge(node1, node2, weight=w)


def create_amino_acid_netwrok():
    os.chdir(str(directory + "\\") + 'Amino_acid')
    AA_G = nx.Graph(name='AA Network')
    for key in amino_acids_kmers_nets_G.keys():
        AA_G.clear()
        for a_node in amino_acids_kmers_nets_G.get(key):
            for b_node in amino_acids_kmers_nets_G.get(key):
                if calculateAminoAcidsSimilarty(a_node.sub_seq, b_node.sub_seq):
                    AA_G.add_edge(a_node, b_node, weight=calculateAminoAcidsSimilarty(a_node.sub_seq, b_node.sub_seq))

        draw_graph3(AA_G, output_filename=str(key) + '_AA_G_graph_output',
                    show_buttons=True,
                    notebook=False)
        show_Nuclotides_netwrok_stats(AA_G, str(key) + '_NC_G_graph_output')
    plot_spikes(str(directory))
    os.chdir('..')



def create_nuclotides_netwrok():
    print('start the nucoltides node positions')
    nucoltides_kmers_nets_G = {}
    print(len(nucoltides_kmers_nets))
    for key in nucoltides_kmers_nets:
        threads = []
        nucoltides_kmers_nets_G.clear()
        nucoltides_kmers_nets_G = nx.Graph(name='nuclotides kmers_of_position ' + str(key))
        nodes_list = nucoltides_kmers_nets.get(key)
        for node1 in nodes_list:
            thread = Thread(target=compare_nucoltides, args=(node1, nodes_list, nucoltides_kmers_nets_G, key))
            threads.append(thread)
        print("starting with the start")
        for x in tqdm(threads):
            x.start()
        print("starting with the join")
        for x in tqdm(threads):
            x.join()
        print("printing network")
        os.chdir(pajek)
        nx.write_pajek(nucoltides_kmers_nets_G, str(key) + '_NC_G_graph_output1.net')
        os.chdir('..')
        make_netwrok_labels(nucoltides_kmers_nets_G, key)
        os.chdir(nucleotides)
        draw_kmers_graph3(nucoltides_kmers_nets_G,
                          output_filename=str(key) + '_NC_G_graph_output',
                          show_buttons=True,
                          notebook=False)

        show_Nuclotides_netwrok_stats(nucoltides_kmers_nets_G, str(key) + '_NC_G_graph_output')
        os.chdir('..')


def create_replaced_nuclotides_netwrok():
    print('start the  replaced nucoltides node positions')
    replaced_nucoltides_kmers_nets_G = {}
    print(len(replaced_nucoltides_kmers_nets))
    for key in replaced_nucoltides_kmers_nets:
        threads = []
        replaced_nucoltides_kmers_nets_G.clear()
        replaced_nucoltides_kmers_nets_G = nx.Graph(name='replaced nuclotides kmers_of_position ' + str(key))
        nodes_list = replaced_nucoltides_kmers_nets.get(key)
        for node1 in nodes_list:
            thread = Thread(target=compare_nucoltides, args=(node1, nodes_list, replaced_nucoltides_kmers_nets_G, key))
            threads.append(thread)
        print("starting with the start")
        for x in tqdm(threads):
            x.start()
        print("starting with the join")
        for x in tqdm(threads):
            x.join()
        print("printing network")
        os.chdir(pajek)
        nx.write_pajek(replaced_nucoltides_kmers_nets_G,
                       str(key) + '_replaced_NC_G_graph_output1.net')
        os.chdir('..')
        make_netwrok_labels(replaced_nucoltides_kmers_nets_G, key)
        os.chdir(nucleotides)
        draw_kmers_graph3(replaced_nucoltides_kmers_nets_G,
                          output_filename=str(key) + '_replaced_NC_G_graph_output.html',
                          show_buttons=True,
                          notebook=False)

        show_Nuclotides_netwrok_stats(replaced_nucoltides_kmers_nets_G, str(key) + '_replaced_NC_G_graph_output')
        os.chdir('..')


def make_netwrok_labels(networkx_graph, key):
    G_labels = nx.Graph(name='AA Netwrok')
    for node, node_attrs in tqdm(networkx_graph.nodes(data=True)):
        G_labels.add_node(node.sub_seq, title=str(str(node.sub_seq) + str(node.subject)),
                          label=str(str(node.sub_seq) + str(node.subject)))
    #         print(node,node_attrs)
    # for each edge and its attributes in the networkx graph
    for source, target, edge_attrs in tqdm(networkx_graph.edges(data=True)):
        # if value/width not specified directly, and weight is specified, set 'value' to 'weight'
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            # place at key 'value' the weight of the edge
            edge_attrs['value'] = edge_attrs['weight']
            G_labels.add_edge(source.sub_seq, target.sub_seq, weight=int(edge_attrs['weight']),
                              title=str(edge_attrs['weight']))
    os.chdir(Amino_acid)
    nx.write_pajek(G_labels, str(key) + '_labels_output1.net')
    os.chdir('..')


global kmers_nets
kmers_nets = {}
global nucoltides_kmers_nets
nucoltides_kmers_nets = {}
nucoltides_kmers_nets_edges_count = {}
replaced_nucoltides_kmers_nets = {}
replaced_nucoltides_kmers_nets_edges_count = {}
amino_acid_kmers_nets = {}

G = nx.Graph(name='CDR3 from DB')
nucoltides_kmers_nets_G = {}
amino_acids_kmers_nets_G = {}

NC_G = nx.Graph(name='AA Netwrok')
aa_count = {"TTT": [0], "CTT": [0], "ATT": [0], "GTT": [0],
            "TTC": [0], "CTC": [0], "ATC": [0], "GTC": [0],
            "TTA": [0], "CTA": [0], "ATA": [0], "GTA": [0],
            "TTG": [0], "CTG": [0], "ATG": [0], "GTG": [0],
            "TCT": [0], "CCT": [0], "ACT": [0], "GCT": [0],
            "TCC": [0], "CCC": [0], "ACC": [0], "GCC": [0],
            "TCA": [0], "CCA": [0], "ACA": [0], "GCA": [0],
            "TCG": [0], "CCG": [0], "ACG": [0], "GCG": [0],
            "TAT": [0], "CAT": [0], "AAT": [0], "GAT": [0],
            "TAC": [0], "CAC": [0], "AAC": [0], "GAC": [0],
            "TAA": [0], "CAA": [0], "AAA": [0], "GAA": [0],
            "TAG": [0], "CAG": [0], "AAG": [0], "GAG": [0],
            "TGT": [0], "CGT": [0], "AGT": [0], "GGT": [0],
            "TGC": [0], "CGC": [0], "AGC": [0], "GGC": [0],
            "TGA": [0], "CGA": [0], "AGA": [0], "GGA": [0],
            "TGG": [0], "CGG": [0], "AGG": [0], "GGG": [0],
            "---": [0]
            }
aa_most_frequent = {}
spike_dic = {}
amino_acid_list = ['ATA', 'ATC', 'ATT', 'ATG',
                   'ACA', 'ACC', 'ACG', 'ACT',
                   'AAC', 'AAT', 'AAA', 'AAG',
                   'AGC', 'AGT', 'AGA', 'AGG',
                   'CTA', 'CTC', 'CTG', 'CTT',
                   'CCA', 'CCC', 'CCG', 'CCT',
                   'CAC', 'CAT', 'CAA', 'CAG',
                   'CGA', 'CGC', 'CGG', 'CGT',
                   'GTA', 'GTC', 'GTG', 'GTT',
                   'GCA', 'GCC', 'GCG', 'GCT',
                   'GAC', 'GAT', 'GAA', 'GAG',
                   'GGA', 'GGC', 'GGG', 'GGT',
                   'TCA', 'TCC', 'TCG', 'TCT',
                   'TTC', 'TTT', 'TTA', 'TTG',
                   'TAC', 'TAT', 'TAA', 'TAG',
                   'TGC', 'TGT', 'TGA', 'TGG', ]

proteinNumber = 64
countrow = 0


