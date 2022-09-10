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
# A LIBRAR0Y2 |TO WORK VIEW PROGRESS IN CONSOLE WINDOW
#
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
        global remainV
        remainV = "{}\{}_id{}_VarRemain.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)
        global firstClusterPath
        firstClusterPath = '{}\{}_id{}_clusters1.csv'.format(str(os.getcwd() + "\\" + directory), DB, ID)
        global secondClusteringPath
        secondClusteringPath = '{}\{}_id{}_clusters2.csv'.format(str(os.getcwd() + "\\" + directory), DB, ID)
        global lastStep
        lastStep = "{}\{}_id{}_final.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)
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
class NetworkNode:
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
        cmd = "select seq.*, coll.*,sammet.* from {}.sequences as seq inner join {}.sequence_collapse as coll on seq.ai=coll.seq_ai inner join {}.sample_metadata as sammet on sammet.sample_id=seq.sample_id WHERE seq.functional=1 AND coll.instances_in_subject !=0 AND coll.copy_number_in_subject > 1  AND seq.deletions is null AND  seq.insertions is null AND LENGTH(seq.cdr3_aa)={} AND(sammet.value like'%spike+%' OR sammet.value like'%spike-%') ORDER BY RAND()".format(
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

                            a_node = NetworkNode(tripleAA, a_sequence, a_subject, a_seq_id, i)
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
                        a_node = NetworkNode(nuclotides_sub_seq, a_sequence, a_subject, a_seq_id, a_position)
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
            ['seq_id', 'sequence', 'TranslatedSeq', 'TranslatedGermline', 'ai', 'subject_id', 'clone_id', 'sample_id',
             'cdr3_aa', 'spike'])
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
                            spkie_count = [0, 0, 0]
                            spike_dic[tripleAA] = spkie_count
                        if "spike+" in str(spike):
                            x = spike_dic.get(tripleAA)
                            x[0] += 1
                            spike_dic[tripleAA] = x
                        elif "spike-" in str(spike):
                            x = spike_dic.get(tripleAA)
                            x[1] += 1
                            spike_dic[tripleAA] = x
                        else:
                            x = spike_dic.get(tripleAA)
                            x[2] += 1
                            spike_dic[tripleAA] = x

                        if protein[dna[i:i + 3]] != "-" and protein[dna[i:i + 3]] != "*":
                            a_sequence = line['sequence']  # protein a node
                            a_subject = line['subject_id']
                            a_seq_id = line['seq_id']
                            # protein a node

                            a_node = NetworkNode(tripleAA, a_sequence, a_subject, a_seq_id, i)
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
                        a_node = NetworkNode(nuclotides_sub_seq, a_sequence, a_subject, a_seq_id, a_position)
                        aa_count[dna[i:i + 3]][0] += 1
                        aa_count[dna[i:i + 3]].append(a_node)
                        old_list = nucoltides_kmers_nets.get(i, [])
                        old_list.append(a_node)
                        nucoltides_kmers_nets[i] = old_list
                        nuclotides_sub_seq = ""
                j = j + 3

                temp_pos = temp_pos + 9
            germProtein_sequence = ""
            for i in range(0, len(germ) - (len(germ) % 3), 3):
                if germ[i] == "N" and germ[i + 1] == "N" and germ[i + 2] == "N":
                    germProtein_sequence += "x"
                elif germ[i] == "N" or germ[i + 1] == "N" or germ[i + 2] == "N" and germ[i] != "-" and germ[
                    i + 1] != "-" and germ[i + 2] != "-":
                    germProtein_sequence += "x"
                elif germ[i] == "-" and germ[i + 1] == "-" and germ[i + 2] == "-":
                    germProtein_sequence += protein[germ[i:i + 3]]
                elif germ[i] == "-" or germ[i + 1] == "-" or germ[i + 2] == "-":
                    i = i + 0;
                else:
                    germProtein_sequence += protein[germ[i:i + 3]]
            csv_writer.writerow(
                [seqID, dna, protein_sequence, germProtein_sequence, ai, subjectID, cloneID, sampleID, cdr3_seq, spike])
            print("Translating DONE!")
            save_fasta.write(dna)
            count += 1
            save_fasta.write("\n>seq" + str(count) + "\n")
        save_fasta.close()

        # Matrix for the mutations
        def build_matrix(rows, cols):
            matrix = []
            for r in range(0, rows):
                matrix.append([0 for c in range(0, cols)])
            return matrix

        # Mutation function
        def mutatedFunc(seqAA, germAA):
            global flag
            flag = 0
            vec = build_matrix(2, len(seqAA))
            if len(seqAA) != len(germAA):
                csv_writer1.writerow([seqAA, germAA])
                flag = 1
            else:
                for i in range(0, len(seqAA), 1):
                    vec[0][i] = i + 1
                    # print(seqAA[i],germAA[i])
                    if seqAA[i] != germAA[i] and seqAA[i] != "x" and seqAA[i] != "-" and germAA[i] != "x" and germAA[
                        i] != "-" and seqAA[i] != "*" and germAA[i] != "*":
                        vec[1][i] = 1
            return vec

        print("AAsequence-mutations starts ....")
        with open(toFile, 'r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            with open(toFile1, 'w', newline='') as new_file, open(removed, 'w', newline='') as nfile:
                csv_writer = csv.writer(new_file)
                csv_writer.writerow(
                    ['ai', 'sequence', 'seq_id', 'translatedSeq', 'translatedGerm', 'vector', 'subject_id', 'clone_id',
                     'sample_id'])
                csv_writer1 = csv.writer(nfile)
                csv_writer1.writerow(['translatedSeq', 'translatedGerm'])
                for line in (csv_reader):
                    seq = (line['sequence'])
                    seqID = (line['seq_id'])
                    ai = (line['ai'])
                    seqAA = (line['TranslatedSeq'])
                    germAA = (line['TranslatedGermline'])
                    cloneID = (line['clone_id'])
                    sampleID = (line['sample_id'])
                    subjectID = (line['subject_id'])
                    vec1 = mutatedFunc(seqAA, germAA)
                    vector = []
                    if flag != 1:
                        for i in range(len(vec1[0])):
                            if vec1[1][i] == 1:
                                vector.append(i + 1)
                        csv_writer.writerow([ai, seq, seqID, seqAA, germAA, vector, subjectID, cloneID, sampleID])
            print("AAsequence-mutations DONE!")

        def kmersFunc(AA, k):
            global start
            start = 0
            p = 0
            kmer = ""
            x = 1
            s = 1
            while (x != 20):
                if AA[-s] != "-":
                    x += 1
                    s += 1
                else:
                    s += 1

            for i in range(0, len(AA), 1):
                if AA[i] == "x" or AA[i] == "-":
                    i += 0
                    start += 1
                else:
                    p1 = i
                    for q in range(i, (len(AA) - s) + 1, 1):
                        for j in range(q, len(AA), 1):
                            if AA[j] == "-" and kmer == "":
                                j = q + 1
                                p1 = j
                                break
                            if AA[j] == "-":
                                j += 0
                            else:
                                p += 1
                                kmer += AA[j]
                                if p == k:
                                    p = 0
                                    p2 = j
                                    pos = (p1 + 1, p2 + 1)
                                    # print(AA[p1:p2+1])
                                    # print(kmer)
                                    csv_writer1.writerow([kmer, pos, seqID, ai, subjectID, cloneID, sampleID])
                                    j = q + 1
                                    p1 = j
                                    kmer = ""
                                    break
                    break

        i = 0
        print("Kmers extraction starts ....")
        with open(toFile1, 'r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            with open(toFile2, 'w', newline='') as new_file1:
                csv_writer1 = csv.writer(new_file1)
                csv_writer1.writerow(['k-mer', 'position', 'seq_id', 'ai', 'subject_id', 'clone_id', 'sample_id'])
                for line in csv_reader:
                    KmerS = (line['translatedSeq'])
                    seqID = (line['seq_id'])
                    ai = (line['ai'])
                    cloneID = (line['clone_id'])
                    sampleID = (line['sample_id'])
                    subjectID = (line['subject_id'])
                    # Function for the k-mers!
                    kmersFunc(KmerS, 20)
                print("Kmers extraction DONE!")
        print("First step Ends\n")



def saif_and_bashara_code():
    col_list = ["k-mer", "ai", "clone_id"]
    kmerss = pd.read_csv(toFile2, usecols=col_list)

    uniqueK(kmerss, u_path)
 ################################################################
    col_list = ["kmer", "id"]
    col_list1 = ["kmer"]
    kmers0 = pd.read_csv(u_path, usecols=col_list)
    kmers1 = pd.read_csv(u_path, usecols=col_list1)
    var(kmers0, kmers1, remainV)
    ######################################################################
    kmers = pd.read_csv(remainV)
    print(len(kmers))
    firstG(kmers, firstClusterPath)

    ######################################################################
    # Second Step Clusters
    dff = pd.read_csv(firstClusterPath)
    secondG(dff, secondClusteringPath)
    ######################################################################
    # prefinal
    ref1 = pd.read_csv(toFile2)
    ref1.rename(columns={'k-mer': 'kmer'}, inplace=True)
    ref2 = pd.read_csv(firstClusterPath)
    df = pd.read_csv(secondClusteringPath)
    linkdata(df, ref1, ref2, lastStep)

#####################################################################################################################

tqdm.pandas()
from collections import Counter


# ID=10
# path=r"/home/saif/metadata/covid/covid{}_seqK.csv".format(ID)
# name="/home/saif/zak/cov1/cov{}_byScore.csv".format(ID)
# col_list = ["k-mer","ai","clone_id"]
# kmerss=pd.read_csv(path,usecols=col_list)

def uniqueK(kmers, name):
    print("Second step (extract the unique k-mers):\n")

    # print("unique starts here")
    kmers = kmers.rename(columns={'Unique-SeqKmer': 'Kmers'}, inplace=False)
    kmers = kmers.rename(columns={'k-mer': 'Kmers'}, inplace=False)

    print("number of k-mers before:" + str(len(kmers)))
    kmers['id'] = kmers.index

    l = kmers.values.tolist()

    d = {}
    for i in tqdm(l):
        d[i[0]] = []
    for j in tqdm(l):
        if (j[3] not in d[j[0]]):
            d[j[0]].append(j[3])

    new_df = pd.DataFrame(kmers['Kmers'])
    l = new_df.values.tolist()
    flat_list = []
    for sublist in l:
        for item in sublist:
            flat_list.append(item)
    l = flat_list

    def count_uniqe(lst):
        new_vals = Counter(l).most_common()
        new_vals = new_vals[::1]  # this sorts the list in scending order
        return new_vals

    new_list = count_uniqe(l)

    df = pd.DataFrame(new_list, columns=['kmer', 'score'])

    print("number of unique k-mers: " + str(len(df)))

    def inx(kmer):
        return d[kmer][0]

    df['id'] = df.kmer.apply(inx)
    df.to_csv(name, index=False)
    print("Second step Ends\n")

def get_spike_for_kmer(key):
    appearances_count = spike_dic.get(key)
    if float(appearances_count[0] + appearances_count[1]) * 0.9 <= float(appearances_count[0]):
        return int(1)
    elif float(appearances_count[0] + appearances_count[1]) * 0.9 <= float(appearances_count[1]):
        return int(-1)
    else:
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
        spike_value = get_spike_for_kmer(x)
        tmplst.append(spike_value)
        spike_plus_sum += get_spike_plus_sum(x)
        spike_minus_sum += get_spike_minus_sum(x)
        if spike_value == 1:
            appearances_count = spike_dic.get(x)
            if appearances_count[0] > 0:
                with open("spike+", 'a+', newline='') as new_file:
                    csv_writer = csv.writer(new_file)
                    csv_writer.writerow(
                        ["tri:" + str(x), "spike+:" + str(appearances_count[0]), "spike-:" + str(appearances_count[1])
                         ])

        elif spike_value == -1:
            appearances_count = spike_dic.get(x)
            if appearances_count[1] > 0:
                with open("spike_-.txt", 'a+', newline='') as new_file:
                    csv_writer = csv.writer(new_file)
                    csv_writer.writerow(
                        ["tri:" + str(x), "spike+:" + str(appearances_count[0]), "spike-:" + str(appearances_count[1])
                         ])
        else:
            appearances_count = spike_dic.get(x)
            with open("spike_neutral.txt", 'a+', newline='') as new_file:
                csv_writer = csv.writer(new_file)
                csv_writer.writerow(
                    ["tri:" + str(x), "spike+:" + str(appearances_count[0]), "spike-:" + str(appearances_count[1])
                     ])
    with open("all_spikes.txt", 'a+', newline='') as new_file:
        csv_writer = csv.writer(new_file)
        for row in spike_dic.keys():
            csv_writer.writerow(
                ["tri:" + row,
                 "spike+:" + str(spike_dic.get(row)[0]),
                 "spike-:" + str(spike_dic.get(row)[1]),
                 "spike_N:" + str(spike_dic.get(row)[2])
                 ])

    data["Y"] = tmplst
    font1 = {'family': 'serif', 'color': 'blue', 'size': 20}
    font2 = {'family': 'serif', 'color': 'darkred', 'size': 15}
    plt.title("number of spike plus =" + str(spike_plus_sum) + " number of spike minus =" + str(spike_minus_sum),
              fontdict=font1)
    plt.scatter(x="X", y="Y", marker='o', data=data)
    plt.xlabel("AA sub sequence 3 mers", fontdict=font2)
    plt.ylabel("'spike+'=1; 'spike-'=0", fontdict=font2)
    plt.savefig("_spike_fig.png")
    # plt.show(block=False)


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
    plt.plot(list(nuclotide_count_dict.keys()), [nuclotide_count_dict.get(x)[0] for x in nuclotide_count_dict.keys()])
    plt.savefig(str(output_filename + "tempfig.png"))
    # plt.show(block=False)

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




def diff_letters(a, b):
    cnt = 0
    for i in range(len(a)):
        if a[i] != b[i] and a[i] != 'x' and b[i] != 'x':
            cnt += 1
            if cnt > 1:
                break
    return cnt


def fiten(row):
    row = row[:10]
    return row


def laten(row):
    row = row[10:]
    return row


def wrap_f(x):
    return fiten(x['kmer'])


def wrap_l(x):
    return laten(x['kmer'])


def var(kmers0, kmers1, name):
    print("Third step (1-mismatch combined k-mers):\n")
    kmers = kmers1.copy()
    kmers1['tmpid'] = kmers1.index
    kmers1['var'] = 0

    print("number of k-mers before combination: " + str(len(kmers0)))

    kmers1['id'] = kmers0['id']

    selected_columns = kmers1[["tmpid", "var"]]
    new_df = selected_columns.copy()

    all_list = new_df.values.tolist()
    All = {}
    st = {}
    for l in tqdm(all_list):
        All[l[0]] = l[1]
        st[l[0]] = l[1]
    print(kmers)
    kmers['first'] = kmers.progress_apply((lambda r: fiten(r['kmer'])), axis=1)
    kmers['last'] = kmers.progress_apply((lambda r: laten(r['kmer'])), axis=1)
    kmers['kmer'] = kmers.index

    lastlist = kmers.values.tolist()
    removed = {}
    for i in tqdm(lastlist):
        removed[i[0]] = []
    last = {}
    for l in tqdm(lastlist):
        last[l[0]] = l[2]
    first = {}
    for l in tqdm(lastlist):
        first[l[0]] = l[1]
    firstd = {}
    for i in tqdm(lastlist):
        firstd[i[1]] = []
    for j in tqdm(lastlist):
        firstd[j[1]].append(j[0])
    lastd = {}
    for i in tqdm(lastlist):
        lastd[i[2]] = []
    for j in tqdm(lastlist):
        lastd[j[2]].append(j[0])

    def mapf():
        r = []
        for l in tqdm(firstd.values()):
            if len(l) > 1:
                for i in range(len(l)):
                    if st[l[i]] != 0:
                        continue
                    for j in range(i + 1, len(l)):
                        if st[l[j]] == 0:
                            if diff_letters(last[l[i]], last[l[j]]) <= 1:
                                All[l[i]] += 1
                                st[l[j]] += 1
                                removed[l[i]].append(l[j])
                                r.append(l[j])

        # second comparision
        for e in tqdm(lastd.values()):
            if len(e) > 1:
                for x in range(len(e)):
                    if st[e[x]] != 0:
                        continue
                    for y in range(x + 1, len(e)):
                        if st[e[y]] == 0:
                            if diff_letters(first[e[x]], first[e[y]]) <= 1:
                                All[e[x]] += 1
                                st[e[y]] += 1
                                removed[e[x]].append(e[y])
                                r.append(e[y])
        return r

    trash = mapf()
    ss = list(set(trash))
    print("\nnumber of kmers to be deleted: {0:d}".format(len(ss)))

    var = pd.DataFrame.from_dict(All, orient='index', columns=['variance'])
    final = kmers1.drop(columns=['tmpid'])
    final['var'] = var['variance']
    final = final.sort_values(by=['var'], ascending=False)
    final = final.drop(trash)
    print("Third step Ends\n")
    final.to_csv(name, index=False)
    print("\nthe final number of kmers: {0:d}".format(len(final)))



def firstG(kmers, name):
    print("first clustering starts...")
    selected_columns = kmers[["kmer"]]
    df = selected_columns.copy()

    ## functions
    def fitri(row):
        row = row[:3]
        return row

    def latri(row):
        row = row[17:]
        return row

    def betw(row):
        row = row[3:17]
        return row

    def wrap_f(x):
        return fitri(x['kmer'])

    def wrap_l(x):
        return latri(x['kmer'])

    def wrap_b(x):
        return betw(x['kmer'])

    def diff_letters(a, b):
        cnt = 0
        for i in range(len(a)):
            if a[i] == b[i] and a[i] != 'x' and b[i] != 'x':
                cnt += 1
                if cnt > 4:
                    return 1
        return 0

    df['first'] = df.progress_apply(wrap_f, axis=1)
    df['last'] = df.progress_apply(wrap_l, axis=1)
    df['k'] = df.index
    df['kmers'] = df.progress_apply(wrap_b, axis=1)
    df['id'] = kmers['id'].copy()

    lastlist = df.values.tolist()
    ids = {}
    for l in tqdm(lastlist):
        ids[l[3]] = l[5]

    st = {}
    for l in tqdm(lastlist):
        st[l[3]] = 0

    km = {}
    for l in tqdm(lastlist):
        km[l[3]] = l[0]

    k = {}
    for l in tqdm(lastlist):
        k[l[3]] = l[4]

    last = {}
    for l in tqdm(lastlist):
        last[l[3]] = l[2]

    first = {}
    for l in tqdm(lastlist):
        first[l[3]] = l[1]

    firstd = {}
    for i in tqdm(lastlist):
        firstd[i[1]] = []

    for j in tqdm(lastlist):
        firstd[j[1]].append(j[3])

    lastd = {}
    for i in tqdm(lastlist):
        lastd[i[2]] = []

    for j in tqdm(lastlist):
        lastd[j[2]].append(j[3])

    cl = -1
    data = []
    for l in tqdm(firstd.values()):
        for i in range(len(l)):
            if st[l[i]] != 0:
                continue
            cl += 1
            data.append([km[l[i]], cl, ids[l[i]]])
            for j in range(i + 1, len(l)):
                if last[l[i]] == last[l[j]]:
                    if diff_letters(k[l[i]], k[l[j]]) == 1:
                        st[l[j]] += 1
                        data.append([km[l[j]], cl, ids[l[j]]])
                    else:
                        st[l[j]] += 1
                        data.append([km[l[j]], -1, ids[l[j]]])

    dff = pd.DataFrame(data, columns=['kmer', 'ClusterID', 'id'])

    dff = dff[dff['ClusterID'] != -1]
    dff = dff.reset_index(drop=True)

    print("number of kmers after first clusterng: " + str(len(dff)))
    print(dff.head())
    dff.to_csv(name, index=False)




def secondG(dff, name):
    print("second clustering starts...")
    maxx = dff.ClusterID.max()

    def prkmer(list1, j):
        AA = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0,
              'F': 0, 'P': 0, 'S': 0, 's': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}
        new = ''
        for i in range(20):
            for k in list1:
                if k[i] != 'x':
                    AA[k[i]] += 1
            new += max(AA.items(), key=operator.itemgetter(1))[0]
            AA = dict.fromkeys(AA, 0)
        new += ':'
        new += str(j)
        return new

    lol = [dff[dff['ClusterID'] == i].kmer.tolist() for i in tqdm(range(maxx + 1))]
    newl = [[prkmer(lol[i], i), i] for i in tqdm(range(len(lol)))]

    newdf = pd.DataFrame(newl, columns=['kmer', 'id'])
    klis = newdf.kmer.tolist()
    klis = list(set(klis))
    st1 = {}
    for l in tqdm(klis):
        st1[l] = 0

    def hamming(a, b):
        cnt = 0
        for i in range(20):
            if a[i] == b[i] and a[i] != 'x' and b[i] != 'x':
                cnt += 1
                if cnt > 10:
                    return 1
        return 0

    cl = -1
    data = []
    for i in tqdm(range(len(klis))):
        if st1[klis[i]] != 0:
            continue
        cl += 1
        data.append([klis[i], cl])
        for j in range(i + 1, len(klis)):
            if st1[klis[j]] == 0:
                if hamming(klis[i], klis[j]) == 1:
                    st1[klis[j]] += 1
                    data.append([klis[j], cl])

    dat = pd.DataFrame(data, columns=['kmer', 'clusterid'])
    print(dat.head())

    dat.to_csv(name, index=False)




def indx(row):
    row = row[21:]
    return row


def wrap_indx(x):
    return indx(x['kmer'])


def linkdata(df, ref1, ref2, name):
    print("start linking data...")
    selected_columns = df[["kmer"]]
    df1 = selected_columns.copy()
    df['kmerid'] = df1.progress_apply(wrap_indx, axis=1)

    maxx = df.clusterid.max()

    bigc = []
    for i in tqdm(range(maxx + 1)):
        newl = []
        flat_list = []
        d = df[df['clusterid'] == i]
        blist = d.kmerid.tolist()
        for i in blist:
            a = ref2.loc[ref2['ClusterID'] == int(i)]
            a = a.id.tolist()
            newl.append(a)
        flat_list = [item for sublist in newl for item in sublist]
        bigc.append(flat_list)

    big = pd.DataFrame()
    for i in tqdm(range(len(bigc))):
        aa = ref1.iloc[bigc[i]]
        l = aa.kmer.tolist()
        aa.insert(7, 'cluster', str(i))
        aa = aa.reset_index(drop=True)
        big = pd.concat([big, aa.copy()])

    x = big.cluster.value_counts()
    x = x.sort_index()
    p = pd.DataFrame(x)
    p = p.reset_index()
    p = p.rename(columns={'index': 'clnum', 'cluster': 'count'}, inplace=False)
    lastlist = p.values.tolist()
    num = {}
    for l in tqdm(lastlist):
        num[l[0]] = l[1]

    def nums(x):
        return num[x]

    def pos(x):
        x = x[1:4]
        if x[2] == ',':
            return int(x[:2])
        else:
            return int(x)

    big['cnt'] = big.apply(lambda row: nums(row['cluster']), axis=1)
    big['pos'] = big.apply(lambda row: pos(row['position']), axis=1)
    print("number of final kmers:" + str(len(big)))
    print(big.head())
    big.to_csv(name, index=False)





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
