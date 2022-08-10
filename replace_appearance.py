from files_handler import toFile
from  prerequest import *
import database_handler as db_connector
from csv import DictReader
import numbers
from files_handler import *
import tqdm
from netwroks_classes import *
from threading import Thread
from distances_and_comparing_functions import *
from netwrok_to_files import *
#####################################################################################################################
# @DESCRIPTION : replace each AA with its most frequented nucleotide way
# @INPUT : DB: DATABASE NAME, ID, cdr_length, selected_analyzed_position, startposition, endposition, rows
# @OUTPUT : A list that holds dicts as rows
#####################################################################################################################
def replace_nucleotide(DB, ID, selected_analyzed_position, startposition, endposition, rows, cdr_length=14):
    print("\nFirst step (translate sequences to amino acids):")
    print("First step strats...\n")
    # CREATING A CONNECTION WITH THE DATABASE
    mycursor=db_connector.connect_to_db()
    create_files(DB, ID, cdr_length, selected_analyzed_position, startposition, endposition, rows)
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

    command = db_connector.select_query(DB, ID, cdr_length)
    mycursor.execute(command)
    global Selected_Number_Of_Rows
    kmers_G = nx.Graph(name='kmers_of_position')
    seq = db_connector.fetch_based_on_row_count(rows, mycursor)
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
