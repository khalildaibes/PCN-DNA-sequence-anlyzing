from files_handler import *
from database_handler import *
from prerequest import *
import tqdm
import csv
from netwroks_classes import *
def startf(DB, ID, selected_analyzed_position, startposition, endposition, rows, cdr_length=14):
    create_files(DB, ID, cdr_length, selected_analyzed_position, startposition, endposition, rows)
    print("\nFirst step (translate sequences to amino acids):")
    print("First step strats...\n")
    sequence_legth = 231
    removed = "trash.csv"


    mycursor = connect_to_db()

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

