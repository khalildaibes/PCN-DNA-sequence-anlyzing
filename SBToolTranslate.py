# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 14:51:44 2021

@author: saifr
"""
import mysql.connector
import csv
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

# This part is for saving the file loacation as csv file.

toFile = filedialog.asksaveasfilename(title="save file", filetypes=(("CSV files", "*.csv"), ("all files", "*.*")),
                                      initialfile=".csv")
"------------------------------------------------------"

# DATABASE connection

mydb = mysql.connector.connect(
    host="132.75.249.34",
    user="guest",
    passwd="Guest2020!",
)
"------------------------------------------------------"
mycursor = mydb.cursor(dictionary=True);

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

# 5 times

# stanard    ===> variblev (avg relatdness)
# i+1        ===> variblev (avg relatdness)
# i+2        ===> variblev (avg relatdness)
# i-1        ===> variblev (avg relatdness)
# i-2        ===> variblev (avg relatdness)

# ----------------------------------


# Command for getting the sequences translated:
# jhasdoinasldknasdjsaopdjsihdasoihndoasijdaop


command = (
    "select seq.*, coll.* from covidpublished.sequences as seq inner join covidpublished.sequence_collapse as coll on seq.ai=coll.seq_ai WHERE seq.subject_id=3 AND seq.functional=0 AND coll.instances_in_subject !=0 AND coll.copy_number_in_subject > 1  AND seq.deletions is null AND  seq.insertions is null"
)
mycursor.execute(command)
# [i:i+3]--->[i:dist-1]
seq = mycursor.fetchall()
with open(toFile, 'w', newline='') as new_file:
    csv_writer = csv.writer(new_file)
    csv_writer.writerow(
        ['seq_id', 'sequence', 'TranslatedSeq', 'TranslatedGermline', 'ai', 'subject_id', 'clone_id', 'sample_id'])
    for line in seq:
        dna = ""
        germ = ""
        protein_sequence = ""
        # fix the germline to match the cdr with N's
        germ = line['germline']
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
        for i in range(0, len(dna) - (len(dna) % 3), 3):
            if dna[i] == "N" and dna[i + 1] == "N" and dna[i + 2] == "N":
                protein_sequence += "x"
            elif dna[i] == "N" or dna[i + 1] == "N" or dna[i + 2] == "N" and dna[i] != "-" and dna[i + 1] != "-" and \
                    dna[i + 2] != "-":
                protein_sequence += "x"
            elif dna[i] == "-" and dna[i + 1] == "-" and dna[i + 2] == "-":
                protein_sequence += protein[dna[i:i + 3]]
            elif dna[i] == "-" or dna[i + 1] == "-" or dna[i + 2] == "-":
                i = i + 0;
            else:
                protein_sequence += protein[dna[i:i + 3]]
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
        csv_writer.writerow([seqID, dna, protein_sequence, germProtein_sequence, ai, subjectID, cloneID, sampleID])
    print("Done!")
























