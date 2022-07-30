# -*- coding: utf-8 -*-

import tkinter as tk
from tkinter import filedialog
import pandas as pd

root = tk.Tk()
root.withdraw()

# This part is for selecting the FASTA file you want to translate.
myfile = filedialog.askopenfilename()

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


## reading the file and parsing the data for SEQ % METADATA
def readFasta(file_name):
    sequences = []
    meta = []
    with open(file_name) as f:
        while True:
            met = f.readline().rstrip()
            seq = f.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            meta.append(met)
    return sequences, meta


## translate the sequences
def translateSeq(seq):
    seqlist = []
    for i, row in seq.iterrows():
        dna = ""
        protein_sequence = ""
        dna = row['sequence']
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
        seqlist.append(protein_sequence)
    return seqlist


## the main function
if __name__ == "__main__":
    seqs, mtas = readFasta(myfile)
    sequences = pd.DataFrame(seqs, columns=['sequence'])
    data = translateSeq(sequences)
    ## here we divide the cdr3 length by 3 to match the amino-acid length
    for i in range(len(mtas)):
        r = mtas[i]
        num = r[len(r) - 2:]
        newnum = int(int(num) / 3)
        mtas[i] = mtas[i].replace(num, str(newnum))

    ## writing the new file
    ofile = open("my_fasta.fasta", "w")
    for i in range(len(data)):
        ofile.write(">" + mtas[i] + "\n" + data[i] + "\n")
    ofile.close()
    print("Done!")
























