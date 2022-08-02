# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 14:36:56 2021

@author: saifr
"""

import csv
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()


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
                            csv_writer1.writerow([kmer, pos, seq, seqID, ai, subjectID, cloneID, sampleID])
                            j = q + 1
                            p1 = j
                            kmer = ""
                            break
            break


def kmersFunc1(AA, k):
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
    p1 = start
    for q in range(start, (len(AA) - s) + 1, 1):
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
                    csv_writer2.writerow([kmer, pos, seq, seqID, ai, subjectID, cloneID, sampleID])
                    j = q + 1
                    p1 = j
                    kmer = ""
                    break


"------------------------------------------------------"
# change your files here
# Take CSV file and get the both sequence and germline kmers and save each one in CSV

fromFile = filedialog.askopenfilename()
toFile1 = filedialog.asksaveasfilename(title="save file KmersFromSequence",
                                       filetypes=(("CSV files", "*.csv"), ("all files", "*.*")), initialfile=".csv")
toFile2 = filedialog.asksaveasfilename(title="save file KmersFromGermline",
                                       filetypes=(("CSV files", "*.csv"), ("all files", "*.*")), initialfile=".csv")
"------------------------------------------------------"
i = 0
with open(fromFile, 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    with open(toFile1, 'w', newline='') as new_file1, open(toFile2, 'w', newline='') as new_file2:
        csv_writer1 = csv.writer(new_file1)
        csv_writer2 = csv.writer(new_file2)
        csv_writer1.writerow(['k-mer', 'position', 'seq', 'seq_id', 'ai', 'subject_id', 'clone_id', 'sample_id'])
        csv_writer2.writerow(['k-mer', 'position', 'seq', 'seq_id', 'ai', 'subject_id', 'clone_id', 'sample_id'])
        for line in csv_reader:
            KmerS = (line['translatedSeq'])
            KmerG = (line['translatedGerm'])
            seq = (line['sequence'])
            seqID = (line['seq_id'])
            ai = (line['ai'])
            cloneID = (line['clone_id'])
            sampleID = (line['sample_id'])
            subjectID = (line['subject_id'])

            # Function for the k-mers!
            kmersFunc(KmerS, 3)
            kmersFunc1(KmerG, 3)
            break
        print("Done!")


