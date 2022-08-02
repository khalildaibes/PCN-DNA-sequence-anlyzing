import csv
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

"------------------------------------------------------"
# browse the source file (FULL AMINO-ACID SEQUENCE)
fromFile = filedialog.askopenfilename()
"------------------------------------------------------"
# Where to save your output file (FULL MUTATIONS FILE)
toFile = filedialog.asksaveasfilename(title="save file", filetypes=(("CSV files", "*.csv"), ("all files", "*.*")),
                                      initialfile=".csv")
"------------------------------------------------------"


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


with open(fromFile, 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    with open(toFile, 'w', newline='') as new_file, open('E:/project/lp15/diffLenSeq.csv', 'w', newline='') as nfile:
        csv_writer = csv.writer(new_file)
        csv_writer.writerow(
            ['ai', 'sequence', 'seq_id', 'translatedSeq', 'translatedGerm', 'vector', 'subject_id', 'clone_id',
             'sample_id'])
        csv_writer1 = csv.writer(nfile)
        csv_writer1.writerow(['translatedSeq', 'translatedGerm'])
        for line in csv_reader:
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
    print("Done!")
