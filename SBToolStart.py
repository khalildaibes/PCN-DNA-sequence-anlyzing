import mysql.connector
import csv
from tqdm import tqdm


def startf(toFile, toFile1, toFile2, DB, ID):
    print("\nFirst step (translate sequences to amino acids):")
    print("First step strats...\n")

    removed = "trash.csv"
    cmd = "select seq.*, coll.* from {}.sequences as seq inner join {}.sequence_collapse as coll on seq.ai=coll.seq_ai WHERE seq.subject_id={} AND seq.functional=1 AND coll.instances_in_subject !=0 AND coll.copy_number_in_subject > 1  AND seq.deletions is null AND  seq.insertions is null".format(
        DB, DB, ID)
    # DATABASE connection
    mydb = mysql.connector.connect(
        host="132.75.249.34",
        user="guest",
        passwd="Guest2020!",
    )
    ########################
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
    # Command for getting the sequences translated:
    command = (cmd)
    mycursor.execute(command)
    seq = mycursor.fetchmany(1000)
    print("Translating starts ....")
    with open(toFile, 'w', newline='') as new_file:
        csv_writer = csv.writer(new_file)
        csv_writer.writerow(
            ['seq_id', 'sequence', 'TranslatedSeq', 'TranslatedGermline', 'ai', 'subject_id', 'clone_id', 'sample_id'])
        for line in tqdm(seq):
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
        print("Translating DONE!")

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







