import os
import datetime
import pickle
#####################################################################################################################
# @DESCRIPTION : DECLARING AND CREATOING THE DIORCTORES NEDDED
# @INPUT : DB : DATABASE NAME , ID: SUBJECT ID
# @OUTPUT : ---
#####################################################################################################################
global toFile
global toFile1
global toFile2
global u_path
global Amino_acid
global nucleotides
global pajek
global aa_count
global nucoltides_kmers_nets
global directory
def create_files(DB, ID, cdr_length, selected_analyzed_position, startposition, endposition, rows):
    # CREATING THE DIRECTORY

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

        toFile = "{}\{}emp_id{}.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)

        toFile1 = "{}\{}emp_AA_id{}.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)

        toFile2 = "{}\{}emp_id{}_seqK.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)

        u_path = "{}\{}_id{}_uniqeKmers.csv".format(str(os.getcwd() + "\\" + directory), DB, ID)

        Amino_acid = "{}\Amino_acid".format(str(os.getcwd() + "\\" + directory))
        isdir = os.path.isdir(Amino_acid)
        if not isdir:
            os.makedirs(Amino_acid)

        nucleotides = "{}\\nucleotides".format(str(os.getcwd() + "\\" + directory))
        isdir = os.path.isdir(nucleotides)
        if not isdir:
            os.makedirs(nucleotides)

        pajek = "{}\pajek".format(str(os.getcwd() + "\\" + directory))
        isdir = os.path.isdir(pajek)
        if not isdir:
            os.makedirs(pajek)
    except OSError as error:
        print("Directory '%s' can not be created" % directory)

def export_preproccessed_data():
    pickle.dump(nucoltides_kmers_nets, open('nucoltides_kmers_nets.pkl', 'wb'))
    pickle.dump(aa_count, open('protein_count.pkl', 'wb'))


def readTrueData(name):
    fName = str(name + '.pkl')
    with open(fName, 'rb') as f:  # with statement avoids file leak
        # Match load with file object, and provide encoding for Py2 str
        return pickle.load(f)

def import_preproccessed_data():

    nucoltides_kmers_nets = readTrueData('nucoltides_kmers_nets')

    aa_count = readTrueData('protein_count')


#TODO(): add write html,.net,.fasta files to use in the creation of the files