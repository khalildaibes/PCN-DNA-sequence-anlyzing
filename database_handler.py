# -*- coding: utf-8 -*-
import mysql.connector


def create_connection():
    mydb = mysql.connector.connect(
        host="132.75.249.34",
        user="guest",
        passwd="Guest2020!",
    )
    ########################
    return mydb

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

def connect_to_db():
    mydb = create_connection()
    mycursor = mydb.cursor(dictionary=True)
    return mycursor