# -*- coding: utf-8 -*-
"""
Created on Sun Jan  2 17:53:17 2022
@author: saifr
"""
import os
import pandas as pd
from SBToolStart import startf
from SBToolunique import uniqueK
from SBToolVar import var
from SBFirstClustering import firstG
from SBToolSecondClustering import secondG
from SBToolPreFinal import linkdata

if __name__ == "__main__":
    def db(n):
        switcher = {
            1: "covid_vaccine_new",
            2: "covidpublished",
        }
        return switcher.get(n, "there is no such DB")


    def dbp(n):
        switcher = {
            1: "\nID options: 3  4  5  6  7",
            2: "\nID options: 3  4  5  6  7  8  9  10  14  15  16  17  18",
        }
        return switcher.get(n, "there is no such ID")


    while True:
        n = input("1 - covid_vaccine_new\n2 - covidpublished\nEnter number of DataBase you like to use: ")
        if int(n) not in (1, 2):
            print("there is no such DataBase number. choose again")
        else:
            break
    DB = db(int(n))
    print(dbp(int(n)))
    while True:
        if int(n) == 1:
            ID = input("Enter person/subjectID: ")
            if int(ID) not in (3, 4, 5, 6, 7):
                print("there is no such ID. choose again")
                print(dbp(int(n)))
            else:
                break
        elif int(n) == 2:
            ID = input("Enter person/subjectID:")
            if int(ID) not in (3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18):
                print("there is no such ID. choose again")
                print(dbp(int(n)))
            else:
                break
    #####################################

    # os.makedirs() method will raise
    # an OSError if the directory
    # to be created already exists
    # But It can be suppressed by
    # setting the value of a parameter
    # exist_ok as True
    # Directory
    directory = "{}_id{}_Results".format(DB, ID)
    try:
        os.makedirs(directory, exist_ok=True)
        print("Directory '%s' created successfully" % directory)
    except OSError as error:
        print("Directory '%s' can not be created" % directory)

    toFile = "{}\{}_id{}.csv".format(directory, DB, ID)
    toFile1 = "{}\{}_AA_id{}.csv".format(directory, DB, ID)
    toFile2 = "{}\{}_id{}_seqK.csv".format(directory, DB, ID)
    u_path = "{}\{}_id{}_uniqeKmers.csv".format(directory, DB, ID)
    remainV = "{}\{}_id{}_VarRemain.csv".format(directory, DB, ID)
    firstClusterPath = '{}\{}_id{}_clusters1.csv'.format(directory, DB, ID)
    secondClusteringPath = '{}\{}_id{}_clusters2.csv'.format(directory, DB, ID)
    lastStep = "{}\{}_id{}_final.csv".format(directory, DB, ID)

    ######################################################################

    startf(toFile, toFile1, toFile2, DB, ID)

    ######################################################################
    col_list = ["k-mer", "ai", "clone_id"]
    kmerss = pd.read_csv(toFile2, usecols=col_list)

    uniqueK(kmerss, u_path)
    ######################################################################
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



