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
