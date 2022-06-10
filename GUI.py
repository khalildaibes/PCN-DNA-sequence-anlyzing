from tkinter import *
#  THIS IS AN IMPORT TO IMPORT THE MAIN MODULE WHICH IS IMPORTANT TO HANDLE THE GUI
import main
# DEFINING GLOBAL VARIABLES THAT WILL HOLD THE DATA FOR THE ANALYZING PROGRESS

global fetched_flag
fetched_flag = False
from tkinter import messagebox

#####################################################################################################################
# THE CLASS DIGRAM THIS WILL HOLD THE PROGRAM PREPROGRESS METADATA
# SUCH AS THE DATABASE AND THE APPROACH AND ALSO THE
# INPUT: THE CONSTRICTOR RECIVES AS AN INPUT THE CLASS META DATA
#####################################################################################################################
class user:
    def __init__(self, username=None, password=None, subject=None, database=None, cdr_length=None,
                 selected_analyzed_position=0, startposition=0, endposition=0, rows='all'):
        self.username = username
        self.password = password
        self.subject = subject
        self.database = database
        self.cdr_length = cdr_length
        self.selected_analyzed_position = selected_analyzed_position
        self.startposition = startposition
        self.endposition = endposition
        self.rows = rows


global loading

#####################################################################################################################
# @DESCRIPTION : SETTING THE ROW_COUNT VARIABLE WHICH HOLD HOW MANY ROWS WE WANT TO FETCH
# @INPUT : HOW MANY ROWS WE WANT TO FETCH
# @OUTPUT : ---
#####################################################################################################################
def set_row_count(rows):
    user.rows = rows

#####################################################################################################################
# @DESCRIPTION : START FETCHING THE DATA AND CONSTRUCTING HE NETWORKS OBJECTS
# @INPUT : ---- THE DATA IS STORED IN THE USER OBJECT
# @OUTPUT : ---
#####################################################################################################################
def start_fetching_data():
    # main.select_query(user.database, user.subject,user.cdr_length)
    main.startf(user.database, user.subject, user.cdr_length, user.selected_analyzed_position, user.startposition,
                user.endposition, user.rows)
    fetched_flag = True

#####################################################################################################################
# @DESCRIPTION : THE REGISTER WINDOWS
# @INPUT : ---- MAIN SCREEN IN THE WINDOWS FOR THE STARTING OF THE APPLICATION
# @OUTPUT : ---
#####################################################################################################################
def register(main_screen):
    # The Toplevel widget work pretty much like Frame,
    # but it is displayed in a separate, top-level window.
    # Such windows usually have title bars, borders, and other “window decorations”.
    # And in argument we have to pass global screen variable

    register_screen = Toplevel(main_screen)
    register_screen.title("Register")
    register_screen.geometry("300x250")

    # Set text variables
    username = StringVar()
    password = StringVar()

    # Set label for user's instruction
    Label(register_screen, text="Please enter details below", bg="blue").pack()
    Label(register_screen, text="").pack()

    # Set username label
    username_lable = Label(register_screen, text="Username * ")
    username_lable.pack()

    # Set username entry
    # The Entry widget is a standard Tkinter widget used to enter or display a single line of text.

    username_entry = Entry(register_screen, textvariable=username)
    username_entry.pack()

    # Set password label
    password_lable = Label(register_screen, text="Password * ")
    password_lable.pack()

    # Set password entry
    password_entry = Entry(register_screen, textvariable=password, show='*')
    password_entry.pack()

    Label(register_screen, text="").pack()

    # Set register button
    Button(register_screen, text="Register", width=10, height=1, bg="blue",
           command=lambda: register_user(username_entry, password_entry)).pack()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE VALUE TO THE START POSITION
# @INPUT : THE VALUE REQUESTED TO BE AT THE START POSITION
# @OUTPUT : ---
#####################################################################################################################
def set_startposition(startposition, variable=None):
    # SET THE CDR LENGTH AS IF WE DECLARE START AND END THATS NOT CDR OR FW
    if str(startposition).isnumeric():
        user.startposition = startposition
        set_cdr_length(0)
    else:
        if str(startposition) != "":
            messagebox.showerror("error", "try again add only numbers")
            variable.set('0')

    print(user.startposition)
#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE VALUE TO THE END POSITION
# @INPUT : THE VALUE REQUESTED TO BE AT THE END POSITION
# @OUTPUT : ---
#####################################################################################################################
def set_endposition(endposition, variable=None):
    if str(endposition).isnumeric():
        user.endposition = endposition
    # SET THE CDR LENGTH AS IF WE DECLARE START AND END THATS NOT CDR OR FW
        set_cdr_length(0)
    else:
        if str(endposition)!="" and variable:
            messagebox.showerror("error", "try again add only numbers")
            variable.set('0')
    print(user.endposition)

#####################################################################################################################
# @DESCRIPTION : ASSIGNING VALUE TO THE SUBJECT/PERSON
# @INPUT : THE VALUE REQUESTED TO BE AT THE SUBJECT VALUE
# @OUTPUT : ---
#####################################################################################################################
def set_subject(subject):
    user.subject = subject
    print(user.subject)

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA  THE ANALYZING APPROACH
# @INPUT : THE VALUE REQUESTED TO BE AT THE SUBJECT VALUE
# @OUTPUT : ---
#####################################################################################################################
def set_selected_analyzed_position(selected_analyzed_position):
    user.selected_analyzed_position = selected_analyzed_position
    print(selected_analyzed_position)


#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE START AND END POSITIONS TO BE THE REQUESTED VALUES
# @INPUT : THE VALUE REQUESTED TO BE AT THE START AND END POSITIONS
# @OUTPUT : ---
#####################################################################################################################
def set_position(start, end):
    # SET THE CDR LENGTH AS IF WE DECLARE START AND END THATS NOT CDR OR FW
    set_cdr_length(0)
    user.startposition = start
    user.endposition = end
    print(user.startposition, user.endposition)

#####################################################################################################################
# @DESCRIPTION : SET THE REQUESTED CDR LENGTH
# @INPUT : THE VALUE REQUESTED TO BE AT CDR LENGTH
# @OUTPUT : ---
#####################################################################################################################
def set_cdr_length(cdr_length):
    user.cdr_length = cdr_length
    print(user.cdr_length)
#####################################################################################################################
#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED DATABASE
# @INPUT : THE VALUE REQUESTED TO BE AT DATABSE
# @OUTPUT : ---
#####################################################################################################################
def set_databse(database):
    user.database = database
    print(user.database)

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA A REQUEST TO ADD USERS
# @INPUT : THE VALUES REQUESTED TO BE AT USERS PASSWORD AND USERNAME
# @OUTPUT : ---
#####################################################################################################################
def register_user(username_entry, password_entry):
    # get username and password
    username_info = username_entry.get()
    password_info = password_entry.get()

    # Open file in write mode
    file = open(username_info, "w")

    # write username and password information into file
    file.write(username_info + "\n")
    file.write(password_info)
    file.close()

    username_entry.delete(0, END)
    password_entry.delete(0, END)

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def choose_subject(choosing_DB_screen, db):
    # DECLARE TOP LEVEL WINDOW
    choose_subject_screen = Toplevel(choosing_DB_screen)
    # choosing_DB_screen.withdraw()
    # SET TITLE
    choose_subject_screen.title("choose subject")
    # SETTINGT HE SCREEN HIGHT AND WIDTH
    choose_subject_screen.geometry("300x250")
    # DECLARING VARIABLE TO HOLD THE VALUE
    variable = StringVar(choose_subject_screen)
    # SETTING A TRACE CALLBACK TO THE FUNCTION CLICK EVENT
    variable.trace("w", lambda *args: set_subject(variable.get()))
    # SEETING THE DEFULT VALUE TO THE OPTION MENUE
    variable.set("3")  # default value
    # CALLING THE SET OIN TEH DEFULT VALUE
    set_subject(variable.get())
    # DECLARING A LABEL
    Label(choose_subject_screen, text="Choose subject ", bg="blue", width="300", height="2",
          font=("Calibri", 13)).pack()
    # VEIW OOPTIONS BASED ON SELECTED DB
    if (db == "covid_vaccine_new"):
        w = OptionMenu(choose_subject_screen, variable, "3", "4", "5", "6", "7", "8", "9", "all")
        w.pack()
        print("the db is", db)
    else:
        w = OptionMenu(choose_subject_screen, variable, "7", "8", "9", "all")
        w.pack()
        print("the db is", db)
    # DECLARE THE BUTTON TO PRESS
    Button(choose_subject_screen, text="Next", height="2", width="30",
           command=lambda: select_position(choose_subject_screen)).pack()
#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED CDR LENGTH WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def choose_cdr_length(choosing_DB_screen, db):
    # DECLARE TOP LEVEL WINDOW
    cdr_length_screen = Toplevel(choosing_DB_screen)
    # SET TITLE
    cdr_length_screen.title("choose CDR length")
    # SETTINGT HE SCREEN HIGHT AND WIDTH
    cdr_length_screen.geometry("300x250")
    # DECLARING VARIABLE TO HOLD THE VALUE
    variable = StringVar(cdr_length_screen)
    # SETTING A TRACE CALLBACK TO THE FUNCTION CLICK EVENT
    variable.trace("w", lambda *args: set_cdr_length(variable.get()))
    # SEETING THE DEFULT VALUE TO THE OPTION MENUE
    variable.set("288")  # default value
    # CALLING THE SET OIN TEH DEFULT VALUE
    set_cdr_length(variable.get())
    # DECLARING A LABEL
    Label(cdr_length_screen, text="Choose CDR Length ", bg="blue", width="300", height="2", font=("Calibri", 13)).pack()
    # DECLKARE OPTION MENUE OPTIONS
    w = OptionMenu(cdr_length_screen, variable, "14", "20", "24", "30", "not specific")
    w.pack()
    # DECLARE THE BUTTON TO PRESS
    Button(cdr_length_screen, text="Next", height="2", width="30",
           command=lambda: choose_subject(choosing_DB_screen, db, variable.get())).pack()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def choose_database(main_screen):
    # DECLARE TOP LEVEL WINDOW
    choosing_DB_screen = Toplevel(main_screen)
    # main_screen.withdraw()
    # SET TITLE
    choosing_DB_screen.title("choose database")
    # SETTINGT HE SCREEN HIGHT AND WIDTH
    choosing_DB_screen.geometry("300x250")
    # DECLARING VARIABLE TO HOLD THE VALUE
    variable = StringVar(choosing_DB_screen)
    # DECLARING VARIABLE TO HOLD THE VALUE
    variable.set("covid_vaccine_new")  # default value
    # SETTING A TRACE CALLBACK TO THE FUNCTION CLICK EVENT
    variable.trace("w", lambda *args: set_databse(variable.get()))
    # CALLING THE SET OIN TEH DEFULT VALUE
    set_databse(variable.get())
    # DECLARING A LABEL
    Label(choosing_DB_screen, text="Choose database ", bg="blue", width="300", height="2", font=("Calibri", 13)).pack()
    # DECLKARE OPTION MENUE OPTIONS
    w = OptionMenu(choosing_DB_screen, variable, "covid_vaccine_new", "covidpublished")
    w.pack()
    # DECLARE THE BUTTON TO PRESS
    Button(choosing_DB_screen, text="choose Subject", height="2", width="30",
           command=lambda: choose_subject(choosing_DB_screen, variable.get())).pack()


#####################################################################################################################
# @DESCRIPTION : THIS METHOD INVOKE THE EXPORT METHOD
# @INPUT : ---
# @OUTPUT : ---
#####################################################################################################################
def export_preproccessed_data():
    main.export_preproccessed_data()

#####################################################################################################################
# @DESCRIPTION : THIS METHOD INVOKE THE IMPORT METHOD
# @INPUT : ---
# @OUTPUT : ---
#####################################################################################################################
def import_preproccessed_data():
    main.import_preproccessed_data()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def change_analyzing_points(chosen_analyzed_sector, select_position_screen, w):
    # CLEAR ALL PREV DECLARES ON THE WINDOW
    w.pack_forget()
    # SET THE CHOSEN ANALYZING APPROACH
    user.selected_analyzed_position = chosen_analyzed_sector
    # VEIW OOPTIONS BASED ON SELECTED DB
    if ('Select Region') in str(chosen_analyzed_sector):
        # ID WE SELECTED REGION THEN WE MUST SET THE CDR LENGTH AND THE START AND END POSITION
        # CALLING A METHOD TO ASSIGN THE CDR LENGTH
        set_cdr_length(0)
        # CALLING A METHOD TO ASSIGN THE END POINT
        set_endposition(0)
        # CALLING A METHOD TO ASSIGN THE START POINT
        set_startposition(0)
        # DECLARING THE VARIABLE
        variable = StringVar(select_position_screen)
        variable.set("'Select Region'")  # default value
        variable.trace("w", lambda *args: set_cdr_length(variable.get()))
        Label(select_position_screen, text="Select Region", bg="blue", width="300", height="2",
              font=("Calibri", 13)).pack()
        Label(select_position_screen, text="Please choose only CDR3", bg="blue", width="300", height="2",
              font=("Calibri", 6)).pack()
        w = OptionMenu(select_position_screen, variable, "CDR1", "CDR2", "CDR3", "FWR1", "FWR2", "FWR3", )
        w.pack()

    if ('Select Start And End Position') in str(chosen_analyzed_sector):
        set_cdr_length(0)
        set_endposition(0)
        variable2 = StringVar(select_position_screen)
        set_startposition(0)
        Label(select_position_screen, text="Select Start Position", bg="blue", width="20", height="2",
              font=("Calibri", 10)).pack()
        variable1 = StringVar(select_position_screen)
        variable1.set("0")  # default value
        variable1.trace("w", lambda *args: set_startposition(variable1.get(),variable1))
        e1 = Entry(select_position_screen, textvariable=variable1)
        e1.pack()
        Label(select_position_screen, text="Select End Position", bg="blue", width="20", height="2",
              font=("Calibri", 10), ).pack()
        variable2.set("0")  # default value
        variable2.trace("w", lambda *args: set_endposition(variable2.get(),variable2))
        e2 = Entry(select_position_screen, textvariable=variable2)
        e2.pack()

    if ('All The Sequence') in str(chosen_analyzed_sector):
        set_cdr_length(0)
        set_endposition(0)
        set_startposition(0)
        Label(select_position_screen, text="analyzing the whole sequence", bg="blue", width="300", height="2",
              font=("Calibri", 13)).pack()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def select_position(choose_subject_screen):
    select_position_screen = Toplevel(choose_subject_screen)
    # choose_subject_screen.withdraw()
    select_position_screen.title("select position")
    select_position_screen.geometry("300x250")
    Label(select_position_screen, text="Select Position", bg="blue", width="300", height="2",
          font=("Calibri", 13)).pack()
    variableStart = StringVar(select_position_screen)
    variableStart.set("Select Region")  # default value
    variableStart.trace("w", lambda *args: set_selected_analyzed_position(variableStart.get()))
    variable = StringVar(select_position_screen)
    variable.set("Select Region")  # default value
    Label(select_position_screen, text="Choose analyzing approach ", bg="blue", width="300", height="2",
          font=("Calibri", 13)).pack()
    w = OptionMenu(select_position_screen, variable, "Select Region", "Select Start And End Position",
                   "All The Sequence")
    w.pack()
    variable.trace("w", lambda *args: change_analyzing_points(variable.get(), select_position_screen, w))
    Button(select_position_screen, text="Import Networks Data", height="2", width="30",
           command=lambda: select_row_count(select_position_screen)).pack()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def select_row_count(select_position_screen):
    select_row_count_screen = Toplevel(select_position_screen)
    # choose_subject_screen.withdraw()
    select_row_count_screen.title("select rows count")
    select_row_count_screen.geometry("300x250")
    variable = StringVar(select_row_count_screen)
    variable.set("all")  # default value
    Label(select_row_count_screen, text="select rows count ", bg="blue", width="300", height="2",
          font=("Calibri", 13)).pack()
    w = OptionMenu(select_row_count_screen, variable, "all", "1000", "2000", "5000", "10000", "20000", "35000", "50000",
                   "80000", "100000")
    w.pack()
    variable.trace("w", lambda *args: set_row_count(variable.get()))
    Button(select_row_count_screen, text="fetch data", height="2", width="30",
           command=lambda: fetch_data(select_row_count_screen)).pack()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def fetch_data(choose_subject_screen):
    fetch_data_screen = Toplevel(choose_subject_screen)
    choose_subject_screen.withdraw()
    fetch_data_screen.title("Fetch Data")
    fetch_data_screen.geometry("300x250")
    Label(fetch_data_screen, text="Fetch Data", bg="blue", width="300", height="2",
          font=("Calibri", 13)).pack()
    Button(fetch_data_screen, text="Fetch", height="2", width="30", command=lambda: start_fetching_data()).pack()
    Button(fetch_data_screen, text="build netwroks", height="2", width="30",
           command=lambda: build_network(fetch_data_screen)).pack()
    Button(fetch_data_screen, text="Export Networks Data", height="2", width="30",
           command=lambda: export_preproccessed_data()).pack()
    Button(fetch_data_screen, text="Import Networks Data", height="2", width="30",
           command=lambda: import_preproccessed_data()).pack()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def build_aa_netwrok():
    try:
        main.create_amino_acid_netwrok()
    except Exception as ex:
        print("ops" + ex)

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def build_nucleotides_network():
    main.create_nuclotides_netwrok()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def build_replaced_nucleotides_network():
    main.replace_nucleotide(user.database, user.subject, user.cdr_length, user.selected_analyzed_position,
                            user.startposition, user.endposition, user.rows)
    main.create_replaced_nuclotides_netwrok()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def show_netwrpk_analyzing_process_results():
    main.show_Nuclotides_netwrok_stats(main.nucoltides_kmers_nets)

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def build_network(fetch_data_screen):
    choosing_network_screen = Toplevel(fetch_data_screen)
    # main_screen.withdraw()
    choosing_network_screen.title("build network")
    choosing_network_screen.geometry("300x250")
    Label(choosing_network_screen, text="Choose network ", bg="blue", width="300", height="2",
          font=("Calibri", 13)).pack()
    Button(choosing_network_screen, text="build AA network", height="2", width="30",
           command=lambda: build_aa_netwrok()).pack()
    Button(choosing_network_screen, text="Build Nucleotides network", height="2", width="30",
           command=lambda: build_nucleotides_network()).pack()
    Button(choosing_network_screen, text="Build most frequented Nucleotides network", height="2", width="30",
           command=lambda: build_replaced_nucleotides_network()).pack()

#####################################################################################################################
# @DESCRIPTION : ASSIGN TO THE USER METADATA THE REQUESTED SUBJECT WINDOW
# @INPUT : choosing_DB_screen: THE EARLIER WINDOW, db: CHOSEN DATABASE
# @OUTPUT : ---
#####################################################################################################################
def main_screen():
    main_screen = Tk()  # create a GUI window
    main_screen.geometry("300x250")  # set the configuration of GUI window
    main_screen.title("Start")  # set the title of GUI window

    # create Login Button
    Button(text="Start analyzing", height="2", width="30", command=lambda: choose_database(main_screen)).pack()
    Label(text="").pack()

    # create a register button
    Button(text="Creators", height="2", width="30", command=lambda: register(main_screen)).pack()

    main_screen.mainloop()


main_screen()
