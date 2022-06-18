# PCN-DNA-sequence-anlyzing
analyzing DNA sequences using PCN and graph theory 


Introduction:
Getting started:
System requirements
•	 Python3 installed
•	 Internet connection
•	 install python pip
Third party components
•	 Pajek
	

Cloning the project:
1.First clone the repo frpm this github repo.
2.Open the editor programming software that you prefer ( we recommend pycharm ).
3.Open the terminal and write this command
	pip install -r requirements.txt

4.Start altering the program code

Upgrading libraries:
First you have to upgrade the library then indicate the version of the library in the requirements file.
recommended use  pip3 freeze > requirements.txt  # Python3
		    pip freeze > requirements.txt  # Python2


Herichia of the program:
the herichia of the program is a simple gui db connection and a logic(main) packages
each one is responable for its meanted naming pourpse.

The gui file is for showing the user interface.
The main file is for creating the networks.
The db copnnection file ios for fetching data from th remote database.
The application starts from the main_screen function in GUI.



Note:

to see the database you can look at the png files in the developer guide attached to this project.

For analyzing the result networks using Pajek:
Download Pajek from: http://mrvar.fdv.uni-lj.si/pajek/
