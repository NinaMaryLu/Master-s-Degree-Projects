# Reference istable: I used a following link to solve how to check if the table already exists in the schema: https://stackoverflow.com/questions/1601151/how-do-i-check-in-sqlite-whether-a-table-exists
# Reference 

# Importint modules
#import time

import logging
import argparse
import sqlite3

# importing classes from other .py files
from CW3032206_database_methods import DatabaseManager
from CW3032206_Subject_entity import Subject
from CW3032206_Annotation_entity import Annotation
from CW3032206_Transcriptome_entity import Transcriptome
from CW3032206_Sample_entity import Sample
from CW3032206_Peak_entity import Peak


# Creating a parser for the script
parser = argparse.ArgumentParser(description= "Arguments needed for the script:")
parser.add_argument('--createdb', required=False, action="store_true", help='Initialises creating the database when put at the command line')
parser.add_argument('--loaddb', required=False, action="store_true", help='Loades data into the database when put at the command line')
parser.add_argument('--querydb', required=False, action="store", help='Performs query of a set number when put at the command line. --querydb=n performs the n-query. There are 9 queries (numbers 1-9).')
parser.add_argument("databasefile", help = "Name of the existing database file / for the file to be created")         # I assume that handles missing file
args = parser.parse_args()

# setting a logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()

## record start time
#start = time.time()


# Functions for the main
#1: Loading data into the subject table
def LoadData_Subject(input):
        # Compiling a list of all subjects called param
    sql = "INSERT INTO SUBJECT VALUES(?, ?, ?, ?, ?, ?, ?)"
    param = []
    with open(input) as source:
        for record in Subject.gather_subject(source):
            row = record.give_subject()
            param.append(row)
    database.db_insertmany(sql, param)

#2: Loading data into the annotation table
def LoadData_ANNOTATION(input):
        # Compiling a list of all subjects called param
    sql = "INSERT INTO ANNOTATION VALUES(?, ?, ?, ?, ?, ?)"
    param = []
    with open(input) as source:
        for record in Annotation.gather_annotation(source):
            row = record.give_annotation()
            param.append(row)
    database.db_insertmany(sql, param)

#4: Loading data into the annotation table
def LoadData_TRANSCRIPTOME(input):
        # Compiling a list of all subjects called param
    sql = "INSERT INTO TRANSCRIPTOME VALUES(?, ?)"
    param = []
    with open(input) as source:
        for record in Transcriptome.gather_transcript(source):
            row = record.give_transcript()
            SampleID_list.append(row[0])
            param.append(row)
    database.db_insertmany(sql, param)
    return SampleID_list

#5: Loading data into the SAMPLE table
def LoadData_SAMPLE(SampleID_list):
    SampleID_UniqueList = Sample.make_unique(SampleID_list)
    sql = "INSERT INTO SAMPLE VALUES(?, ?, ?)"
    param = []
    for item in SampleID_UniqueList:
        Sample.SampleID = item
        item = item.split("-")
        Sample.SubjectID, Sample.VisitID = item
        row = [Sample.SubjectID, Sample.VisitID, Sample.SampleID]
        param.append(row)
    database.db_insertmany(sql, param)

#6: Loading data into the PEAK table
def LoadData_PEAK(input):
    sql = "INSERT INTO PEAK VALUES(?, ?)"
    param = []
    with open(input) as source:
        for record in Peak.gather_peak(source):
            row = record.give_peak()
            param.append(row)
    database.db_insertmany(sql, param)

#7: Query - getting the sql command for specific question
def define_query(question):
    if question == "1":
        sql = """SELECT SubjectID, Age
        FROM SUBJECT
        WHERE Age > 70"""
    elif question == "2":
        sql = """SELECT SubjectID
        FROM SUBJECT
        WHERE Sex = "F" 
        AND BMI BETWEEN 18.5 AND 24.9
        ORDER BY BMI DESC"""
    elif question == "3":
        sql = """SELECT VisitID
        FROM SAMPLE
        WHERE SubjectID ="ZNQOVZV"
        """
    elif question == "4":
        sql = """SELECT DISTINCT SubjectID
        FROM SUBJECT
        WHERE Insulin = "IR" AND SubjectID IN (
        SELECT SubjectID
        FROM SAMPLE
        WHERE SAMPLE.SampleID IN (
        SELECT SampleID
        FROM PEAK
        WHERE PEAK.metabolite1 IS NOT NULL));
        """
    elif question == "5":
        sql = """SELECT DISTINCT KEGG
        FROM ANNOTATION
        WHERE PeakID IN ("nHILIC_121.0505_3.5", "nHILIC_130.0872_6.3", "nHILIC_133.0506_2.3", "nHILIC_133.0506_4.4")
        """

    elif question == "5a":
        sql="""SELECT DISTINCT KEGG
        FROM ANNOTATION
        WHERE PeakID = "nHILIC_121.0505_3.5"
        """
    elif question == "5b":
        sql = """SELECT DISTINCT KEGG
        FROM ANNOTATION
        WHERE PeakID = "nHILIC_130.0872_6.3"
        """
    elif question == "5c":
        sql = """SELECT DISTINCT KEGG
        FROM ANNOTATION
        WHERE PeakID = "nHILIC_133.0506_2.3"
        """
    elif question == "5d":
        sql = """SELECT DISTINCT KEGG
        FROM ANNOTATION
        WHERE PeakID = "nHILIC_133.0506_4.4"
        """
    elif question == "6":
        sql = """SELECT MIN(Age) as SmallestAge
        MAX(Age) as BiggestAge
        AVG(Age) as AverageAge
        FROM SUBJECT;
        """
    elif question == "7":
        sql = """SELECT COUNT(PeakID) AS AnnotationNumber, Pathway
        FROM ANNOTATION
        WHERE AnnotationNumber >=10
        GROUP BY Pathway
        ORDER BY AnnotationCount DESC"""
    elif question == "8":
        sql = """SELECT MAX(A1BG)
        FROM TRANSCRIPT
        WHERE SampleID IN(
        SELECT SampleID
        FROM SAMPLE
        WHERE SubjectID = ZOZOW1T)
        """
    elif question == "9":
        sql = """SELECT Age, BMI
        FROM SUBJECT
        WHERE SubjectID NOT IN (
        SELECT SubjectID
        FROM SUBJECT
        WHERE Age IS NULL OR BMI IS NULL)
        """
    return sql

#8: Query - actually performing the query and getting the answer
def answer_query(sql):
    output = []
    results = database.db_query(sql)
    output = "\n".join([("\t".join([str(x) for x in result_row])) for result_row in results]) # I got help from https://stackoverflow.com/questions/35465468/create-tab-separated-print-output-from-a-list-of-tuples-python
    return output

# path of the database - the positional argument
db_file = args.databasefile   # will be inputted from the command line

# I'm leaving it before the if-else so that it exists for the load, even if the create statements are not passed (assumes that the file of that name exists then)
database = DatabaseManager(db_file)

# if the --createdb present in the command line, go to "else"
if args.createdb == False:
    pass
else:
    # Creating tables
    table_names = ["SUBJECT", "SAMPLE", "PEAK", "PeakAnnotated", "ANNOTATION", "TRANSCRIPTOME"] # no proteome table because it is not needed
    for name in table_names:
        try:
            database.create_tables(name)
        except: # if the table already exists, pass and continue creating other tables
            if database.db_istable(name) != None:
                logger.warning(f"Table {name} already exists")
                continue

# if the --loaddb present in the command line, go to "else"
if args.loaddb == False:
    pass
else:
    # Loading the data from Subject
    Subject_Data = "Subject.csv"
    try:
        LoadData_Subject(Subject_Data)
    except FileNotFoundError:
        logger.error(f"File {Subject_Data} not found. Please check whether the file placement and the file name given are correct")
        raise SystemExit
    except sqlite3.OperationalError:
        logger.error(f"The table SUBJECT doesn't exist or database is locked (Operational Error). Please create the table first and run again.")

    # Loading the data for PeakAnnotated
    # data will be loaded while the Annotated table is being compiled. 
    # row_PeakAnnotated - list of rows for the table
    param_PeakAnnotated = []

    # Loading the data for Annotation
    Annotation_Data = "HMP_metabolome_annotation.csv"
    try:
        LoadData_ANNOTATION(Annotation_Data)
    except FileNotFoundError:
        logger.error(f"File {Annotation_Data} not found. Please check whether the file placement and the file name given are correct")
        raise SystemExit
    except sqlite3.OperationalError:
        logger.error(f"The table ANNOTATION doesn't doesn't exist or database is locked (Operational Error). Please create the table first and run again.")

    # Parsing SubjectID and VisitID from the Sample ID
    # will happen simultaniously to loading the data for the Transcriptome entity: SampleID will be added to a list and then converted to a set to make values unique
    SampleID_list = []

    # Loading the data for Transcriptome
    Transcriptome_Data = "HMP_transcriptome_abundance.tsv"
    try:
        LoadData_TRANSCRIPTOME(Transcriptome_Data)
    except FileNotFoundError:
        logger.error(f"File {Transcriptome_Data} not found. Please check whether the file placement and the file name given are correct")
        raise SystemExit
    except sqlite3.OperationalError:
        logger.error(f"The table TRANSCRIPTOME doesn't exist or database is locked (Operational Error). Please create the table first and run again.")
    # parse SampleID into SubjectID and VisitID and loading the data into the table
    try:
        LoadData_SAMPLE(SampleID_list)
    except sqlite3.OperationalError:
        logger.error(f"The table SAMPLE doesn't exist. Please create the table first.")
    # Loading the data for PEAK
    Peak_Data = "HMP_metabolome_abundance.tsv"
    try:
        LoadData_PEAK(Peak_Data)
    except FileNotFoundError:
        logger.error(f"File {Peak_Data} not found. Please check whether the file placement and the file name given are correct")
        raise SystemExit
    except sqlite3.OperationalError:
        logger.error(f"The table PEAK doesn't exist or database is locked (Operational Error). Please create the table first and run again.")

# if the --querydb is present in the command line
if args.querydb == False:
    pass
else:
    sql = define_query(args.querydb)
    output = answer_query(sql)
    print(output)

# record end time
#end = time.time()

# print the difference between start 
# and end time in milli. secs
#print("The time of execution of above program is :",
#      (end-start) * 10**3, "ms")
