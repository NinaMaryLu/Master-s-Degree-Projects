import sqlite3
import logging
# setting a logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()

# Class: methods for database manipulation ie. creation, updating, query, data input
class DatabaseManager:
    def __init__(self, db_file):
        """  db_path - a path for the database file. In a working directory, this should be merely a file name ending with .db """
        self._db_file = db_file       # Name has to be a string

    @property
    def db_file(self):
        return self._db_file
    
    def db_query(self, sql):
        """
        sql - sql SELECT command to be executed, str
        """
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        cur.execute(sql)
        # Retrieves all results of a query
        results = cur.fetchall()
        cur.close()
        connection.close()
        return results

    def db_istable(self, name):
        """
        sql - sql SELECT command to be executed, str
        """
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        sql = "SELECT name FROM sqlite_master WHERE type='table' AND name='{0}'"
        cur.execute(sql.format(name))
        # Retrieves all results of a query
        results = cur.fetchall()
        cur.close()
        connection.close()
        return results

    def db_insertmany(self, sql, params):
        """
        sql - sql INSERT command to be executed, str
        params - list of the parameters to be used as attributes for the query
        """
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        try:
            cur.executemany(sql, params)
            connection.commit()
        except sqlite3.IntegrityError:
            logger.warning("There is already data with the unique attribute that you are trying to insert.")
    
    def db_insert(self, sql, params):
        """
        sql - sql INSERT command to be executed, str
        params - list of the parameters to be used as attributes for the query
        """
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        try:
            cur.execute(sql, params)
            connection.commit()
        except sqlite3.IntegrityError:
            logger.warning("There is already data with the unique attribute that you are trying to insert.")
    
    def create_tables(self, entity):
        """
        sql - sql CREATE command to be executed, str
        """
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        if entity == "SUBJECT":
            sql = """CREATE TABLE SUBJECT(
            SubjectID VARCHAR(255), 
            Race VARCHAR(255), 
            Sex VARCHAR(255), 
            Age INT,
            BMI INT, SSPG INT, 
            Insulin VARCHAR(255) CHECK (Insulin in ('IR','IS', NULL)), 
            PRIMARY KEY (SubjectID)
            );"""

        elif entity == "SAMPLE":
            sql = """CREATE TABLE SAMPLE(
            SubjectID VARCHAR(255), 
            VisitID VARCHAR(255), 
            SampleID VARCHAR(255), 
            PRIMARY KEY (SampleID),
			FOREIGN KEY (SubjectID) REFERENCES SUBJECT(SubjectID)
            )"""

        elif entity == "PEAK":  
            sql = """CREATE TABLE PEAK(
            SampleID VARCHAR(255),
            metabolite1 INT,
            PRIMARY KEY (SampleID)
            )"""
        # To import the whole table, one has to either convert it to the long format OR make each metabolite ID an attribute. 
        # To do that, make a list of all the IDs from the column headers, create an object similarly to the Subject. 
        # Then follow the same logic as for INPUT for Subject and pass it to excutemany as a list of parameters
        # (smilarly to the INPUT for Subject) for the CREATE statement, for the names of attributes.
        # The whole table is not used later, so it's not imported as a whole neither.
        elif entity == "PeakAnnotated":
            sql = """CREATE TABLE PeakAnnotated(
            PeakID VARCHAR(255), 
            MetaboliteName VARCHAR(255), 
            PRIMARY KEY (PeakID, MetaboliteName),
			FOREIGN KEY (PeakID) REFERENCES PEAK(PeakID),
			FOREIGN KEY (MetaboliteName) REFERENCES ANNOTATION(MetaboliteName)
            )"""
        elif entity == "ANNOTATION":
            sql = """CREATE TABLE ANNOTATION(
            PeakID VARCHAR(255), 
            MetaboliteName VARCHAR(255),
			KEGG VARCHAR(255),
			HMDB VARCHAR(255),
			ChemicalClass VARCHAR(255),
			Pathway VARCHAR(255),
            PRIMARY KEY (PeakID, MetaboliteName),
			FOREIGN KEY (PeakID) REFERENCES PEAK(PeakID)
            )"""
        elif entity == "TRANSCRIPTOME":
            sql = """CREATE TABLE TRANSCRIPTOME(
            SampleID VARCHAR(255), 
            A1BG VARCHAR(255),
            PRIMARY KEY (SampleID)
            )"""
        # Data not needed for the query
        #elif table == "Proteome":
        #    sql =
        # print(sql)
        cur.execute(sql)
        cur.close()
        connection.close()