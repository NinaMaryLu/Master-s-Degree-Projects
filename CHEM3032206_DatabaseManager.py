
# Database methods amended from the file used in the Programming and Databases module, written by me for the purpose of the assignment in that course.

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
    
    def db_query_conditions(self, sql, params):
        """
        sql - sql SELECT command to be executed, str
        """
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        cur.execute(sql, params)
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
    
    def create_tables(self, sql, name):
        """
        sql - sql CREATE command to be executed, str
        """
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        try:
            cur.execute(sql)
            status = "Table created."
        except: # if the table already exists, pass and continue creating other tables
            if self.db_istable(name) != None:
                status = f"Table {name} already exists."
                logger.warning(status)
                
            else:
                status = f"Table {name} doesn't exist and yet - it failed to be created."
                logger.warning(status)
        cur.close()
        connection.close()
        return status