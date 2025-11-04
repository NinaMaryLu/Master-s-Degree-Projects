import logging
import sqlite3

from CHEM3032206_Molecule import Molecule
from CHEM3032206_DatabaseManager import DatabaseManager

# setting a logger
logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()

# fetch all

def filter_Lipinski(filter, database):
    if filter == True:
        #conditions = sql_filter_byvalue("Lipinski", True)
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        """
    if filter == False:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        """
    results = database.db_query(sql)
    return results

def filter_LeadLikeness(filter, database):
    if filter == True:
        sql = """SELECT *
        FROM MOLECULE
        WHERE LeadLikeness = True
        """
    if filter == False:
        sql = """SELECT *
        FROM MOLECULE
        WHERE LeadLikeness = False
        """
    results = database.db_query(sql)
    return results

def filter_Bioavailability(filter, database):
    if filter == True:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Bioavailability = True
        """
    if filter == False:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Bioavailability = False
        """
    results = database.db_query(sql)
    return results

def sql_request(request):
    """gets a sql request based on the filters"""
    
    if request == ["All", "All", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        """
    # for Lipinski rule only
    elif request == ["Only TRUE", "All", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        """
    elif request == ["Only FALSE", "All", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        """
    # for lead likeness only
    elif request == ["All", "Only TRUE", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE LeadLikeness = True
        """
    elif request == ["All", "Only FALSE", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE LeadLikeness = False
        """
    # for bioavailability only
    elif request == ["All", "All", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Bioavailability = True
        """
    elif request == ["All", "All", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Bioavailability = False
        """
    # for Lipinski + lead likeness same direction
    elif request == ["Only TRUE", "Only TRUE", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        AND LeadLikeness = True
        """
    elif request == ["Only FALSE", "Only FALSE", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        AND LeadLikeness = False
        """
    # for Lipinski + lead likeness opposite direction
    elif request == ["Only TRUE", "Only FALSE", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        AND LeadLikeness = False
        """
    elif request == ["Only FALSE", "Only TRUE", "All"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        AND LeadLikeness = True
        """
    # for Lipinski + bioavailability same direction
    elif request == ["Only TRUE", "All", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        AND Bioavailability = True
        """
    elif request == ["Only FALSE", "All", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        AND Bioavailability = False
        """
    # for Lipinski + bioavailability opposite direction
    elif request == ["Only TRUE", "All", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        AND Bioavailability = False
        """
    elif request == ["Only FALSE", "All", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        AND Bioavailability = True
        """
        # for lead likeness + bioavailability same direction
    elif request == ["All", "Only TRUE", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE LeadLikeness = True
        AND Bioavailability = True
        """
    elif request == ["All", "Only FALSE", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE LeadLikeness = False
        AND Bioavailability = False
        """
    # for lead likeness + bioavailability opposite direction
    elif request == ["All", "Only TRUE", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE LeadLikeness = True
        AND Bioavailability = False
        """
    elif request == ["All", "Only FALSE", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE LeadLikeness = False
        AND Bioavailability = True
        """
    # combinations of 3: same direction
    elif request == ["Only TRUE", "Only TRUE", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        AND LeadLikeness = True
        AND Bioavailability = True
        """
    elif request == ["Only FALSE", "Only FALSE", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        AND LeadLikeness = False
        AND Bioavailability = False
        """
    # combinations of 3: one different than others
    elif request == ["Only FALSE", "Only TRUE", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        AND LeadLikeness = True
        AND Bioavailability = True
        """
    elif request == ["Only TRUE", "Only False", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        AND LeadLikeness = False
        AND Bioavailability = True
        """
    elif request == ["Only TRUE", "Only TRUE", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        AND LeadLikeness = True
        AND Bioavailability = False
        """
    # combinations of 3: two different than others
    elif request == ["Only FALSE", "Only FALSE", "Only TRUE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        AND LeadLikeness = False
        AND Bioavailability = True
        """
    elif request == ["Only TRUE", "Only False", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = True
        AND LeadLikeness = False
        AND Bioavailability = False
        """
    elif request == ["Only FALSE", "Only TRUE", "Only FALSE"]:
        sql = """SELECT *
        FROM MOLECULE
        WHERE Lipinski = False
        AND LeadLikeness = True
        AND Bioavailability = False
        """
    return sql

def compound_filter_query(db_file, sql):
    """gets full query based on the sql_request result"""
    connection = sqlite3.connect(db_file)
    cur = connection.cursor()
    cur.execute(sql)
    # Retrieves all results of a query
    results = cur.fetchall()
    cur.close()
    connection.close()
    return results

