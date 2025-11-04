# Import modules
# big libraries
import logging
import sqlite3
import pandas as pd
import os

# for chem
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors
#from rdkit.Chem.Draw import IPythonConsole
#from rdkit.Chem.Draw import rdMolDraw2D

# others
import io
from io import BytesIO
from PIL import Image

# for the app
import customtkinter
from customtkinter import CTk
from customtkinter import filedialog
from CustomTkinterMessagebox import CTkMessagebox
from tkinter import *
from tkinter import ttk

# from other py files
from CHEM3032206_Molecule import Molecule
from CHEM3032206_DatabaseManager import DatabaseManager
import CHEM3032206_SQLcommands
import CHEM3032206_TableManager

# colours
button_bg_colour = "#105714"
title_colour = "#88AC4D"
heading_bg_colour = "#445626"

# function to open the app in the centre based on: https://www.geeksforgeeks.org/how-to-center-a-window-on-the-screen-in-tkinter/
def screen_centre(window):
    my_width = window.winfo_screenwidth()
    my_height = window.winfo_screenheight()
    pos_x = (my_width - window.winfo_reqwidth()) // 3           # calculating the position in the middle of the window
    pos_y = (my_height - window.winfo_reqheight()) // 3         # calculating the position on the y axis in the middle of the window
    window.geometry(f"+{pos_x}+{pos_y}")

def get_filter(combobox):
        """ retriving the filter conditions"""
        filterchoice = combobox.get()
        return filterchoice

# setting a logger
logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()

customtkinter.set_appearance_mode("dark")  # Modes: system (default), light, dark
customtkinter.set_default_color_theme("green") 

#DbName = StringVar()

class App(CTk):
    def __init__(self, title):       
        # main setup
        super().__init__()
        self.title(title)
        self.geometry('700x450')
        self.minsize(700, 450)
        # widgets
        self.main = Main(self)  # runs the menu class, thus opening all the widgets
        screen_centre(self)
        self.mainloop()
        self.second_window = None
    
class Main(customtkinter.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent)
        self.place(x = 0, y = 0, relwidth=1, relheight=1)   # rel width and height to adjust the window
        self.mainframe() # will always run the following function when the class is called
        
        frame_upload = customtkinter.CTkFrame(self)
        frame_create = customtkinter.CTkFrame(self)
        frame_filter = customtkinter.CTkFrame(self)
        
        self.uploadframe(frame_upload) # will always run the following function when the class is called
        self.createframe(frame_create)
        self.filterframe(frame_filter)

    def preFilter(self, db_file, combobox_Lipinski, combobox_leadlikeness, combobox_bioavailability):
        """refreshes the table view based on the filters"""
        LipinskiRequest = combobox_Lipinski
        LikenessRequest = combobox_leadlikeness
        BioavRequest = combobox_bioavailability

        if (LipinskiRequest == None) or (LikenessRequest == None) or (BioavRequest == None):
            CTkMessagebox.messagebox(title='Error', text="Can't filter without conditions. Please choose 'All' if you wish to ignore any of them", sound="on", button_text="Oh, OK")
        else:
            request = [LipinskiRequest, LikenessRequest, BioavRequest]
            sql = CHEM3032206_SQLcommands.sql_request(request)
            dataframe = CHEM3032206_SQLcommands.compound_filter_query(db_file, sql)
            return dataframe

    def show_SecondWindow(self, db_file, combobox_Lipinski, combobox_leadlikeness, combobox_bioavailability):
        dataframe = self.preFilter(db_file, combobox_Lipinski, combobox_leadlikeness, combobox_bioavailability)
        #self.app2_window = App2(self,"CHEM5042 Databases, 2", (1100, 400), dataframe, dataframe)
        self.second_window = App2(self,"CHEM5042 Databases", db_path, dataframe, dataframe)
    
    # db_file
    def create_database(self, db_file, sdfFile, label_CreateStatus):
        
        """creating the database"""
        # check if the user listened and actually input the name of the database with .db
        if (db_file.endswith(".db")==False) and (db_file!=""):
            db_file = db_file + ".db"
        if db_file.endswith(".db")==True:
            pass
        elif db_file == "":
            CTkMessagebox.messagebox(title='Error', text="Please enter the name for the database and try again.", sound="on", button_text="Oh, OK")

        if db_file.endswith(".db"):
            DbName = StringVar()
            DbName.set(db_file)
            # "create table is not exists"
            sql_create = """CREATE TABLE if not exists MOLECULE(
            CdId VARCHAR(255),
            Mol_ID VARCHAR(255),
            Image VARCHAR(255),
            SMILES VARCHAR(255),
            Formula VARCHAR(255),
            Name VARCHAR(255),
            MolWeight DECIMAL(12, 3),
            LogP DECIMAL(12, 3),
            LogD DECIMAL(12, 3),
            HydroDonor INT,
            HydroAccept INT,
            RotatableBonds INT,
            NumberRings INT,
            Lipinski BOOLEAN,
            LeadLikeness BOOLEAN,
            Bioavailability BOOLEAN,
            PRIMARY KEY (CdID)
            );"""

            # setting up a database
            try:
                database = DatabaseManager(db_file)
            except:
                logging.error("Couldn't create the database. Python out.")
                CTkMessagebox.messagebox(title='Error', text="Couldn't create the database.", sound="on", button_text="Oh, OK")
            status = database.create_tables(sql_create, "MOLECULE")  # try / error loop is within that function.
            label_CreateStatus.configure(text = status)

            # gathering the moelcules from the sdf file and importing them to the database
            try:
                with Chem.SDMolSupplier(sdfFile) as input:
                    params = []
                    for molecule in input:
                        if molecule is None: continue
                        #print(list(input.GetPropNames())) # gets a list of all available properties

                        mol_input = Molecule()
                        # assigning attributes that are available
                        mol_input.CdId = molecule.GetProp("CdId")
                        mol_input.MolWeight = float(molecule.GetProp("Mol Weight"))
                        mol_input.Formula = molecule.GetProp("Formula")
                        mol_input.Name = molecule.GetProp("Name")
                        mol_input.Mol_ID = molecule.GetProp("Mol_ID")
                        mol_input.LogD = float(molecule.GetProp("LogD"))

                        # assigning attributes that are derived by methods
                        mol_input.Image = AllChem.Compute2DCoords(molecule)
                        mol_input.SMILES = Chem.MolToSmiles(molecule)
                        mol_input.LogP = round(Crippen.MolLogP(molecule), 3)
                        mol_input.HydroDonor = Lipinski.NumHDonors(molecule)
                        mol_input.HydroAccept = Lipinski.NumHAcceptors(molecule)
                        mol_input.RotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)
                        mol_input.NumberRings = rdMolDescriptors.CalcNumRings(molecule)
                        mol_input.Lipinski = Molecule.validate_Lipinski(mol_input)

                        # print(mol_input.NumberRings)
                        # visually checking if the Lipinski value seems correct 
                        # print(mol_input.Lipinski)
                        # print(f"\t{mol_input.HydroDonor} <= 5\n\t{mol_input.HydroAccept} <=10\n\t{mol_input.MolWeight} <=500\n\t{mol_input.LogP} <=5")

                        mol_input.leadlikeness = Molecule.validate_LeadLikeness(mol_input)
                        mol_input.bioavailability = False

                        params.append(mol_input.give_attributes())

                        mol_image = Chem.Draw.MolToImage(molecule, returnPNG=True)
                        temp = BytesIO()
                        mol_image.save(temp, format="PNG")
                        temp.seek(0) # enables reading the file later in tkinter
                        binary_string = temp.read()  
                        mol_input.Image = binary_string

                    # populating the database
                    try:
                        sql_insert = "INSERT INTO MOLECULE VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
                        database.db_insertmany(sql = sql_insert, params = params)
                    except sqlite3.OperationalError:
                        logger.error(f"The table MOLECULE doesn't exist or database is locked (Operational Error). Please create the table first and run again.")
                        CTkMessagebox.messagebox(title='Error', text=status, sound="on", button_text="Oh, OK")
                        label_CreateStatus.configure(text = "Error: data could not be inserted into the database.")

                    # set the status of database creation as done
                    status = "Creation completed successfully"
                    label_CreateStatus.configure(text = status)
                    CTkMessagebox.messagebox(title='Success', text=status, sound="on", button_text="Fantastic")
                    logger.info(f"Creation completed.")
                    global db_path
                    db_path = os.path.join(os. getcwd(), db_file)

            except FileNotFoundError:
                logger.error(f"File {sdfFile} not found. Please check whether the file placement and the file name given are correct")
                label_CreateStatus.configure(text = "Error: sdf file not found.")
                raise SystemExit

    def select_file_db(self, folderPath_db, label_FilePath_db):
        Selected = filedialog.askopenfilename()         # allows you to open the file using the file directory     
        folderPath_db.set(Selected)                     # sets the full path as chosen above
        FileName_db = os.path.basename(Selected)           # extracts just the file name from the path to display nicely
        FileName_db = "Chosen file: \t"+ FileName_db          # returns a string of a file directory
        label_FilePath_db.configure(text = FileName_db)
        global db_path
        db_path = Selected

    def select_file_sdf(self, folderPath_sdf, label_FilePath_sdf):
        Selected = filedialog.askopenfilename()
        folderPath_sdf.set(Selected)                        # sets the folder path as chosen above
        FileName_sdf = os.path.basename(Selected)             # extracts just the file name from the path to display nicely
        FileName_sdf = "Chosen file: \t"+ FileName_sdf          # returns a string of a file directory
        label_FilePath_sdf.configure(text = FileName_sdf)
    
    # frame 0: main frame of the window
    def mainframe(self):
        # defining all the elements in the main frame
        #button_show = customtkinter.CTkButton(self, text = "Show database", font=("Verdana", 10), fg_color= "#105714")
        label_welcome = customtkinter.CTkLabel(self, text = "CHEM5042", font=("Verdana", 20, "bold"), text_color = "#88AC4D")
        label_welcome2 = customtkinter.CTkLabel(self, text = "Chem database viewer", font=("Verdana", 15), text_color = "#88AC4D")
        # defining the grid of the main frame
        self.columnconfigure((0,1), weight = 1, uniform = "a")
        self.rowconfigure((1,2,3,4), weight = 1, uniform = "a")
        # placing all the elements within that frame
        label_welcome.grid(row = 0, column = 0, sticky = "n", rowspan = 1, columnspan=1)
        label_welcome2.grid(row = 0, column = 1, sticky = "n", rowspan = 1, columnspan=1)
    
    # frame 1: upload the database you already have
    def uploadframe(self, frame):
        # defining all the elements in the frame
        label_UploadTitle = customtkinter.CTkLabel(frame, text = "Open an existing database.", font=("Verdana", 15, "bold"), text_color = "#88AC4D")
        label_FilePath_db = customtkinter.CTkLabel(frame, text = " ", font=("Verdana", 8), text_color = "#D9D9D9") #Chosen file:\t
        button_file_db = customtkinter.CTkButton(frame, text="Choose database file", font=("Verdana", 10, "bold"), text_color= "#88AC4D", fg_color="#F1F5FF", command=lambda: self.select_file_db(folderPath_db, label_FilePath_db))
        folderPath_db = StringVar()    # tkinter variable that is being updated
        DbName = StringVar()
        DbName.set(folderPath_db)
        # defining the grid of the frame
        frame.grid(row = 1, column = 0, sticky = "nwes", rowspan = 2, columnspan = 1, ipadx=10, ipady=10, padx=(10,5), pady=(10,5))
        # placing all the elements within that frame
        frame.columnconfigure((0), weight = 1, uniform = "a")
        frame.rowconfigure((2), weight = 1, uniform = "a")
        label_UploadTitle.grid(row = 0, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=10, pady=(10,5))
        button_file_db.grid(row = 1, column = 0, sticky = "we", rowspan = 1, columnspan = 1, padx=20, pady=(5,5))
        label_FilePath_db.grid(row = 2, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=20, pady=(5,10))
    
    # frame 2: create a database
    def createframe(self, frame):
        # defining all the elements in the frame
        folderPath_sdf = StringVar()    # tkinter variable that is being updated
        FileName_sdf = StringVar()
        DbName = StringVar()
        
        label_CreateTitle = customtkinter.CTkLabel(frame, text = "Create a new database.", font=("Verdana", 15, "bold"),text_color = "#88AC4D")
        label_FilePath_sdf = customtkinter.CTkLabel(frame, text = " ", font=("Verdana", 8), text_color = "#D9D9D9")    #Chosen file: \t
        button_file_sdf = customtkinter.CTkButton(frame, text="Choose sdf file", font=("Verdana", 10, "bold"), text_color= "#88AC4D", fg_color="#F1F5FF", command=lambda: self.select_file_sdf(folderPath_sdf, label_FilePath_sdf))
        label_DbName = customtkinter.CTkLabel(frame, text = "Name for the\ndatabase(with .db)", font=("Verdana", 10), padx = 2)
        entry_DbName = customtkinter.CTkEntry(frame)

        button_StartCreate = customtkinter.CTkButton(frame, text = "Start creating the database.", font=("Verdana", 10), fg_color= title_colour, command=lambda: self.create_database(entry_DbName.get(), folderPath_sdf.get(), label_CreateStatus))
        label_CreateStatus = customtkinter.CTkLabel(frame, text = "Status: ", font=("Verdana", 8), text_color = "#D9D9D9")
        # defining the grid of the frame
        frame.grid(row = 3, column = 0, sticky = "nwes", rowspan = 2, columnspan = 1, ipadx=10, ipady=10, padx=(10,5), pady=(5,10))
        # placing all the elements within that frame
        frame.columnconfigure((0,1), weight = 1, uniform = "a")
        frame.rowconfigure((3,5), weight = 1, uniform = "a")
        label_CreateTitle.grid(row = 0, column = 0, sticky = "nw", rowspan = 1, columnspan = 2, padx=10, pady=(10,5))
        button_file_sdf.grid(row = 1, column = 0, sticky = "nwe", rowspan = 1, columnspan = 2, padx=20, pady=(5,5))
        label_DbName.grid(row = 2, column = 0, sticky = "w", rowspan = 1, columnspan = 2, padx=10, pady=(5,5))
        entry_DbName.grid(row = 2, column = 1, sticky = "e", rowspan = 1, columnspan = 2, padx=20, pady=(5,5))
        label_FilePath_sdf.grid(row = 3, column = 0, sticky = "nw", rowspan = 1, columnspan = 2, padx=20, pady=(5,5))
        button_StartCreate.grid(row = 4, column = 0, sticky = "we", rowspan = 1, columnspan = 2, padx=20, pady=(5,5))
        label_CreateStatus.grid(row = 5, column = 0, sticky = "nw", rowspan = 1, columnspan = 2, padx=20, pady=(5,10))
    
    # frame 3: filter the database to show
    def filterframe(self, frame):
        options = ["Only TRUE", "Only FALSE", "All"]

        # defining all the elements in the frame
        
        label_FilterTitle = customtkinter.CTkLabel(frame, text = "Filter the database before browsing.", font=("Verdana", 15, "bold"), text_color = "#88AC4D")
        
        label_Lipinski = customtkinter.CTkLabel(frame, text = "Lipinski rule", font=("Verdana", 9))
        label_LeadLikeness = customtkinter.CTkLabel(frame, text = "Lead-likeness", font=("Verdana", 9))
        label_Bioavailability = customtkinter.CTkLabel(frame, text = "bioavailability", font=("Verdana", 9))

        combobox_Lipinski_var=customtkinter.StringVar()
        combobox_Lipinski = customtkinter.CTkComboBox(frame, values=options, state="readonly", font=("Verdana", 9), variable=combobox_Lipinski_var)
        combobox_Lipinski.set(options[2])
        combobox_Lipinski.bind("<<ComboboxSelected>>", get_filter(combobox_Lipinski))
        label_Lipinski.grid(row = 1, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=(10,5), pady=10)
        combobox_Lipinski.grid(row = 1, column = 1, sticky = "nw", rowspan = 1, columnspan = 2, padx=(5,10), pady=10)

        combobox_leadlikeness_var=customtkinter.StringVar()
        combobox_leadlikeness = customtkinter.CTkComboBox(frame, values=options, state="readonly", font=("Verdana", 9), variable=combobox_leadlikeness_var)
        combobox_leadlikeness.set(options[2])
        combobox_leadlikeness.bind("<<ComboboxSelected>>", get_filter(combobox_leadlikeness))
        label_LeadLikeness.grid(row = 2, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=(10,5), pady=10)
        combobox_leadlikeness.grid(row = 2, column = 1, sticky = "nw", rowspan = 1, columnspan = 2, padx=(5,10), pady=10)
        
        combobox_bioavailability_var=customtkinter.StringVar()
        combobox_bioavailability = customtkinter.CTkComboBox(frame, values=options, state="readonly", font=("Verdana", 9), variable=combobox_bioavailability_var)
        combobox_bioavailability.set(options[2])
        combobox_bioavailability.bind("<<ComboboxSelected>>", get_filter(combobox_bioavailability))
        label_Bioavailability.grid(row = 3, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=(10,5), pady=10)
        combobox_bioavailability.grid(row = 3, column = 1, sticky = "nw", rowspan = 1, columnspan = 2, padx=(5,10), pady=10)

        #toggle_SaveSelection = customtkinter.CTkCheckBox(frame, text = "Export subset as csv: ", font=("Verdana", 10), onvalue = 1, offvalue = 0)
        #entry_SelectionCsvName = customtkinter.CTkEntry(frame)
        #button_Export = customtkinter.CTkButton(frame, text = "Export Selection as csv.", font=("Verdana", 10), fg_color= title_colour, command=lambda: self.export_selection(entry_SelectionCsvName.get()))
        #label_ExportStatus = customtkinter.CTkLabel(frame, text = "Status: ", font=("Verdana", 8), text_color = "#D9D9D9")
        button_show = customtkinter.CTkButton(frame, text = "Show database", font=("Verdana", 10), fg_color= "#105714", command=lambda: self.show_SecondWindow(db_path, get_filter(combobox_Lipinski), get_filter(combobox_leadlikeness), get_filter(combobox_bioavailability)))

        # defining the grid of the frame
        frame.grid(row = 1, column = 1, sticky = "nwes", rowspan = 4, columnspan = 1, ipadx=10, ipady=10, padx=(5,10), pady=10)
        # placing all the elements within that frame
        frame.columnconfigure((0,1), weight = 1, uniform = "a")
        frame.rowconfigure((1,2,3,4,5,6,7), weight = 1, uniform = "a")
        label_FilterTitle.grid(row = 0, column = 0, sticky = "nw", rowspan = 1, columnspan = 2, padx=10, pady=(10,5))

        #toggle_SaveSelection.grid(row = 4, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=(20,5), pady=(5,5))
        #entry_SelectionCsvName.grid(row = 4, column = 1, sticky = "nw", rowspan = 1, columnspan = 1, padx=(5,20), pady=(5,5))
        #button_Export.grid(row = 5, column = 0, sticky = "wse", rowspan = 1, columnspan = 2, padx=20, pady=(5,5))
        #label_ExportStatus.grid(row = 6, column = 0, sticky = "nw", rowspan = 1, columnspan = 2, padx=20, pady=(5,10))
        button_show.grid(row = 7, column = 0, sticky = "ew", rowspan = 1, columnspan = 2, padx=10, pady=10)

class App2(customtkinter.CTkToplevel):
    def __init__(self, parent, title, path, dataframe, subset):   
        super().__init__(parent)    
        # main setup
        #self = ttk.Toplevel(parent)
        self.title(title)
        self.geometry('1100x400')
        self.minsize(1100, 400)
        self.attributes("-topmost", True)   # window on top of the first
        # widgets
        self.main2 = Main2(self, dataframe, subset)  # runs the menu class, thus opening all the widgets
        #self.mainloop()

class Main2(customtkinter.CTkFrame):
    def __init__(self, parent, dataframe, subset):
        super().__init__(parent)
        self.place(x = 0, y = 0, relwidth=1, relheight=1)   # rel width and height to adjust the window
        self.mainframe() # will always run the following function when the class is called
        #screen_centre(self)
        dataframe_source = dataframe
        database_opened = subset
        
        # defining the frames within
        frame_filter = customtkinter.CTkFrame(self)   # add the scroll if needed (ScrollableFrame)
        frame_filter.grid(row = 1, column = 0, sticky = "nwes", rowspan = 4, columnspan = 1, ipadx=10, ipady=10, padx=(10,5), pady=10)
        frame_view = customtkinter.CTkFrame(self) 
        frame_view.grid(row = 1, column = 1, sticky = "nwes", rowspan = 4, columnspan = 3, ipadx=10, ipady=10, padx=(5,10), pady=10)

        headings = ["CdId", "Mol_ID", "Image", "SMILES", "Formula", "Name", "MolWeight", "LogP", "LogD",
                    "Hydrogen Donors", "Hydrogen Acceptors", "RotatableBonds", "NumberRings", "Lipinski", "leadlikeness", "bioavailability"]
        frame_tree = ttk.Treeview(self, columns=headings, show='headings')
        frame_tree.grid(row = 2, column = 1, sticky = "nwse", rowspan = 3, columnspan = 3, padx=10, pady=(5,10))

        # adding the widgets
        self.filterframe(dataframe_source, frame_filter, frame_tree) #the second is needed to reset the table view
        self.viewframe(frame_view)
        self.tableview(database_opened, frame_tree)

    def refresh(self, db_path, combobox_Lipinski, combobox_leadlikeness, combobox_bioavailability, frame):
        """refreshes the table view based on the filters."""
        LipinskiRequest = combobox_Lipinski.get()
        LikenessRequest = combobox_leadlikeness.get()
        BioavRequest = combobox_bioavailability.get()

        if (LipinskiRequest == None) or (LikenessRequest == None) or (BioavRequest == None):
            CTkMessagebox.messagebox(title='Error', text="Can't filter without conditions. Please choose 'All' if you wish to ignore any of them", sound="on", button_text="Oh, OK")
        else:
            request = [LipinskiRequest, LikenessRequest, BioavRequest]
            sql = CHEM3032206_SQLcommands.sql_request(request)
            dataframe = CHEM3032206_SQLcommands.compound_filter_query(db_path, sql)
            #frame.grid_forget() can't hide the frame (because the latter doesn't appear) but can be covered
            self.tableview(dataframe, frame)

    # frame 0: main frame of the window
    def mainframe(self):
        # defining all the elements in the main frame
        label_welcome = customtkinter.CTkLabel(self, text = "CHEM5042", font=("Verdana", 20, "bold"), text_color = "#88AC4D")
        label_welcome2 = customtkinter.CTkLabel(self, text = "Chem database viewer", font=("Verdana", 15), text_color = "#88AC4D")
        # defining the grid of the main frame
        self.columnconfigure((1,2,3), weight = 1, uniform = "a")
        self.rowconfigure((1,2,3,4), weight = 1, uniform = "a")
        # placing all the elements within that frame
        label_welcome.grid(row = 0, column = 0, sticky = "n", rowspan = 1, columnspan=1)
        label_welcome2.grid(row = 0, column = 1, sticky = "n", rowspan = 1, columnspan=1)
    
    def filterframe(self, db_file, frame, frame2):
        # options for the drop down filters
        options = ["Only TRUE", "Only FALSE", "All"]
        
        # defining the grid of the frame
        frame.columnconfigure((0,1,2), weight = 1, uniform = "a")
        frame.rowconfigure((0,1,2,3), weight = 1, uniform = "a")

        # placing all the elements within that frame
        label_FilterTitle = customtkinter.CTkLabel(frame, text = "Filter the database.", font=("Verdana", 12, "bold"), text_color = "#88AC4D")
        label_FilterTitle.grid(row = 0, column = 0, sticky = "nw", rowspan = 1, columnspan = 2, padx=10, pady=(10,5))

        label_Lipinski = customtkinter.CTkLabel(frame, text = "Lipinski rule", font=("Verdana", 9))
        label_LeadLikeness = customtkinter.CTkLabel(frame, text = "Lead-likeness", font=("Verdana", 9))
        label_Bioavailability = customtkinter.CTkLabel(frame, text = "bioavailability", font=("Verdana", 9))

        combobox_Lipinski_var=customtkinter.StringVar()
        combobox_Lipinski = customtkinter.CTkComboBox(frame, values=options, state="readonly", font=("Verdana", 9), variable=combobox_Lipinski_var)
        combobox_Lipinski.set(options[2])
        combobox_Lipinski.bind("<<ComboboxSelected>>", get_filter(combobox_Lipinski))
        label_Lipinski.grid(row = 1, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=(10,5), pady=10)
        combobox_Lipinski.grid(row = 1, column = 1, sticky = "nw", rowspan = 1, columnspan = 2, padx=(5,10), pady=10)

        combobox_leadlikeness_var=customtkinter.StringVar()
        combobox_leadlikeness = customtkinter.CTkComboBox(frame, values=options, state="readonly", font=("Verdana", 9), variable=combobox_leadlikeness_var)
        combobox_leadlikeness.set(options[2])
        combobox_leadlikeness.bind("<<ComboboxSelected>>", get_filter(combobox_leadlikeness))
        label_LeadLikeness.grid(row = 2, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=(10,5), pady=10)
        combobox_leadlikeness.grid(row = 2, column = 1, sticky = "nw", rowspan = 1, columnspan = 2, padx=(5,10), pady=10)
        
        combobox_bioavailability_var=customtkinter.StringVar()
        combobox_bioavailability = customtkinter.CTkComboBox(frame, values=options, state="readonly", font=("Verdana", 9), variable=combobox_bioavailability_var)
        combobox_bioavailability.set(options[2])
        combobox_bioavailability.bind("<<ComboboxSelected>>", get_filter(combobox_bioavailability))
        label_Bioavailability.grid(row = 3, column = 0, sticky = "nw", rowspan = 1, columnspan = 1, padx=(10,5), pady=10)
        combobox_bioavailability.grid(row = 3, column = 1, sticky = "nw", rowspan = 1, columnspan = 2, padx=(5,10), pady=10)
    
        button_refresh = customtkinter.CTkButton(frame, text = "Refresh view.", font=("Verdana", 10), fg_color= "#105714", command=lambda: self.refresh(db_path, combobox_Lipinski, combobox_leadlikeness, combobox_bioavailability, frame2))
        button_refresh.grid(row = 4, column = 0, sticky = "nwes", rowspan = 1, columnspan = 3, padx=10, pady=10)

    def viewframe(self, frame):
        # defining the grid of the frame
        frame.columnconfigure((0,1), weight = 1, uniform = "a")
        frame.rowconfigure((1), weight = 1, uniform = "a")

        # placing all the elements within that frame
        label_ViewTitle = customtkinter.CTkLabel(frame, text = "View the data.", font=("Verdana", 12, "bold"), text_color = "#88AC4D")
        label_ViewTitle.grid(row = 0, column = 0, sticky = "nw", rowspan = 1, columnspan = 2, padx=10, pady=(10,5))

    def tableview(self, data, frame):
        headings = ["CdId", "Mol_ID", "Image", "SMILES", "Formula", "Name", "MolWeight", "LogP", "LogD",
                    "Hydrogen Donors", "Hydrogen Acceptors", "RotatableBonds", "NumberRings", "Lipinski", "leadlikeness", "bioavailability"]
        
        # styling the table based on the https://github.com/TomSchimansky/CustomTkinter/discussions/524
        mystyle = ttk.Style()
        mystyle.theme_use('default')
        bg_colour = self._apply_appearance_mode(customtkinter.ThemeManager.theme["CTkFrame"]["fg_color"])
        text_colour = self._apply_appearance_mode(customtkinter.ThemeManager.theme["CTkLabel"]["text_color"])
        selected_colour = self._apply_appearance_mode(customtkinter.ThemeManager.theme["CTkButton"]["fg_color"])
        mystyle.configure("Treeview", background=bg_colour, foreground=text_colour, fieldbackground=bg_colour, borderwidth=0)
        mystyle.map('Treeview', background=[('selected', bg_colour)], foreground=[('selected', selected_colour)])
        
        # creating the headings
        mystyle.configure("Treeview.Heading", background=heading_bg_colour, foreground="white")
        CHEM3032206_TableManager.create_headings(frame, headings)
        for col in frame["columns"]:
            frame.column(col, anchor=CENTER)

        # populating the table
        for row in data:
            frame.insert(parent="", index = 0, values = row)

        # horizontal scroll bar
        HorizontalScroll =  ttk.Scrollbar(self, orient="horizontal", command=frame.xview)
        HorizontalScroll.grid(row = 2, column = 1, sticky = "swe", rowspan = 3, columnspan = 3, padx=10, pady=(0,10))
        frame.configure(xscrollcommand=HorizontalScroll.set)
        # change colour if you have time

        # vertical scroll bar
        VerticalScroll =  ttk.Scrollbar(self, orient="vertical", command=frame.yview)
        VerticalScroll.grid(row = 2, column = 1, sticky = "sne", rowspan = 3, columnspan = 3, padx=(0,10), pady=10)
        frame.configure(yscrollcommand=VerticalScroll.set)
        
        # open selection in PubChem
        frame.bind("<<TreeviewSelect>>", lambda e: CHEM3032206_TableManager.TheChosen(frame))

App("CHEM5042 Databases")
