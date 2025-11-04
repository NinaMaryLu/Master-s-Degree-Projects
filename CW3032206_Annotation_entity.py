class Annotation:
    # the constructor for the object
    def __init__(self, PeakID ="", MetaboliteName="", KEGG="", HMDB="", ChemicalClass="", Pathway=""):
        # To avoid unintentional replacing of data in the database, copies of attributes are used in all data manipulations
        self._PeakID = PeakID
        self._MetaboliteName = MetaboliteName
        self._KEGG = KEGG
        self._HMDB = HMDB
        self._ChemicalClass = ChemicalClass
        self._Pathway = Pathway

    @property
    def PeakID(self):
        return self._PeakID
    
    @property
    def MetaboliteName(self):
        return self._MetaboliteName
    
    @property
    def KEGG(self):
        return self._KEGG
    
    @property
    def HMDB(self):
        return self._HMDB
    
    @property
    def ChemicalClass(self):
        return self._ChemicalClass
    
    @property
    def Pathway(self):
        return self._Pathway
    
    # function for returning the values for the subject as a list of values
    def give_annotation(self):
        row_values = [self.PeakID, self.MetaboliteName, self.KEGG, self.HMDB, self.ChemicalClass, self.Pathway]
        return row_values
    
    """
    SETTERS for the object attributes
    To validate if the values match certain criteria
    """
    @PeakID.setter
    def PeakID(self, value):
        # Value has to be a string consisting of anything really. Has to be unique. Can't be empty.
        if isinstance(value, str) and value is not None:
            self._SubjectID = value
        else:
            raise ValueError(f"Error found for the PeakID value of ID {self.PeakID}:{self.MetaboliteName}. PeakID must be a string.")

    @MetaboliteName.setter
    def MetaboliteName(self, value):
        # Value has to be a string consisting of anything really.
        if isinstance(value, str):
            value = value.replace("(1)","").replace("(2)","").replace("(3)","").replace("(4)","").replace("(5)","")
            self._MetaboliteName = value
        else:
            raise ValueError(f"Error found for the MetaboliteName of ID {self.PeakID}:{self.MetaboliteName}. Metabolite name must be a string.")

    @KEGG.setter
    def KEGG(self, value):
        # Value has to be a string
        if isinstance(value, str):
            self._KEGG = value
        else:
            raise ValueError(f"Error found for the KEGG of ID {self.PeakID}:{self.MetaboliteName}. KEGG must be a string.")

    @HMDB.setter
    def HMDB(self, value):
        # Value has to be a str
        if isinstance(value, str):
            self._HMDB = value
        else:
            raise ValueError(f"Error found for the HMDB of ID {self.PeakID}:{self.MetaboliteName}. HMDB must be a string.")

    @ChemicalClass.setter
    def ChemicalClass(self, value):
        # Value has to be a str
        if isinstance(value, str):
            self._ChemicalClass = value
        else:
            raise ValueError(f"Error found for the ChemicalClass of ID {self.PeakID}:{self.MetaboliteName}. ChemicalClass must be a string.")

    @Pathway.setter
    def Pathway(self, value):
        # Value has to be a str
        if isinstance(value, str):
            self._Pathway = value
        else:
            raise ValueError(f"Error found for the Pathway of ID {self.PeakID}:{self.MetaboliteName}. Pathway must be a string.")

    """
    Gather_subject(): getter for the attributes not yet had - to source it from the csv file
        working with open ("Subject.csv") file:
            reading the csv file line by line:
                unifying the missing values as "" an empty string,
                extracting the information from consequtive columns separated by "," (split),
                removing any white spaces at the end (rstrip),
                returning the information as subject attributes
    """
    def gather_annotation(input):
        line = next(input)          # to skip the header
        for line in input:
            list = ["NA", "Unknown", "unknown", "None", "NULL"]
            line = line.rstrip().split(",")
            new_line = [x if x not in list else None for x in line]
            PeakID, MetaboliteName, KEGG, HMDB, ChemicalClass, Pathway = new_line
            MetaboliteName = MetaboliteName.replace("(1)", "").replace("(2)", "").replace("(3)", "").replace("(4)", "").replace("(5)", "")
            yield Annotation(PeakID, MetaboliteName, KEGG, HMDB, ChemicalClass, Pathway)


#dbFile = "HMP_metabolome_annotation(1).csv"
##Annotations = Annotation(dbFile)      # _input = input = Subject_Data
#assemble_Annotations = []
#with open(dbFile) as source:
#    for record in Annotation.gather_annotation(source):
#        row = record.give_annotation()
#        print(row)
