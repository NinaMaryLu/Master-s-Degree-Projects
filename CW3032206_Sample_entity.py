class Sample:
    # the constructor for the object
    def __init__(self, SubjectID ="", VisitID="", SampleID=""):
        # To avoid unintentional replacing of data in the database, copies of attributes are used in all data manipulations
        self._SubjectID = SubjectID
        self._VisitID = VisitID
        self._SampleID = SampleID

    @property
    def SubjectID(self):
        return self._SubjectID
    
    @property
    def VisitID(self):
        return self._VisitID
    
    @property
    def SampleID(self):
        return self._SampleID
    
    # function for returning the values for the subject as a list of values
    def give_sample(self):
        row_values = [self.SubjectID, self.VisitID, self.SampleID]
        return row_values
    
    """
    SETTERS for the object attributes
    To validate if the values match certain criteria
    """
    @SubjectID.setter
    def SubjectID(self, value):
        # Value has to be a string consisting of anything really. Has to be unique. Can't be empty.
        if isinstance(value, str) and value is not None:
            self._SubjectID = value
        else:
            raise ValueError(f"Error found for the PeakID value of ID {self.SampleID}:{self.VisitID}. PeakID must be a string.")

    @VisitID.setter
    def VisitID(self, value):
        # Value has to be a string consisting of anything really.
        if isinstance(value, str):
            self._VisitID = value
        else:
            raise ValueError(f"Error found for the MetaboliteName of ID {self.SampleID}:{self.VisitID}. Metabolite name must be a string.")

    @SampleID.setter
    def SampleID(self, value):
        # Value has to be a string
        if isinstance(value, str):
            self._SampleID = value
        else:
            raise ValueError(f"Error found for the KEGG of ID {self.SampleID}:{self.VisitID}. KEGG must be a string.")

    def make_unique(input):   # input should be a list of SampleID as str
        MySet = set(input)
        uniqueList = list(MySet)
        return uniqueList
