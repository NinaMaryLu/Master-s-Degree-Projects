class Peak:
    # the constructor for the object
    def __init__(self, SampleID = "", metabolite1 =""):
        # To avoid unintentional replacing of data in the database, copies of attributes are used in all data manipulations
        self._SampleID = SampleID
        self._metabolite1 = metabolite1
    
    @property
    def SampleID(self):
        return self._SampleID
    
    @property
    def metabolite1(self):
        return self._metabolite1
    
    # function for returning the values for the subject as a list of values
    # left it as a function for consistency and possibility to scale it up
    def give_peak(self):
        row_values = [self.SampleID, self.metabolite1]
        return row_values
    
    """
    SETTERS for the object attributes
    To validate if the values match certain criteria
    """
    @SampleID.setter
    def SampleID(self, value):
        # Value has to be a string consisting of anything really. Has to be unique. Can't be empty.
        if isinstance(value, str) and value is not None:
            self._SampleID = value
        else:
            raise ValueError(f"Error found for the SampleID value of ID {self.SampleID}. SampleID must be a string.")
        
    @metabolite1.setter
    def metabolite1(self, value):
        # Value has to be a string consisting of anything really. Has to be unique. Can't be empty.
        if isinstance(value, float) and value is not None:
            self._metabolite1 = value
        else:
            raise ValueError(f"Error found for the metabolite1 value of ID {self.SampleID}. metabolite1 must be a float.")

    """
    Gather_peak(): getter for the attributes not yet had - to source it from the csv file
        working with open(input) file:
            reading the csv file line by line:
                unifying the missing values as "" an empty string,
                extracting the information from consequtive columns separated by "," (split),
                removing any white spaces at the end (rstrip),
                returning the information as subject attributes
    """
    def gather_peak(input):
        line = next(input)          # to skip the header
        for line in input:
            line = line.rstrip().split("\t")
            SampleID = line[0]
            metabolite1 = line[1]
            yield Peak(SampleID, metabolite1)