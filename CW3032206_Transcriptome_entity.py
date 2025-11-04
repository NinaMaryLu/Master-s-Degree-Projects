class Transcriptome:
    # the constructor for the object
    # participates in the one-to-many relationship so the FK is on the many side
    def __init__(self, SampleID = "", A1BG =""):
        # To avoid unintentional replacing of data in the database, copies of attributes are used in all data manipulations
        self._SampleID = SampleID
        self._A1BG = A1BG
    
    @property
    def SampleID(self):
        return self._SampleID
    
    @property
    def A1BG(self):
        return self._A1BG
    
    # function for returning the values for the subject as a list of values
    def give_transcript(self):
        row_values = [self.SampleID, self.A1BG]
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
            raise ValueError(f"Error found for the SampleID value of sample {self.SampleID}. SampleID must be a string.")

    @A1BG.setter
    def A1BG(self, value):
        # Value has to be a string consisting of anything really. Has to be unique. Can't be empty.
        if isinstance(value, float):
            self._A1BG = value
        else:
            raise ValueError(f"Error found for the A1BG value of sample {self.SampleID}. A1BG must be a float.")

    """
    Gather_transcript(): getter for the attributes not yet had - to source it from the csv file
        working with open ("Subject.csv") file:
            reading the csv file line by line:
                unifying the missing values as "" an empty string,
                extracting the information from consequtive columns separated by "," (split),
                removing any white spaces at the end (rstrip),
                returning the information as subject attributes
    """
    def gather_transcript(input):
        line = next(input)          # to skip the header
        for line in input:
            list = ["NA", "Unknown", "unknown", "None", "NULL"]
            line = line.rstrip().split("\t")
            line = line[0:2]
            new_line = [x if x not in list else None for x in line]
            PeakID, A1BG = new_line
            yield Transcriptome(PeakID, A1BG)
