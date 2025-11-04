# Creating the entity Subject as an object with the set attributes as a superclass
class Subject:
    # the constructor for the object
    def __init__(self, SubjectID ="", Race="", Sex="",Age="", BMI="", SSPG="", Insulin=""):
        # To avoid unintentional replacing of data in the database, copies of attributes are used in all data manipulations
        self._SubjectID = SubjectID
        self._Race = Race
        self._Sex = Sex
        self._Age = Age
        self._BMI = BMI
        self._SSPG = SSPG
        self._Insulin = Insulin

    """
    GETTERS for the object attributes
    Those copies are being returned when the attribute is requested
    _ marks an internal implementaton
    """
    
    @property
    def SubjectID(self):
        return self._SubjectID
    
    @property
    def Race(self):
        return self._Race
    
    @property
    def Sex(self):
        return self._Sex
    
    @property
    def Age(self):
        return self._Age
    
    @property
    def BMI(self):
        return self._BMI
    
    @property
    def SSPG(self):
        return self._SSPG
    
    @property
    def Insulin(self):
        return self._Insulin
    
    # function for returning the values for the subject as a list of values
    def give_subject(self):
        subject_values = [self.SubjectID, self.Race, self.Sex, self.Age, self.BMI, self.SSPG, self.Insulin]
        return subject_values

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
            raise ValueError(f"Error found for the SubjectID value of subject {self.SubjectID}. SubjectID must be a string.")

    @Race.setter
    def Race(self, value):
        # Value has to be a string consisting of anything really.
        if isinstance(value, str):
            self._Race = value
        else:
            raise ValueError(f"Error found for the race of subject {self.SubjectID}. Race must be a string.")
        
    @Sex.setter
    def Sex(self, value):
        # Value has to be a string consisting of anything really.
        if isinstance(value, str):
            self._Sex = value
        else:
            raise ValueError(f"Error found for the sex of subject {self.SubjectID}. Sex must be a string.")

    @Age.setter
    def Age(self, value):
        # Value has to be a float higher than or equal 0
        if isinstance(value, float) and value >= 0:
            self._Age = value
        else:
            raise ValueError(f"Error found for the age of subject {self.SubjectID}. Age must be a float and be higher or equal 0.")

    @BMI.setter
    def BMI(self, value):
        # Value has to be a float higher than 0
        if isinstance(value, float) and value > 0:
            self._BMI = value
        else:
            raise ValueError(f"Error found for the BMI of subject {self.SubjectID}. BMI must be a float and be higher than 0.")

    @SSPG.setter
    def SSPG(self, value):
        # Value has to be a float higher than 0
        if isinstance(value, float) and value >0:
            self._SSPG = value
        else:
            raise ValueError(f"Error found for the SSPG of subject {self.SubjectID}. SSPG must be a float and be higher than 0.")

    @Insulin.setter
    def Insulin(self, value):
        # Value has to be a string of either IR or IS. 
        # Empty values can be described differently, so have to be unified into None before inputting into database
        if isinstance(value, str) and value in ["IR", "IS", None]:
            self._Insulin = value
        else:
            raise ValueError(f"Error found for the insulin value of subject {self.SubjectID}. Insulin must equal IR or IS.")

    """
    Gather_subject(): getter for the attributes not yet had - to source it from the csv file
        working with open ("Subject.csv") file:
            reading the csv file line by line:
                unifying the missing values as "" an empty string,
                extracting the information from consequtive columns separated by "," (split),
                removing any white spaces at the end (rstrip),
                returning the information as subject attributes
    """
    def gather_subject(input):
        line = next(input)          # to skip the header
        for line in input:
            list = ["NA", "Unknown", "unknown", "None", "NULL"]
            line = line.rstrip().split(",")
            new_line = [x if x not in list else None for x in line]
            SubjectID, Race, Sex, Age, BMI, SSPG, Insulin = new_line
            yield Subject(SubjectID, Race, Sex, Age, BMI, SSPG, Insulin)

#Subject_Data = "Subject(1).csv"
#Subjects = SubjectSetter(Subject_Data)      # _input = input = Subject_Data
#for record in Subjects.load_subject():
#    print(record.give_subject())
