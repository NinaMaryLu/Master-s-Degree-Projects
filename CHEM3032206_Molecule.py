class Molecule:
    # the constructor for the object
    def __init__(self, CdId = "", Mol_ID = "", Image = "", SMILES="", Formula = "", Name="", 
                MolWeight = "", LogP = "", LogD ="", HydroDonor="", HydroAccept = "", RotatableBonds = "", NumberRings="",
                Lipinski ="", leadlikeness = "", bioavailability = ""):
        # To avoid unintentional replacing of data in the database, copies of attributes are used in all data manipulations
        self._CdId = CdId
        self._Mol_ID = Mol_ID
        self._Image = Image
        self._SMILES = SMILES
        self._Formula = Formula
        self._Name = Name
        self._MolWeight = MolWeight
        self._LogP = LogP
        self._LogD = LogD
        self._HydroDonor = HydroDonor
        self._HydroAccept = HydroAccept
        self._RotatableBonds = RotatableBonds
        self._NumberRings = NumberRings
        self._Lipinski = Lipinski
        self._leadlikeness = leadlikeness
        self._bioavailability = bioavailability

    @property
    def CdId(self):
        return self._CdId
    
    @property
    def Mol_ID(self):
        return self._Mol_ID
    
    @property
    def SMILES(self):
        return self._SMILES
    
    @property
    def Image(self):
        return self._Image
    
    @property
    def Formula(self):
        return self._Formula
    
    @property
    def Name(self):
        return self._Name

    @property
    def MolWeight(self):
        return self._MolWeight

    @property
    def LogP(self):
        return self._LogP

    @property
    def LogD(self):
        return self._LogD
    
    @property
    def HydroDonor(self):
        return self._HydroDonor
    
    @property
    def HydroAccept(self):
        return self._HydroAccept
    
    @property
    def RotatableBonds(self):
        return self._RotatableBonds
    
    @property
    def NumberRings(self):
        return self._NumberRings
    
    @property
    def Lipinski(self):
        return self._Lipinski
    
    @property
    def leadlikeness(self):
        return self._leadlikeness
    
    @property
    def bioavailability(self):
        return self._bioavailability
    
    # function for returning the values for the subject as a list of values
    def give_attributes(self):
        row_values = [self.CdId, self.Mol_ID, self.Image, self.SMILES, self.Formula, self.Name, self.MolWeight,
                    self.LogP, self.LogD, self.HydroDonor, self.HydroAccept, self.RotatableBonds, self.NumberRings, self.Lipinski, self.leadlikeness, self.bioavailability]
        return row_values
    
    def validate_Lipinski(self):
        Lipinski_check = 0
        if self.HydroDonor <= 5:
            Lipinski_check += 1
        
        if self.HydroAccept <=10:
            Lipinski_check +=1
        
        if self.MolWeight <=500:
            Lipinski_check +=1
        
        if self.LogP <=5:
            Lipinski_check +=1
        
        if Lipinski_check >=3:
            Lipinski_result = True
        else:
            Lipinski_result = False
        return Lipinski_result

    def validate_LeadLikeness(self):
        LeadLikeness_check = 0
        if self.MolWeight <= 450:
            LeadLikeness_check += 1
        
        if (self.LogD >= -4) and (self.LogD <=4):
            LeadLikeness_check +=1
        
        if self.NumberRings <=4:
            LeadLikeness_check +=1
        
        if self.RotatableBonds <=10:
            LeadLikeness_check +=1
        
        if self.HydroDonor <= 5:
            LeadLikeness_check +=1
        
        if self.HydroAccept <= 8:
            LeadLikeness_check +=1
        
        if LeadLikeness_check == 6:
            LeadLikeness_result = True
        else:
            LeadLikeness_result = False
        return LeadLikeness_result

    #def validate_bioavailability(self):




    """
    SETTERS for the object attributes
    To validate if the values match certain criteria
    """
    @CdId.setter
    # Value has to be an int / str?. Has to be unique. Can't be empty.
    def CdId(self, value):
        self._CdId = value

    @Mol_ID.setter
    # Value has to be an int / str?. Has to be unique. Can't be empty.
    def Mol_ID(self, value):
        self._Mol_ID = value

    @Image.setter
    # Value has to be an int / str?. Has to be unique. Can't be empty.
    def Image(self, value):
        self._Image = value
    
    @SMILES.setter
    def SMILES(self, value):
        # Value has to be a string consisting of anything really. Has to be unique. Can't be empty.
        if isinstance(value, str) and value is not None:
            self._SMILES = value
        else:
            raise ValueError(f"Error found for the SMILES value of {self.Name}:{self.SMILES}. SMILES must be a string and can't be empty.")

    @Formula.setter
    # Value has to be an int / str?. Has to be unique. Can't be empty.
    def Formula(self, value):
        self._Formula = value

    @Name.setter
    def Name(self, value):
        # Value has to be a string consisting of anything really.
        if isinstance(value, str):
            self._Name = value
        else:
            raise ValueError(f"Error found for the name of {self.Name}. Name must be a string.")

    @MolWeight.setter
    def MolWeight(self, value):
        # Value has to be a float
        if isinstance(value, float):
            self._MolWeight = value
        else:
            raise ValueError(f"Error found for the molecular weight of {self.Name}: {self.MolWeight}. Molecular weight must be a float.")
        
    @LogP.setter
    def LogP(self, value):
        # Value has to be a float
        if isinstance(value, float):
            self._LogP = value
        else:
            raise ValueError(f"Error found for the LogP of {self.Name}: {self.LogP}. LogP must be a float.")
        self._LogP = value

    @LogD.setter
    # Value has to be  a float
    def LogD(self, value):
        # Value has to be a float
        if isinstance(value, float):
            self._LogD = value
        else:
            raise ValueError(f"Error found for the LogD of {self.Name}: {self.LogD}. LogD must be a float.")
        self._LogD = value

    @HydroDonor.setter
    def HydroDonor(self, value):
        # Value has to be an int
        if isinstance(value, int):
            self._HydroDonor = value
        else:
            raise ValueError(f"Error found for the number of hydrogen bond donors of {self.Name}: {self.HydroDonor}. HydroDonor must be an integer.")
        self._HydroDonor = value

    @HydroAccept.setter
    # Value has to be an int
    def HydroAccept(self, value):
        if isinstance(value, int):
            self._HydroAccept = value
        else:
            raise ValueError(f"Error found for the hydrogen bond acceptors of {self.Name}: {self.HydroAccept}. HydroAccept must be an integer.")
        self._HydroAccept = value

    @RotatableBonds.setter
    # Value has to be an int
    def RotatableBonds(self, value):
        if isinstance(value, int):
            self._RotatableBonds = value
        else:
            raise ValueError(f"Error found for the number of rotatable bonds {self.Name}: {self.RotatableBonds}. RotatableBonds must be an integer.")
        self._RotatableBonds = value

    @NumberRings.setter
    # Value has to be an int
    def NumberRings(self, value):
        if isinstance(value, int):
            self._NumberRings = value
        else:
            raise ValueError(f"Error found for the number of rings {self.Name}: {self.NumberRings}. RotatableBonds must be an integer.")
        self._NumberRings = value

    @Lipinski.setter
    # Value has to be a Bolean or none
    def Lipinski(self, value):
        self._Lipinski = value

    @leadlikeness.setter
    # Value has to be a Bolean or none
    def leadlikeness(self, value):
        self._leadlikeness = value

    @bioavailability.setter
    # Value has to be a Bolean or none
    def bioavailability(self, value):
        self._bioavailability = value

