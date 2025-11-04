import re
import argparse
import logging

def CIGARintoLIST(CIGAR):                   # separating the CIGAR into list of (bases + type) elements
    CIGARoptions = ["M","I","D","N","S"]
    for character in CIGAR:
        if character in CIGARoptions:
            newCharacter = character + ","      
            CIGAR = CIGAR.replace(character,newCharacter)
            # if there is more than one instance of the same character, one of them is going to have more than one instance of ","
    CIGAR = CIGAR.split(",")
    # the "" produced by multiple "," are now removed from the list
    while "" in CIGAR:          
        CIGAR.remove("")
    return CIGAR

# Creating a parser for the use of positional arguments from the command line for the script
parser = argparse.ArgumentParser(description= "Arguments needed for the script:")
parser.add_argument("reads", nargs='?', type=str, help = "Input file for the script: sam file with reads", default = "None")             # without "--" for positional arguments # I assume that handles missing file
parser.add_argument("genes", nargs='?', type=str, help = "Input file for the script: txt file with gene locations", default = "None")       # I assume that handles missing file
args = parser.parse_args()

# creating a logger for hanfling errors nicely
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter("%(levelname)s - %(asctime)s - %(message)s"))
logger.addHandler(ch)

# Creating variables for the input files. The script takes 2 input files: 1) .sam file 2) .txt file 
input_sam = args.reads            #"Toxo_chr8_subset.sam"            # path of the file will be given in a command line
input_table = args.genes          #"Toxo_Gene_Locations.txt"          # path of the file will be given in a command line

"""
1. READING THE SAM FILE
      RNAME : 3rd column      chromosome that it alignes to
      POS : 4th column        mapping position - start position of the alignment
      CIGAR : 6th column      describes the alignment
      NH:i:x : last column    how many times the read aligned
"""
RNAME_intron = []       # list of tuples (strings + int) to store position name (RNAME) and positions of introns
try:    
    with open(input_sam) as my_file:
        for line in my_file:
            # reads every row of the txt file as a separate line, cuts out the finishing "\n" and separates the columns based on tabs   
            line = line.rstrip().split("\t")       

            # skip all the header rows that start with @. 
            if line[0].startswith("@"):            
                pass

            # Read Name, POS, CIGAR and NH:i:x from columns 3, 4, 6 and the last, respectively.
            else:
                RNAME = line[2]                                 # 3rd column in homo sapiens terms
                POS = int(line[3])                              # 4th column in homo sapiens terms
                CIGAR = line[5]                                 # 6th column in homo sapiens terms
                NH = line[-1]                                   # last column
                NH_index = NH.split(":")                        # separates NH into variables to extract the number of alignments

                # check if there read is only aligned once. If more than once or none, ignore
                if NH_index[2] == "1":

                    # check if the read is split by any introns
                    if CIGAR.count("N") != 0:
                        CIGAR_list = CIGARintoLIST(CIGAR)       # separating the CIGAR into list of (bases + type) elements

                        # Creats two lists with matching order - for the type of the read (match, deletion etc.) and its length.
                        CIGAR_list_TYPE = []
                        CIGAR_list_LENGTH = []
                        for item in CIGAR_list: 
                            CIGAR_list_TYPE.append(item[-1])
                            CIGAR_list_LENGTH.append(int(item[:-1]))

                        # creating a tuple of CIGAR readings (type + length)
                        # based on: https://www.geeksforgeeks.org/python-create-a-list-of-tuples/
                        CIGAR_TYPE_LENGTH = []
                        for i in range(len(CIGAR_list_TYPE)):
                                CIGAR_TYPE_LENGTH.append((CIGAR_list_TYPE[i], CIGAR_list_LENGTH[i]))

                        sum = 0
                        # detecting the position of an intron
                        for object in CIGAR_TYPE_LENGTH:
                            if object[0] == "N":
                                N_position = CIGAR_TYPE_LENGTH.index(object)
                                N_number = CIGAR_TYPE_LENGTH[N_position][1]

                                # calculating the number of bases preceeding the intron.
                                for i in range(N_position):
                                    if (CIGAR_TYPE_LENGTH[i][0] == "M" or CIGAR_TYPE_LENGTH[i][0] == "N"):      # only counting matches and preceeding junctions
                                        sum = sum + CIGAR_TYPE_LENGTH[i][1]

                                # Position of the intron on the chromosome based on POS and CIGAR
                                N_start = int(POS) + sum
                                N_end = N_start + N_number

                                # Adding to a tuple of (gene name, intron start, intron end)
                                RNAME_intron.append((RNAME, N_start, N_end))
                    else:    # when there are no intron sequences
                        pass
                else:    # ignores all reads with more than 1 match
                    pass
except FileNotFoundError:
    if input_sam == "None":
        # I wanted to make an exception loop as below but I can't make it work. I could assume the extensions of the files are ".txt" and ".sam" but I can't know that. 
        # otherwise as long as the table is present, the default value here will not be read, but the input_table.
        logger.error(f"File name for the input file with readings was not given (positional argument 1). Please input a file name after the name of the python script.")
    else:        
        logger.error(f"File {input_sam} not found. Please check whether the file placement and the file name given are correct")
    raise SystemExit
                
# finding all the unique introns = unique values in the tuple (same name, start and end of the intron)
introns = set(RNAME_intron)     # using a set to gather unique junctions (same gene name, start and stop)
introns = list(introns)         # converting to list for further use

# adding a 4th position to the tuple - counting the number a specific reading exists on the list
RNAME_intron_reads = []  # list of tuples (strings + int) to store position name (RNAME), positions of introns and number of supporting reads.
for intron in introns:
    intron_index = RNAME_intron.index(intron)
    read_count = RNAME_intron.count(intron)
    # tuple can't be amended - creating an alternative tuple: (gene name, start, stop, read count)
    RNAME_intron_reads.append((RNAME_intron[intron_index][0], RNAME_intron[intron_index][1], RNAME_intron[intron_index][2], read_count)) 

"""
2. READING THE TAB-SEPARATED FILE
    has a header + 3 columns: geneID, transcriptID, GeneLocation
    Genome location read as redex?
        chromosome name: "{name}_chr{RomanNumeral}"
        start position: number separated with "," for every thousand
        ".."
        end position: number separated with "," for every thousand
        strand: "(-)"
"""
try:
    with open(input_table) as my_file2:
        with open("3032206.txt", "w") as out:
            geneMAP = {}
            next(my_file2)        # skips the header
            for line in my_file2:
                # reads every row of the txt file as a separate line, cuts out the finishing "\n" and separates the columns based on tabs   
                line = line.rstrip().split("\t")
                GeneID, sourceID, GenomLoc = line       # sourceID not accessed
                GenomLoc = GenomLoc.replace(",","")     # leaves numbers uninterrupted by commas
                geneMAP[GeneID] = GenomLoc

                for GeneID, GenomLoc in geneMAP.items():
                    for GeneID in geneMAP:
                        # using redex to read the Gene Location to extract gene name, its start and end
                        match = re.search(r"(^.+_.+)(:)(\d+)([^\d]{2})(\d+)([^\d]{3})", GenomLoc)
                        # Group1: 
                        #   Gene name: starting with any character repeated at least once
                        #   Chromosome number: chr + roman numeral. Longest Roman numeral below 1440 has 14 characters, but let's assume any. 
                        #   After all, this does not have to be a chromosome number, but has to match.
                        # Group2: separator
                        # Group3: start of the gene: any digit, repeated at least one (includes position 0 for simplicity, but I admit that's impossible)
                        # Group3: separator ".." but can be any character that is not a digit (not stated in instructions as "..")
                        # Group4: end of the gene: any digit, repeated at least one

                        if match:
                            geneName = match.groups()[0]         # Group1: Name_chromosome
                            geneStart = int(match.groups()[2])   # Group3: Gene start position
                            geneStop = int(match.groups()[4])    # Group4: Gene stop position
                            # bring both together: read the junction info and find all that are within genes
                            for junction in RNAME_intron_reads:
                                junction_index = RNAME_intron_reads.index(junction)  # I need the index for further iteration
                                # if the names match...
                                if geneName == RNAME_intron_reads[junction_index][0]:
                                        # if the junction starts within the gene. As I understand it can logically either end within it, or exceed it. 
                                        # But if it exceeds it, then the gene ends when the intron starts, so such situation can't exist. 
                                    if ((RNAME_intron_reads[junction_index][1] > geneStart) and (RNAME_intron_reads[junction_index][1] < geneStop)):
                                        out.write(f'{GeneID}\t{RNAME_intron_reads[junction_index][1]}\t{RNAME_intron_reads[junction_index][2]}\t{RNAME_intron_reads[junction_index][3]}\n')
                                        # (f'{GeneID}\t{Junction start}\t{Junction end}\t{o	Number of reads supporting the junction}\n')
                                    else:
                                        pass
                                else:
                                    pass
                        out.write(f"\n")        # writes an empty line in between genes
except FileNotFoundError: 
    if input_table == "None":
        logger.error(f"File name for the input file with gene locations was not given (positional argument 2). Please input a file name after the name of the .sam file.")
    else:        
        logger.error(f"File {input_table} not found. Please check whether the file placement and the file name given are correct")
    raise SystemExit