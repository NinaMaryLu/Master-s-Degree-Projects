# imports
import sqlite3
import vcf
import logging
import gffutils
import os
import argparse
from Bio.Seq import Seq
import math
import matplotlib.pyplot as plt

# functions
def isfile(file, overwrite):
    # input: file - file to be checked for the existence (default: in the working directory)
    # overwrite - boolean value whether the file should be overwrritten
    if os.path.exists(file) == True and overwrite==False:
        logger.error(f"File {file} already exists. Choose a different prefix and try again.")
        raise SystemExit
    elif os.path.exists(file) == True and overwrite==True:
        logging.warning(f"File {file} already exists and will be overwritten.")
        return False
    else:
        return False

def assertinput(file, filetype):
    # file is a path to the file - here from the command line input
    # filetypes that are checked in that function: "fasta", "database", "gff", "vcf"
    if filetype == "fasta": extension = (".fasta", ".fa") # accepted file extensions as a tuple
    elif filetype == "database": extension = (".db")
    elif filetype == "gff": extension = (".gff", ".gff3")
    elif filetype == "vcf": extension = (".vcf", ".vcf.gz")
    else: logging.warning(f"The filetype input for the file {file} is not supported by the assert function for script's input. It's extension's validity couldn't have been checked.")

    try: 
    # input file exists in that location and file extension is correct
        assert os.path.exists(file)==True and file.endswith(extension)
        logging.info(f"Input file {os.path.basename(file)} is located in the directory: {os.path.join(os.path.abspath(file), os.path.basename(file))}")
    except: #only two cases possible here
        if os.path.exists(file)==False:
            logging.error(f"File {file} wasn't found. Confirm the file path or working directory and try again.")
            raise SystemExit
        else:
            logging.error(f"File {file} has wrong extension for that input. Confirm the placement in the command and try again.")
            raise SystemExit

def assertoutput(file):
    try:
        assert isfile(file, False)==False
    except FileExistsError:
        logging.warning(f"File {file} already exists. Choose a different name and try again.")
        raise SystemExit

def assertprotein(protein):
    # input: protein sequence
    try: 
        assert protein.startswith("M") and protein.endswith("*") and protein.count("*") == 1
    except AssertionError: 
        if protein.startswith("M") == False:
            logging.warning(f"Transcript {parent.id} does not produce a valid protein - it does not start with Methionine.")
        if protein.endswith("*") == False:
            logging.warning(f"Transcript {parent.id} does not produce a valid protein - it does not end with a stop codon.")
        if protein.count("*") != 1:
            logging.warning(f"Transcript {parent.id} does not produce a valid protein - it has more than one stop codon.")

class outputRecord:
    def __init__(self, Chrom ="", Pos ="", Ref ="", Alt ="", Type ="", Transcript = "", ProtLocation ="", RefAA ="", AltAA =""):
        self._Chrom = Chrom
        self._Pos = Pos
        self._Ref = Ref
        self._Alt = Alt
        self._Type = Type
        self._Transcript = Transcript
        self._ProtLocation = ProtLocation
        self._RefAA = RefAA
        self._AltAA = AltAA

    @property
    def Chrom(self):
        return self._Chrom
    
    @property
    def Pos(self):
        return self._Pos
    
    @property
    def Ref(self):
        return self._Ref
    
    @property
    def Alt(self):
        return self._Alt
    
    @property
    def Type(self):
        return self._Type
    
    @property
    def Transcript(self):
        return self._Transcript
    
    @property
    def ProtLocation(self):
        return self._ProtLocation
    
    @property
    def RefAA(self):
        return self._RefAA
    
    @property
    def AltAA(self):
        return self._AltAA
    
    def return_attributes(self):
        attributes_records = [self.Chrom, self.Pos, self.Ref, self.Alt, self.Type, self.Transcript, self.ProtLocation, self.RefAA, self.AltAA]
        return attributes_records
    
    @Chrom.setter
    def Chrom(self, value):
        self._Chrom = value
    
    @Pos.setter
    def Pos(self, value):
        self._Pos = value
    
    @Ref.setter
    def Ref(self, value):
        self._Ref = value

    @Alt.setter
    def Alt(self, value):
        self._Alt = value

    @Type.setter
    def Type(self, value):
        self._Type = value

    @Transcript.setter
    def Transcript(self, value):
        self._Transcript = value

    @ProtLocation.setter
    def ProtLocation(self, value):
        self._ProtLocation = value  

    @RefAA.setter
    def RefAA(self, value):
        self._RefAA = value
    
    @AltAA.setter
    def AltAA(self, value):
        self._AltAA = value

# examplary command line command: 
# python Assignm2.py --Input_vcf testData.vcf.gz --Input_gff PlasmoDB-54_Pfalciparum3D7.gff --Input_fasta PlasmoDB-54_Pfalciparum3D7_Genome.fasta

# python Assignm2.py --Input_vcf assessmentData.vcf.gz --Input_gff PlasmoDB-54_Pfalciparum3D7.gff --Input_fasta PlasmoDB-54_Pfalciparum3D7_Genome.fasta

# logging: logging informations, warning and errors
    # Output: Stream in the command line
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()        # That would put the log into a stream bc StreamHandler. Default is the standard error stream
ch.setFormatter(logging.Formatter("%(levelname)s - %(asctime)s - %(message)s"))
logger.addHandler(ch)

    # Output: Log file documentation, new logs automatically added 
log_output = "Log_error.txt"
fh = logging.FileHandler(log_output)
fh.setFormatter(logging.Formatter("%(levelname)s - %(asctime)s - %(message)s"))
fh.setLevel(logging.INFO)
logger.addHandler(fh)

# arg-parse: arguments for the command line
parser = argparse.ArgumentParser(description='Analyse the VCF file to find unique variations in the protein-coding regions')
parser.add_argument('--Input_vcf', required=True, help='VCF input file path with genomic variations between samples: .vcf file')
parser.add_argument('--Input_gff', required=True, help='GFF input file path with genomic feature data: .gff file')
parser.add_argument('--Input_fasta', required=True, help='Fasta input  path file with genomic sequences: .fasta / .fa file')
args = parser.parse_args()

# checking if the input statements are correct
assertinput(args.Input_vcf, "vcf")
assertinput(args.Input_gff, "gff")
assertinput(args.Input_fasta, "fasta")

# setting output files
output_tsv = "output.tsv"
output_graph = "output.png"

assertoutput(output_tsv)
assertoutput(output_graph)

# Creating a database for the gff file handling
    # creating a database file from gff
DatabaseFile = "output_gff-db.db"
logging.info(f"Database file created from the gff file: {DatabaseFile}")    #change it into the path with os module

if os.path.isfile(DatabaseFile)==True:
    # file exists and we do not want to overwrite it: open the database
    #logger.info(f"Loading the database file {DatabaseFile}...\n")
    try:
        db = gffutils.FeatureDB(dbfn=DatabaseFile, keep_order=True)
    except:
        logger.error(f"Database file couldn't be read from {DatabaseFile}.\n")
        raise SystemExit(1)
else:
    # file exists and we want to overwrite it: create the database
    # file doesn't exist: create the database
    # fn -> source; dbfn -> database file; force=True -> overwrites; keep_order -> attribute order; merge_strategy -> behaviour when PK identical;
    try:
        #logger.info(f"Creating database file of {DatabaseFile}") 
        db = gffutils.create_db(args.Input_gff, dbfn=DatabaseFile, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    except sqlite3.OperationalError:
        logger.error(f"Cannot create database file {DatabaseFile}.\n")

# handling the vcf file
logging.info("Exploring the vcf file and its attributes...\n")
try:
    vcf_reader = vcf.Reader(filename = args.Input_vcf)
except:
    logging.error("vcf file couldn't be opened. Check if the file is actively handled elsewhere.")

logging.info(f"Samples from the data: {vcf_reader.samples}")
count_total = 0             # counting all recods
count_quality = 0           # counting records of quality >20
count_coding = 0            # count records occuring within the coding region
count_noncoding = 0         # count records outside of the coding regions
count_synonymous = 0        # count records that did cause a codon change
count_nonsynonymous = 0     # count records that did not cause a codon change

try:
    with open(output_tsv, "wt") as out:
        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format("Chrom", "Pos", "Ref", "Alt", "Type", "Transcript", "Protein Location", "RefAA", "AltAA"))
        #logging.info("Iterating through the features in the vcf file...")
        # iterate through all the features in the vcf file. Count all features
        for feature in vcf_reader:
            count_total += 1
            # assess how many have a quality of over 20.
            if feature.QUAL >20:
                count_quality += 1
                # checking if the list we are iterating is not empty:
                # the list are all the gff features that are coding, and at the same location as the SNPs (chromosome, start and end)
                if len(list(db.region(seqid=feature.CHROM, start=feature.POS, end=feature.POS, featuretype="CDS"))) > 0:
                    count_coding += 1
                    # iterating through that list
                    for record in db.region(seqid=feature.CHROM, start=feature.POS, end=feature.POS, featuretype="CDS"):
                        for parent in db.parents(record.id, featuretype='mRNA'): # a parent of a CDS that is a mRNA transcript
                            # gathering the necessarily data for the output tsv file
                            f1 = outputRecord()
                            f1.Type = "Coding"
                            f1.Chrom, f1.Pos = feature.CHROM, feature.POS
                            f1.Transcript=parent.id
                            # converting the list to str for consistent output
                            f1.Ref = str(feature.REF[0])
                            f1.Alt = str(feature.ALT[0])
                            # different ordering for the positive and negative strands
                            preceedingCDS_ttl = 0        # sum of the length of the proceeding CDS
                            distance_inCDS = 0           # distance between the change and the start of the affected CDS
                            distanceStart =0             # distance between the change and the start of the transcript
                            distanceEnd=0                # distance between the change and the start of the reverse transcript
                            # for the positive strand
                            try:
                                if parent.strand == "+":
                                    seq=""
                                    for child in db.children(parent.id, featuretype="CDS", order_by="start", reverse=False):
                                        seq = seq + child.sequence(args.Input_fasta, use_strand=True)

                                    for child in db.children(parent.id, featuretype="CDS", order_by="start", reverse=False):
                                        CDS_length = child.stop - child.start + 1
                                        # if SNP is within this CDS
                                        if (f1.Pos >= child.start and f1.Pos <= child.stop) == True:
                                            distance_inCDS = f1.Pos - child.start + 1       # +1 because length calc from positions
                                            distanceStart = preceedingCDS_ttl + f1.Pos - child.start +1
                                            distanceEnd = preceedingCDS_ttl + child.stop - f1.Pos -1
                                            break
                                        else:
                                            preceedingCDS_ttl = preceedingCDS_ttl + CDS_length
                                    transcript_position = distanceStart - 1         # pythonic
                                    #establishing the transcript seq of ref and alt
                                    ref_seq = Seq(seq)      # making a Seq object of the compound sequence
                                    # introducing a base change per vcf file
                                    #print(f"changeg base by ref is: {f1.Ref}, by calculation {ref_seq[transcript_position]}")

                                    try: assert (ref_seq[transcript_position] == f1.Ref)
                                    except AssertionError: logging.warning("Indexing mistake. Expected base is not present at that location.")
                                    alt_seq = Seq((ref_seq[:transcript_position] + f1.Alt + ref_seq[(transcript_position+1):]))
                                    # assertion for ref and alt transcript sequences
                                    #try: assert len(ref_seq) == len(alt_seq)
                                    #except AssertionError: logging.warning("The length of the transcript sequence changed with single base mutation. Impossible.")

                                    ref_protein = Seq.translate(ref_seq)
                                    alt_protein = Seq.translate(alt_seq)

                                    assertprotein(ref_protein)
                                    assertprotein(alt_protein)

                                    f1.ProtLocation = math.ceil(distanceStart / 3) # 1 indexed which codon
                                    f1.RefAA = ref_protein[f1.ProtLocation - 1]    # -1 to make it pythonic 0 indexed
                                    f1.AltAA = alt_protein[f1.ProtLocation - 1]

                                    if f1.RefAA == f1.AltAA:
                                        f1.Type = "Synonymous"
                                        count_synonymous +=1
                                    else:
                                        f1.Type = "Non-Synonymous"
                                        count_nonsynonymous +=1
                                    out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(f1.Chrom, f1.Pos, f1.Ref, f1.Alt, f1.Type, f1.Transcript, f1.ProtLocation, f1.RefAA, f1.AltAA))

                                # for the negative strand
                                # for the negative strand, I have a bug I cannot locate. With assert statements, I see that the base I expect at that position is not 
                                # the same as the one from the VFC file. It works for the testing data, so I suspect the difference to be small. Though, when I check the surrounding 
                                # bases around the index I expect it to be, I see that it varies from feature to feature, and thus have trouble solving it that way. 
                                # I reverse the order of CDS, add the lengths, add the distance between the end of the CDS and position of the change.
                                # I'm leaving some of my assertion statements, prints, comments, in case they are helpful in following my logic.
                                # for the same reason, I'm leaving the whole pipeline from negative strand to the protein (even though some of it is shared between the positive and negative strands)

                                elif parent.strand == "-":
                                    seq=""
                                    for child in db.children(parent.id, featuretype="CDS", order_by="start", reverse=True): # reverse = True for negative strand reversed order
                                        #print(f"CDS seq:\n{child.sequence(args.Input_fasta, use_strand=True)}")
                                        seq = seq + child.sequence(args.Input_fasta, use_strand=True)           # later try with False if everything else fails checked - did not help

                                    for child in db.children(parent.id, featuretype="CDS", order_by="start", reverse=True):
                                        CDS_length = child.stop - child.start + 1

                                        # if SNP is within this CDS
                                        if (f1.Pos >= child.start and f1.Pos <= child.stop) ==True:
                                            distanceStart = preceedingCDS_ttl + f1.Pos - child.start -1
                                            distanceEnd = preceedingCDS_ttl + child.stop - f1.Pos -1
                                            break
                                        else:
                                            preceedingCDS_ttl = preceedingCDS_ttl + CDS_length

                                    transcript_position = distanceEnd -1
                                    #establishing the transcript seq of ref and alt
                                    ref_seq = Seq(seq)
                                    alt_seq = Seq((seq[:transcript_position] + f1.Alt + seq[(transcript_position+1):]))
                                    # first one is exclusive. Then the new base. Then an inclusive index
                                    # assertion for ref and alt transcript sequences
                                    #try: assert len(ref_seq) == len(alt_seq)
                                    #except AssertionError: logging.warning("The length of the transcript sequence changed with single base mutation. Impossible.")
                                    # the position should be calculated from the end
                                    #print(f"changeg base by ref is: {f1.Ref}, by calculation {ref_seq[transcript_position]}")
                                    #print(f"{ref_seq[(transcript_position-5):transcript_position]}--{ref_seq[transcript_position]}--{ref_seq[(transcript_position+1):(transcript_position+6)]}")

                                    try: assert ref_seq[transcript_position] == f1.Ref
                                    except AssertionError: logging.warning(f"Strand {parent.strand} Indexing mistake. Expected base is not present at that location.")                      

                                    ref_protein = Seq.translate(ref_seq)
                                    alt_protein = Seq.translate(alt_seq)

                                    try: 
                                        assert ref_protein.startswith("M") and ref_protein.endswith("*") and ref_protein.count("*") == 1
                                        f1.ProtLocation = math.ceil(distanceStart / 3) # 1 indexed which codon
                                        f1.RefAA = ref_protein[f1.ProtLocation - 1]    # -1 to make it pythonic 0 indexed
                                        f1.AltAA = alt_protein[f1.ProtLocation - 1]
                                        if f1.RefAA == f1.AltAA:
                                            f1.Type = "Synonymous"
                                            count_synonymous +=1
                                        else:
                                            f1.Type = "Non-Synonymous"
                                            count_nonsynonymous +=1
                                        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(f1.Chrom, f1.Pos, f1.Ref, f1.Alt, f1.Type, f1.Transcript, f1.ProtLocation, f1.RefAA, f1.AltAA))
                                    except AssertionError: 
                                        if ref_protein.startswith("M") == False:
                                            logging.warning(f"Transcript {parent.id} does not produce a valid protein - it does not start with Methionine.")
                                        if ref_protein.endswith("*") == False:
                                            logging.warning(f"Transcript {parent.id} does not produce a valid protein - it does not end with a stop codon.")
                                        if ref_protein.count("*") != 1:
                                            logging.warning(f"Transcript {parent.id} does not produce a valid protein - it has more than one stop codon.")
                            except:
                                logging.error(f"Couldn't open {args.Input_fasta}")

                        for parent in db.parents(record.id, featuretype='pseudogene_transcript'):
                            logging.warning(f"{record.id} is not bound to a mRNA, but to a pseudogene and will not be included in the final tsv output")

                else: # The are no matching coding regions for that SNP.
                    count_noncoding += 1
                    if len(list(db.region(seqid=feature.CHROM, start=feature.POS, end=feature.POS, featuretype="mRNA"))) > 0:
                        # in case the lack of protein comes from invalid mRNA or unsuccessful translation == if transcripts exist
                        f1.Transcript = []
                        for record in db.region(seqid=feature.CHROM, start=feature.POS, end=feature.POS):
                            for parent in db.parents(record.id, featuretype='mRNA'):
                                f1 = outputRecord()
                                f1.Type = "Non-Coding"
                                f1.Chrom, f1.Pos = feature.CHROM, feature.POS
                                # converting the list to str for consistent output
                                if len(feature.REF) == 1: f1.Ref = str(feature.REF[0]) # for type inconsistency
                                else: f1.Ref = feature.REF
                                if len(feature.ALT) == 1: f1.Alt = str(feature.ALT[0]) # for type inconsistency
                                else: f1.Alt = feature.ALT
                                f1.Transcript=parent.id                                # has to have a separate line for each
                                f1.ProtLocation = "NA"
                                f1.RefAA = "NA"
                                f1.AltAA = "NA"
                                out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(f1.Chrom, f1.Pos, f1.Ref, f1.Alt, f1.Type, f1.Transcript, f1.ProtLocation, f1.RefAA, f1.AltAA))
                    else: # is there is no CDS and no transcript
                        f1 = outputRecord()
                        f1.Type = "Non-Coding"
                        f1.Chrom, f1.Pos = feature.CHROM, feature.POS
                        # converting the list to str for consistent output
                        if len(feature.REF) == 1: f1.Ref = str(feature.REF[0]) # for type inconsistency
                        else: f1.Ref = feature.REF
                        if len(feature.ALT) == 1: f1.Alt = str(feature.ALT[0]) # for type inconsistency
                        else: f1.Alt = feature.ALT
                        f1.Transcript="NA"
                        f1.ProtLocation = "NA"
                        f1.RefAA = "NA"
                        f1.AltAA = "NA"
                        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(f1.Chrom, f1.Pos, f1.Ref, f1.Alt, f1.Type, f1.Transcript, f1.ProtLocation, f1.RefAA, f1.AltAA))
        logging.info(f"Out of {count_total} features, {count_quality} have quality score of more than 20. {count_total - count_quality} features had a subpar quality.")

        # producing a bar plot
        names = ["Non-coding", "Synonymous", "Non-synonymous"]
        values = [count_noncoding, count_synonymous, count_nonsynonymous]   # I realise those are incorrect due to the negative strand issues described above
        # I understood the instructions as meaning proportion = split = how many are synonymous, non coding etc. 
        # If I was to calculate the actual proportion, I'd be dividing those values by count_quality
        plt.bar(names, values, color="#6DB2FC")
        plt.title('Types of variations in the file')
        plt.xlabel("Type")
        plt.ylabel('# of variations')
        plt.show()
        plt.savefig(output_graph, dpi='figure', format="png")

        # produced at the end - only if the script runs correctly.
        logging.info(f"Tab-separated output file was saved as {output_tsv} (working directory)")
        logging.info(f"Graph output file was saved as {output_graph} (working directory")
        out.close()
except FileExistsError:
    logging.error(f"File {output_tsv} couldn't be created.")


# producing the bar graph with fraction instead of absolute numbers
# total = count_quality ie. number of variants with quality of over 20.
#fraction_synonymous = count_synonymous / count_quality *100
#fraction_nonsynonymous = count_nonsynonymous / count_quality *100
#fraction_noncoding = count_noncoding / count_quality *100




