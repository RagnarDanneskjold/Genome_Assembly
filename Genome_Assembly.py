#!/usr/bin/python
import sys, os, re, getopt, glob, subprocess, os.path, numpy as np
import timeit
import itertools

now = timeit.default_timer()

Usage = """
Usage:  ./Genome_Assembly.py -o ASSEMBLED -i READ_FILES_LIST.tsv -A TruSeq3 -F Y -M 2,4

REQUIRED ARGUMENTS:
		-o	output directory
	
                -i	your input file names in the form of a .tsv file with four columns : "Sample", "Paired End File 1" and "Paired End File 2 and Phylogenetice Gene Marker Sequence"
			(accepted compression formats: .gz & .bz2 - NOTE: compressed formats must have the same name as read file + ".gz" or ".bz2")

DEPENDENCIES:
		-FastX Toolkit (built with v.0.7)
		-khmer (specifically the commands: interleave-reads.py and extract-paired-reads.py)
		-KmerGenie
		-MeGAMerge (built with v1.1)
		-Any of : 	-Velvet (built with v.1.2.02)
				-IDBA (built with v.1.1.1)
				-RAY (built with v.v2.3.1)
				-ABySS (built with v1.3.7)
				## NOTE: IDBA does NOT make use of unpaired reads
OPTIONAL SOFTWARE:
		-Trimmomatic (built with v.0.32) 
		-FLASH (built with v.1.2.10)
		-EA-UTILS	
OPTIONAL:
		## Related To Adapter Removal
		-A 	Specify whether to Remove Illumina Adapters using Trimmomatic
			(Must specify which adapter library to use, i.e. "NexteraPE","TruSeq2" or "TruSeq3")

		-E <Y>	Specify to use EA-Utils for Quality Cleaning instead of FastQC
			(Included this as a legacy from previous genome assembling)

		## Related to Quality Trimming
		-Q	Specify Phred Format  [default: Phred33]
			(Usage -Q 64)

		-q	Specify Minimum Quality Score for Quality Trimming [default: 30]
			(Usage -q 25)

		-l	Specify Minimum Read Length to Keep With Above Quality Score [default: 50]
			(Usage -l 75)

		## Related to Improving Assembly
		-F <Y>	Use FLASH to extend paired end reads

		## Related to Assembly
		-M	Specify which assembler to use by either name or number. One can specify multiple assemblers, which will yield separate folders with assembled data in the /ASSEMBLY/ folder.
			(1: IDBA, 2: Velvet, 3: RAY, 4: ABySS)

		-C	Specify whether you would like to input custom parameters for assembly. NECESSARY if one wishes to perform MULTIPLE ITERATIONS of k-mer length
			[default: k-mer length is selected with software tool "KmerGenie"]
			(This will PAUSE the script before assembly and lead the user through whatever changes they would like to make, INCLUDING the use of additional parameters)

		-m 	Perform super assembly with MeGAMerge
			(NOTE: If selected, the user will be prompted to select a range of multiple k-mer lengths for iteration)

UTILITY:
This script will perform the basic quality control (i.e. remove artefacts from Illumina adapters and quality trim) and assemble paired-end genomic data with a number of options.

NOTES:
- File extensions accurately denote what kind of analysis has been done to the file:

.pe - Interleaved Paired end reads file
.se - Unpaired reads file
.fq - Fastq format
.qc - Quality trimmed
.tr - Adapter trimmed
.fl - Flash merged

- Logging is set ON to by default. Editorializing is logged to "run.log", while all commands performed on your files in "command.log". This can be shut off by changing the section of script "LOG = 'ON'" to 'OFF'


Usage:  ./Genome_Assembly.py -o ASSEMBLED -i READ_FILES_LIST.tsv -A TruSeq3 -F Y -M 2,4

"""

if len(sys.argv)<2:
        print Usage
        sys.exit()

# Initialize all containers for input 
LIST=''
ADAPTER=''
QUAL_FORMAT=''
MIN_QUAL=''
MIN_LENGTH=''
FLASH=''
ASSEMBLER=''
MERGE=''
CUSTOM=''
OUTPUT=''
ITERATE=''
EA_UTILS=''

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"o:i:A:Q:q:l:L:F:D:M:s:C:m:E:")

###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-o':
        OUTPUT= a
    if o == '-i':
        LIST= a
    if o == '-A':
        ADAPTER= a
    if o == '-Q':
        QUAL_FORMAT= a
    if o == '-q':
        MIN_QUAL= a
    if o == '-l':
        MIN_LENGTH= a
    if o == '-F':
        FLASH= a
    if o == '-S':
        PROCESSES= a
    if o == '-M':
        ASSEMBLER= a
    if o == '-m':
        MERGE= a
    if o == '-C':
        CUSTOM= a
    if o == '-E':
        EA_UTILS= a

LOG = 'ON'

if LOG == 'ON':
        print "\n-- You are Logging a Description of Analyses and All Commands --"
	log = open("run.log", "w")
	command = open("command.log", "w")

if len(OUTPUT)>0:
        if os.path.exists('./' + OUTPUT):
                print "\n-- Output Folder Exists - Caution: Files May Be Over-written --"
        else:
                os.mkdir(OUTPUT)

####################################
# IMPORT FILE LIST AND HOUSE KEEPING
####################################

## IMPORT FILE NAMES into DICTIONARY
FILE_DICT={}
SEQUENCE_DICT={}

with open(LIST) as f:
	next(f)

	for line in f:
		line = line.strip("\r\n")
		line = line.split("\t")
		FILE_DICT[line[0]] = [line[1], line[2]]
		SEQUENCE_DICT[line[0]] = line[3]

if LOG == 'ON':
	log.write("Your Input Files Were:\n\n")
	log.write(str(FILE_DICT))
	command.write("Your Input Files Were:\n\n")
	command.write(str(FILE_DICT))

## Set Defaults unless user specified otherwise
if len(QUAL_FORMAT) == 0:
	QUAL_FORMAT = "33"

if len(MIN_QUAL) == 0:
	MIN_QUAL = "30"

if len(MIN_LENGTH) == 0:
	MIN_LENGTH = "50"

if LOG == 'ON':
	log.write("\n\nThe following defaults were used:\n\nFastX Quality Trim\nPhred Format: "+QUAL_FORMAT+"\nMinimum Quality Score: "+MIN_QUAL+"\nMinimum Length: "+MIN_LENGTH)

## Initialize Output Dictionary To Store Paired and Unpaired Read File Names Output at Different Stages of the Pipeline 
OUTPUT_DICT={}

## Define Function for Unzipping and Interleaving Virgin Paired Reads At Various Entry Points for User along the Pipeline
def UnZ_InteR():
        #Unzip if necessary
        if re.search(".bz2", value[0]):
                os.system(' '.join([
                        "tar -xjfv",
                        value[0],
                        value[1]
                ]))

                # Remove zip extension from name values
                value[0] = re.sub(".bz2","", value[0])  #"value" corresponds to the name of the the paired ends file names stipulated in the inputted .tsv file and is named as such b/c it is called from a dictionary
                value[1] = re.sub(".bz2","", value[1])

        if re.search(".gz", value[0]):
                os.system(' '.join([
                        "gunzip",
                        value[0],
                        value[1]
                ]))

                # Remove zip extension from name values
                value[0] = re.sub(".gz","", value[0])
                value[1] = re.sub(".gz","", value[1])

        if len(value[1]) > 0:
                os.system(' '.join([
                        "interleave-reads.py ./"+value[0],
                        "./"+value[1],
                        ">",
                        "./"+OUTPUT+"/"+key+".pe.fq"
                        ]))

                if LOG == 'ON':
                        command.write(' '.join([
                                "\ninterleave-reads.py ./"+value[0],
                                "./"+value[1],
                                ">",
                                "./"+OUTPUT+"/"+key+".pe.fq"
                        ]))
        else:
                log.write("\nWe Only Found a Single Reads File and Therefore Assumed the Provided Reads Were Already Interleavened.\n")
                os.system(' '.join([
                        "cp ./"+value[0],
                        "./"+OUTPUT+"/"+key+".pe.fq"
                ]))

        value[0] = "./"+OUTPUT+"/"+key+".pe.fq" #Attribute 1st naming value to PE
        value[1] = ''                           #Attribute 2nd naming value to SE

## Define Function to Strip ".gz" from file names
def Clean_GZ(x):
	#Remove ".gz" file extension from name
	x = re.sub(".gz", "", x)
	return x

def Convert_FQ(x):
	if re.search(".fq", x):
		os.system(' '.join([
			"fastq_to_fasta -n -i",
			x,
			"-o",
			re.sub("fq","fa",x)
		]))
	
		if LOG == 'ON':
			command.write(' '.join([
				"\nfastq_to_fasta -n -i",
                     		x,
	                        "-o",
        	                re.sub("fq","fa",x)
                	]))
			
			log.write("We converted your .fastq files into .fasta format")
			
		x = re.sub("fq","fa",x)

		return x

## Small Function To Re-name ".pe" and ".se" read files for RAY-meta assembly
def EXTENSION_STRIP(READS,EXT):

        os.system(' '.join([
                "cp",
                READS,
                re.sub(".fa"+EXT, EXT+".fa",READS)
        ]))

        READS = re.sub(".fa"+EXT, EXT+".fa",READS)

        return READS
		
def FASTQC(IN,OUT):
	os.system(' '.join([
		"fastq_quality_filter",
		"-Q"+QUAL_FORMAT,
		"-q"+MIN_QUAL,
		"-p"+MIN_LENGTH,
		"-i",
		IN,
		">",
		OUT
	]))

	if LOG == 'ON':
		command.write(' '.join([
                       	"\nfastq_quality_filter",
      	                "-Q"+QUAL_FORMAT,
                        "-q"+MIN_QUAL,
       	       	        "-p"+MIN_LENGTH,
			"-i",
			IN,
			">",
			OUT
 	   	]))

def EAUTILS(IN,OUT):

	dummy = open("adapter-file.fq", "w")
	dummy.close()

	os.system(' '.join([
		"fastq-mcf",
                "adapter-file.fq",
                IN,
                "-o",
		OUT,
                ">>",
                "./" + OUTPUT + "/" + key + "_fastq-mcf.log",
                "2>&1"
	]))

	if LOG == 'ON':
		command.write(' '.join([
			"fastq-mcf",
        	        "adapter-file.fq",
                	IN,
	                "-o",
			OUT,
                	">>",
	                "./" + OUTPUT + "/" + key + "_fastq-mcf.log",
        	        "2>&1"
 	   	]))

	os.system('rm ./adapter-file.fq')

#########################
##QUALITY FILTERING STEPS
#########################

print "\n-- Performing Quality Trimming --\n\n"
	
## Do TRIMMOMATIC If Specified
if ADAPTER:
	print "-- You Have Selected To Trim Illumina Adapters Using the Following Adapter File: "+ADAPTER
	if LOG == 'ON':
		log.write("\nYou performed adapter trimming with"+ADAPTER+"as the adapter type.\n")

	#Make Temporary Trim Folder
       	if os.path.exists("TRIM"):
		pass
	else:
                os.mkdir("TRIM")

	for key, value in FILE_DICT.iteritems():
		
		if LOG == 'ON':
			log.write("\n\nYou are now processing reads associated with: "+key+".\n")
			command.write("\n\nYou are now processing reads associated with: "+key+".\n")

		## Perform Adapter Removal
		os.system(' '.join([
			"java -jar /usr/local/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE",
			value[0],
			value[1]+" ./TRIM/s1_pe ./TRIM/s1_se ./TRIM/s2_pe ./TRIM/s2_se",
			"ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.32/adapters/"+ADAPTER+"-PE.fa:2:30:10"
		]))

		if LOG == 'ON':
			command.write(' '.join([
				"\njava -jar /usr/local/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE",
              	                value[0],
                        	value[1]+"./TRIM/s1_pe ./TRIM/s1_se ./TRIM/s2_pe ./TRIM/s2_se",
       	                        "ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.32/adapters/"+ADAPTER+"-PE.fa:2:30:10"
        	       	]))

		## Interleave Reads from the two paired ends files
		os.system(' '.join([
			"interleave-reads.py ./TRIM/s?_pe",
			"> ./TRIM/combined.fq"
		]))

		if LOG == 'ON':
			command.write(' '.join([
               	                "\ninterleave-reads.py ./TRIM/s?_pe",
       	                        "> ./TRIM/combined.fq"
                       	]))

		############
		#QC Trimming
		############

		if EA_UTILS:
			EAUTILS("./TRIM/combined.fq","./TRIM/combined-trim.fq")
			EAUTILS("./TRIM/s1_se","./TRIM/s1_se.trim")
			EAUTILS("./TRIM/s2_se","./TRIM/s2_se.trim")

		#Use FastQC
		else:
			FASTQC("./TRIM/combined.fq","./TRIM/combined-trim.fq")
			FASTQC("./TRIM/s1_se","./TRIM/s1_se.trim")
			FASTQC("./TRIM/s2_se","./TRIM/s2_se.trim")
	
		## Extract ONLY paired reads from Quality Trimmed File (Creating a third orphaned single pair file)
		os.system(' '.join([
			"extract-paired-reads.py ./TRIM/combined-trim.fq"
		]))

		if LOG == 'ON':
			command.write("\nextract-paired-reads.py ./TRIM/combined-trim.fq")

		## Merge Paired Reads using FLASH if user specified
		if FLASH:
			if LOG == 'ON':
				log.write("\nYou Selected To Merge Paired Reads Using FLASH.\n")

			os.system(' '.join([
				"flash --interleaved -z combined-trim.fq.pe",
			]))

			if LOG == 'ON':
				command.write("\nflash --interleaved -z combined-trim.fq.pe")

			## Merge three other orphan single read files with newly extended single reads from FLASH
			# First create .gz of all single read files up to this point	
			os.system(' '.join([
				"gzip -9c combined-trim.fq.se",
				"./TRIM/s1_se.trim",
				"./TRIM/s2_se.trim",
				"> ./combined-trim.fq.se.gz"
			]))

			os.system(' '.join([
				"cat out.extendedFrags.fastq.gz combined-trim.fq.se.gz",
				"> ./"+OUTPUT+"/"+key+".se.tr.qc.fl.fq.gz"
			]))

			os.system(' '.join([
				"mv out.notCombined.fastq.gz",
				"./"+OUTPUT+"/"+key+".pe.tr.qc.fl.fq.gz"
			]))
		
			## Save output file names
			OUTPUT_DICT[key] = ["./"+OUTPUT+"/"+key+".pe.tr.qc.fl.fq.gz","./"+OUTPUT+"/"+key+".se.tr.qc.fl.fq.gz"]

		else:
			## Re-Compress Files and Move to Output Folder
			os.system(' '.join([
				"gzip -9c combined-trim.fq.pe > combined-trim.fq.pe.gz",
			]))

			os.system(' '.join([
				"mv combined-trim.fq.pe.gz",
				"./"+OUTPUT+"/"+key+".pe.tr.qc.fq.gz"
			]))

			os.system(' '.join([
				"gzip -9c combined-trim.fq.se",
				"./TRIM/s1_se.trim",
				"./TRIM/s2_se.trim",
				" > combined-trim.fq.se.gz"
			]))

			os.system(' '.join([
				"mv combined-trim.fq.se.gz",
				"./"+OUTPUT+"/"+key+".se.tr.qc.fq.gz"
			]))

			## Save output file names
             		OUTPUT_DICT[key] = ["./"+OUTPUT+"/"+key+".pe.tr.qc.fq.gz","./"+OUTPUT+"/"+key+".se.tr.qc.fq.gz"]


	## Remove Temp Folder "TRIM" and clean up
	os.system(' '.join([
		"rm -fr TRIM",
	]))

	os.system(' '.join([
		"rm combined-trim.fq.pe combined-trim.fq.se",
	]))

	if FLASH:
		os.system(' '.join([
			"rm out.extended* out.*hist* combined-trim.fq.se.gz",
		]))

## Do Not Perform Adapter trimming
else:
	if LOG == 'ON':
		log.write("\n\n-- You did NOT perform adapter trimming. --\n")

	#Make Temporary Trim Folder
        if os.path.exists("TRIM"):
		pass
	else:
        	os.mkdir("TRIM")

	for key, value in FILE_DICT.iteritems():

		if LOG == 'ON':
			log.write("\n\nYou are now processing reads associated with: "+key+".\n")
			command.write("\n\nYou are now processing reads associated with: "+key+".\n")

		#Pretreatment of virgin paired end reads
		UnZ_InteR()

		os.system(' '.join([
			"cp",
			value[0],
			"./TRIM/combined.fq"
		]))

		if LOG == 'ON':
			command.write(' '.join([
				"cp",
                               	value[0],
	       	                "./TRIM/combined.fq"
               		]))

		############
		#QC Trimming
		############

		if EA_UTILS:
			EAUTILS("./TRIM/combined.fq","./TRIM/combined-trim.fq")
			EAUTILS("./TRIM/s1_se","./TRIM/s1_se.trim")
			EAUTILS("./TRIM/s2_se","./TRIM/s2_se.trim")

		#Use FastQC
		else:
			FASTQC("./TRIM/combined.fq","./TRIM/combined-trim.fq")
			FASTQC("./TRIM/s1_se","./TRIM/s1_se.trim")
			FASTQC("./TRIM/s2_se","./TRIM/s2_se.trim")

		## Extract ONLY paired reads from Quality Trimmed File
		os.system(' '.join([
			"extract-paired-reads.py ./TRIM/combined-trim.fq"
		]))

		if LOG == 'ON':
			command.write(' '.join([
	                                "\nextract-paired-reads.py ./TRIM/combined-trim.fq"
               	        ]))

		## Merge Paired Reads if users specified
		if FLASH:
			if LOG == 'ON':
				log.write("\nYou Selected To Merge Paired Reads Using FLASH.\n")
	
			os.system(' '.join([
				"flash --interleaved -z combined-trim.fq.pe",
			]))

			if LOG == 'ON':
				command.write("\nflash --interleaved -z combined-trim.fq.pe")

			## Merge other single reads with newly merged paired end reads
			# First create .gz of all single read files up to this point	
			os.system(' '.join([
				"gzip -9c combined-trim.fq.se",
				"> ./combined-trim.fq.se.gz"
			]))
			
			os.system(' '.join([
				"cat out.extendedFrags.fastq.gz combined-trim.fq.se.gz",
				"> ./"+OUTPUT+"/"+key+".se.qc.fl.fq.gz"
			]))

			os.system(' '.join([
				"mv out.notCombined.fastq.gz",
				"./"+OUTPUT+"/"+key+".pe.qc.fl.fq.gz"
			]))

			## Save output file names
			OUTPUT_DICT[key] = ["./"+OUTPUT+"/"+key+".pe.qc.fl.fq.gz","./"+OUTPUT+"/"+key+".se.qc.fl.fq.gz"]

		else:
			## Re-Compress Files and Move to Output Folder
			os.system(' '.join([
				"gzip -9c combined-trim.fq.pe > combined-trim.fq.pe.gz",
			]))

			os.system(' '.join([
				"mv combined-trim.fq.pe.gz",
				"./"+OUTPUT+"/"+key+".pe.qc.fq.gz"
			]))

			os.system(' '.join([
				"gzip -9c combined-trim.fq.se",
				" > combined-trim.fq.se.gz"
			]))

			os.system(' '.join([
				"mv combined-trim.fq.se.gz",
				"./"+OUTPUT+"/"+key+".se.qc.fq.gz"
			]))

			## Save output file names
			OUTPUT_DICT[key] = ["./"+OUTPUT+"/"+key+".pe.qc.fq.gz","./"+OUTPUT+"/"+key+".se.qc.fq.gz"]

	## Remove Temp Folder "TRIM" and clean up
	os.system(' '.join([
		"rm -fr TRIM",
	]))

	os.system(' '.join([
		"rm combined-trim.fq.pe combined-trim.fq.se",
	]))

	if FLASH:
		os.system(' '.join([
			"rm out.extended* out.*hist* combined-trim.fq.se.gz",
		]))

##############################
# K-MER SELECTION w/ KMERGENIE (or user specified)
##############################
GENIE = "Yes"
KMER = "39"
KMER_DICT = {}

if len(CUSTOM) !=0:
	print "\n\nK-mer Genie will be used to select optimum k-mer length for assembly by default. If you would prefer to enter the k-mer length, or later, specify a range of k-mer lengths, Enter \"No\"\n"
	GENIE = raw_input()		

if GENIE != "No" or GENIE != "NO" or GENIE != "no" or GENIE != "N" or GENIE != "n":
	for key, value in OUTPUT_DICT.iteritems():

		#Perform KMERGENIE on the paired end file (assuming the single end reads file will be comparable)
		#KMERGENIE prints to standard output, which we re-direct to a text file
		os.system(' '.join([
			"kmergenie",
			value[0],
	                ">",
	                "./bestk.log",
        	        "2>&1"
		]))

		if LOG:
			command.write(' '.join([
				"kmergenie",
				value[0],
	                	">",
		                "./bestk.log",
        		        "2>&1"
			]))

		for line in open("bestk.log", "r"):
			if re.search("best k:", line):
				line = line.strip("\n\r")
				line = line.split()
				KMER_DICT[key] = line[2]

	if LOG:
		log.write("You based your k-mer lenght off of KMERGENIE's output and are using k-mer length: "+str(KMER_DICT)+"for your assembly\n\n")

	os.system(' '.join([
		"rm",
                "histograms*"
	]))

if GENIE == "No" or GENIE == "NO" or GENIE == "no" or GENIE == "N" or GENIE == "n":
	print "\n\nWould you like to iterate over multiple k-mer lengths? [YES/NO]\n"
	ITERATE = raw_input()		
	
	if ITERATE == "Yes" or ITERATE == "YES" or ITERATE == "yes" or ITERATE == "y" or ITERATE == "Y":
		MIN_KMER = ''
		MAX_KMER = ''
		KMER_STEP = ''

		print "\nSelect the range of k-mer lengths. ENTER mininum k-mer length now: "
		MIN_KMER = raw_input()
		print "\nSelect the range of k-mer lengths. ENTER maximum k-mer length now: "
		MAX_KMER = raw_input()
		print "\nSelect the range of k-mer lengths. ENTER iterative step length now: "
		KMER_STEP = raw_input()
		KMER = MIN_KMER+","+MAX_KMER+","+KMER_STEP
	
		if LOG:
			log.write("You will iterate over multiple k-mer lengths according to the following information: "+str(KMER)+"\n\n")
	else:
		print "ENTER k-mer value to use for assembly: "
		KMER = raw_input()

		if LOG:
			log.write("You manually selected your k-mer lenght and are using k-mer length: "+str(KMER)+"for your assembly\n\n")

if len(ITERATE) == 0 and MERGE or ITERATE == "NO|no|No" and MERGE:
	print "\n\nYou must perform an iteration over multiple k-mer lengths to perform a super assembly\n"
	MIN_KMER = ''
	MAX_KMER = ''
	KMER_STEP = ''

	print "\nSelect the range of k-mer lengths. ENTER mininum k-mer length now: "
	MIN_KMER = raw_input()
	print "\nSelect the range of k-mer lengths. ENTER maximum k-mer length now: "
	MAX_KMER = raw_input()
	print "\nSelect the range of k-mer lengths. ENTER iterative step length now: "
	KMER_STEP = raw_input()
	KMER = MIN_KMER+","+MAX_KMER+","+KMER_STEP

	if LOG:
		log.write("You were forced to iterate over multiple k-mer lengths b/c you plan on using MeGAMerge and will iterate based off of the folloing information: "+str(KMER)+"\n\n")
	
###########
# ASSEMBLY
###########

print "\n\n-- Performing Assembly --"

if LOG == 'ON':
	log.write("\n\n-- Performing Assembly --\n\n")
	command.write("\n\n-- Performing Assembly --\n\n")

############################
# LOADING CORRECT FILE NAMES
############################

FILE_DICT = OUTPUT_DICT

## Log Names of Files Used in Assembly
if LOG == 'ON':
	log.write("\n\nUsing the Following Files For Assembly: \n"+str(FILE_DICT)+"\n")
	command.write("\n\nUsing the Following Files For Assembly: \n"+str(FILE_DICT)+"\n")

#########
# IDBA-UD
#########

## NOTE: IDBA-UD does NOT make use of unpaired reads

if re.search("1|IDBA", ASSEMBLER):
	print "\n\n-----Using IDBA-----\n"

	if LOG == 'ON':
		log.write("\nNOTE: IDBA-UD does NOT make use of unpaired reads. FURTHER, by default it relies on abundance information of reads. If you are using Digital Normalization, be sure to run the script in CUSTOM mode (-C Y) and specify \"-- no_coverage\".\n\nFurther, IDBA is an iterative assembler and thus you do not need to specify a k-mer length, though you can select the minimum and maximum lengths to try.")

	####
	## Create Assembly Directory	
	####
	if os.path.exists('./' + OUTPUT + "/ASSEMBLY/IDBA/"):
        	print "\n\n-- IDBA-UD Assembly Folder Exists - Caution: Files May Have Been Over-written --"

	        if LOG == 'ON':
			log.write("\n-- IDBA-UD Assembly Folder Exists - Caution: Files May Have Been Over-written --\n")
	else:
       	        os.makedirs('./' + OUTPUT + "/ASSEMBLY/IDBA/")

	############################
	# Preparing Input Parameters
	############################

	## Enable User To Tailor Assembly Parameters
	if CUSTOM:
		## Clear Some Space and Print the List of Parameters Offered by IDBA
		print "\n\n\n"
		os.system(' '.join([			
			"idba"
		]))

		print "\n\nAbove Are The Options Afforded by IDBA_UD.\n\nYou will now be given the opportunity to specify argument for either long read (-l) or short read (-r) data, not both (unless you want to modify this script. Example input might be: \"-- mink 45 --maxk 95 --min_support 2\".\n\nIMPORTANT: IDBA offers few options that do not involve a specific read set. If you ARE changing parameters, be sure you are running the appropriate files through the pipeline. Whatever arguments you specify will be ADDED to the names of the read files already stored based on your input .tsv. NOTE: The Output Directory Cannot Be Specified.\n\nENTER input(s) for short read data \"-r\": "
		R = raw_input()

		print "ENTER input(s) for long read data \"-l\": "
		L = raw_input()

	for key, value in FILE_DICT.iteritems():

		###################################
		## Create Sample Specific Directory
		###################################
		if os.path.exists('./' + OUTPUT + "/ASSEMBLY/IDBA/" + key):
			print "\n\n-- IDBA Assembly Folder Exists - Caution: Files May Have Been Over-written --"

			if LOG == 'ON':
				log.write("\n-- IDBA Assembly Folder Exists - Caution: Files May Have Been Over-written --\n")
		else:
			os.makedirs('./' + OUTPUT + "/ASSEMBLY/IDBA/" + key)

		## Set-up Custom Arguments (if none, specify empty arrays)
		try:
			R
		except NameError:
			R = ''

		try:
			L
		except NameError:
			L = ''

		if R or L:
			if LOG == 'ON':
				log.write("You Ran IDBA with the following additional parameters: \n"+R+"\n"+L+"\n")

		##################
		# Run the Analysis
		##################	

		## Fork to Assemble Reads for one file only, two files or multiple sets of files.
		if np.size(value) == 2:
			# Case 1
			if np.size(value[0]) == 2:
				if LOG == 'ON':
					log.write("\n\nCase: 1\n")

				## Makes files are in correct format
				if re.search(".gz", value[0][0]):		
			                os.system(' '.join([
		        	                "gunzip",
		                	        value[0][0],
			                ]))

					PE_READS = Clean_GZ(value[0][0])

					if LOG == 'ON':
						log.write("\n\nWe unzipped your files: "+str(PE_READS)+"\n")
						command.write(' '.join(["gunzip",value[0][0]])+"\n")
				else:
					PE_READS = value[0][0]

			else:
				# Case 2
				if LOG == 'ON':
					log.write("\n\nCase: 2\n")

				## Makes files are in correct format
				if re.search(".gz", value[0]):		
		        	        os.system(' '.join([
		                	        "gunzip",
		                        	value[0],
			                ]))

					PE_READS = Clean_GZ(value[0])

					if LOG == 'ON':
						log.write("\n\nWe unzipped your files: "+str(PE_READS)+"\n")
						command.write(' '.join(["gunzip",value[0]])+"\n")

					else:
						PE_READS = value[0]

			##If necessary convert from .fastq to .fasta
			PE_READS = Convert_FQ(PE_READS)

			if len(L) > 0:
				if LOG == 'ON':
					command.write(' '.join(["\nidba_ud","-o","./"+OUTPUT+"/ASSEMBLY/IDBA/","-l",PE_READS,"--pre_correction",L])+"\n")		

				os.system(' '.join([
					"idba_ud",
					"-o",
					"./"+OUTPUT+"/ASSEMBLY/IDBA/",
					"-l",
					PE_READS,
					L
				]))		

			else:
				if LOG == 'ON':
					command.write(' '.join(["\nidba_ud","-o","./"+OUTPUT+"/ASSEMBLY/IDBA/","-r",PE_READS,"--pre_correction",R])+"\n")		

				os.system(' '.join([
					"idba_ud",
					"-o",
					"./"+OUTPUT+"/ASSEMBLY/IDBA/",
					"-r",
					PE_READS,
					R
				]))		

		## Perform Assembly On Multiple Different Files Associated with One Sample (i.e. Partitions) 
		else:				
			# Case 3
			if LOG == 'ON':
				log.write("\n\nCase: 3\n")

			for Pair in value:
				## Makes files are in correct format
				if re.search(".gz", Pair[0]):		
			                os.system(' '.join([
			                        "gunzip",
			                        Pair[0],
			                ]))

					PE_READS = Clean_GZ(Pair[0])

					if LOG == 'ON':
						log.write("\n\nWe unzipped your files: "+str(PE_READS)+"\n")
						command.write(' '.join(["gunzip",Pair[0]])+"\n")

				else:
					PE_READS = Pair[0]

				##If necessary convert from .fastq to .fasta
				PE_READS = Convert_FQ(PE_READS)

				if LOG == 'ON':
					log.write("\n\nNow Assembling: "+str(PE_READS)+"\n")
	
				if len(L) > 0:
					if LOG == 'ON':
						command.write(' '.join(["\nidba_ud","-o","./"+OUTPUT+"/ASSEMBLY/IDBA/","-l",PE_READS,"--pre_correction",L])+"\n")		

					os.system(' '.join([
						"idba_ud",
						"-o",
						"./"+OUTPUT+"/ASSEMBLY/IDBA/",
						"-l",
						PE_READS,
						"--pre_correction",
						L
					]))		

				else:
					if LOG == 'ON':
						command.write(' '.join(["\nidba_ud","-o","./"+OUTPUT+"/ASSEMBLY/IDBA","-r",PE_READS,"--pre_correction",R])+"\n")		

					os.system(' '.join([
						"idba_ud",
						"-o",
						"./"+OUTPUT+"/ASSEMBLY/IDBA/",
						"-r",
						PE_READS,
						"--pre_correction",
						R
					]))		

		###################################
		## Move Files to Sample Specific Directory
		###################################

		os.system(' '.join([
			"mv",
			'./' + OUTPUT + "/ASSEMBLY/IDBA/align*",
			'./' + OUTPUT + "/ASSEMBLY/IDBA/*.fa",
			'./' + OUTPUT + "/ASSEMBLY/IDBA/kmer",
			'./' + OUTPUT + "/ASSEMBLY/IDBA/begin",
			'./' + OUTPUT + "/ASSEMBLY/IDBA/end",
			'./' + OUTPUT + "/ASSEMBLY/IDBA/log",
			'./' + OUTPUT + "/ASSEMBLY/IDBA/" + key
		]))

	print "\n\n--- Assembly Completed --- \n\nYou've Completed Assembly with IDBA.The files you are most interested in are entitled \"contig.fa\".\n" 
	print "\nThe run statistics can be found in the \"log\" file.\n\n" 
	if LOG == 'ON':
		log.write("\n\n--- Assembly Completed ---\n\nYou've Completed Assembly with IDBA.The files you are most interested in are entitled \"contig.fa\".\n")
		log.write("\nThe run statistics can be found in the \"log\" file.\n\n")

############
# Velvet
############

if re.search("2|Velvet", ASSEMBLER):
	## Set Variable Names
	DIRECTORY_LISTING = []

	print "\n\n-----Using Velvet-----\n"

	if LOG == 'ON':
		log.write("\nNOTE: Velvet uses unpaired reads and so we shall. FURTHER, by default Velvet relies on abundance information of reads. It is therefore an open question as to how Digital Normalization will effect results. It is possible, according to khmer to re-inflate partitions which have been digital normalized, but the script they reference was missing in the latest release.\n\nFurther, Velvet uses a single k-mer length, we will implement an interative process, which will require user input for selecting the best length.\n\n")

	############################
	## Create Main Assembly Directory	
	############################
	if os.path.exists('./' + OUTPUT + "/ASSEMBLY/VELVET/"):
		print "\n\n-- VELVET Assembly Folder Exists - Caution: Files May Have Been Over-written --"
	
		if LOG == 'ON':
			log.write("\n-- VELVET Assembly Folder Exists - Caution: Files May Have Been Over-written --\n")
	else:
		os.makedirs('./' + OUTPUT + "/ASSEMBLY/VELVET/")
	
	############################
	# Preparing Input Parameters
	############################

	## Enable User To Tailor Assembly Parameters
	if CUSTOM:
		print "There are multiple stages requiring input for Velvet"
	
		## Clear Some Space and Print the List of Parameters Offered by IDBA
		print "\n\n\nWe offer assisted guidance for specifying important parameters throughout. However we do provide a place for the use of any parameters. [We DO NOT recommend altering any defaults]\n\n"

		###########
		## Velvet-H
		###########

		print "\n\n-- K-mer Hash Table Construction --\n\n"
			
		os.system(' '.join([			
			"velveth"
		]))
	
		print "\n\nENTER any additional parameters for \"velveth\". NOTE: All of your upstream read files will be automatically specified. Do not specify additional reads, unless they are coming from outside of this analysis.\nENTER now, if no preference, hit ENTER: "
		VELVETH_XTRA_PARAMS = raw_input()
	
		##########
		##Velvet-G
		#########
	
		print "\n\n-- De Bruijn Graph Construction --\n\n"
	
		os.system(' '.join([			
			"velvetg"
		]))
	
		print "\n\nAbove Are All Options Afforded by \"velvetg\"\n\nYou will now have to select an insert lenght. If left blank, the \"velvetg\" default is to treat all paired reads as unpaired. Example: 200. ENTER now: "
		INSERT_LENGTH = raw_input()
				
		print "\n\nENTER any additional parameters. NOTE: All of your upstream read files will be automatically specified. Do not specify additional reads, unless they are coming from outside of this analysis.\nENTER now, if NONE: hit ENTER: "
		VELVETG_XTRA_PARAMS = raw_input()
	
	## Set-up Custom Arguments (if none, specify defaults)
	try:
		INSERT_LENGTH
	except NameError:
		INSERT_LENGTH = '200'

	if len(INSERT_LENGTH) == 0:
		INSERT_LENGTH = '200'
	
	try:
		VELVETH_XTRA_PARAMS
	except NameError:
		VELVETH_XTRA_PARAMS = ''
	
	try:
		VELVETG_XTRA_PARAMS
	except NameError:
		VELVETG_XTRA_PARAMS = ''

	##################
	# Run the Analysis
	##################	
	
	for key, value in FILE_DICT.iteritems():

		###################
		# Set KMER ARGUMENT
		###################

		# In the former case KMER has been set as a constant, the latter due to KMERGENIE, varies with each file
		if ITERATE == "Yes" or ITERATE == "YES" or ITERATE == "yes" or ITERATE == "y" or ITERATE == "Y" or MERGE:
			KMER = KMER
		else:
			KMER = KMER_DICT[key]
	
		###################################
		## Create Sample Specific Directory	
		###################################
		if not ITERATE:
			if os.path.exists('./' + OUTPUT + "/ASSEMBLY/VELVET/" + key):
				print "\n\n-- VELVET Assembly Folder Exists - Caution: Files May Have Been Over-written --"
	
				if LOG == 'ON':
					log.write("\n-- VELVET Assembly Folder Exists - Caution: Files May Have Been Over-written --\n")
			else:
				os.makedirs('./' + OUTPUT + "/ASSEMBLY/VELVET/" + key)
	
		##########
		## Velvet-H
		##########
		if LOG == 'ON':
			command.write(str(VELVETH_XTRA_PARAMS)+"\n")

		if np.size(value[0]) == 2:
			## Makes files are in correct format
			if re.search(".gz", value[0][0]):		
		                os.system(' '.join([
		                        "gunzip",
		                        value[0][0],
		                        value[0][1]
		                ]))

				PE_READS = Clean_GZ(value[0][0])
				SE_READS = Clean_GZ(value[0][1])

				if LOG == 'ON':
					log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
					command.write(' '.join(["gunzip",value[0][0],value[0][1]])+"\n")

			else:
				PE_READS = value[0][0]
				SE_READS = value[0][1]

		else:
			## Makes files are in correct format
			if re.search(".gz", value[0]):		
		                os.system(' '.join([
		                        "gunzip",
		                        value[0],
		                        value[1]
		                ]))

				PE_READS = Clean_GZ(value[0])
				SE_READS = Clean_GZ(value[1])

				if LOG == 'ON':
					log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
					command.write(' '.join(["gunzip",value[0],value[1]])+"\n")

			else:
				PE_READS = value[0]
				SE_READS = value[1]

		##If necessary convert from .fastq to .fasta
		PE_READS = Convert_FQ(PE_READS)
		SE_READS = Convert_FQ(SE_READS)

		if LOG == 'ON':
			log.write("\n\nNow Creating Hash Tables for: "+str(PE_READS)+" and "+str(SE_READS)+"\n")

		if LOG == 'ON':
			command.write(' '.join([
				"\nvelveth",
				'./' + OUTPUT + "/ASSEMBLY/VELVET/" + key,
				KMER,
				"-shortPaired",
				PE_READS,
				"-short",
				SE_READS,
				VELVETH_XTRA_PARAMS+"\n"
			]))
	
		os.system(' '.join([
			"velveth",
			'./' + OUTPUT + "/ASSEMBLY/VELVET/" + key,
			KMER,
			"-shortPaired",
			PE_READS,
			"-short",
			SE_READS,
			VELVETH_XTRA_PARAMS
		]))

		if ITERATE:
			DIRECTORY_LISTING.append('./' + OUTPUT + "/ASSEMBLY/VELVET")
		else:
			DIRECTORY_LISTING.append('./' + OUTPUT + "/ASSEMBLY/VELVET/" + key)

		##########
		## Velvet-G
		##########
		
		if LOG == 'ON':
			log.write("You Used the Following Directories to Construct De Bruijn Graphs: \n"+str(DIRECTORY_LISTING)+"\n")
	
		## Velvet-G takes directories as inputs. Here we have to fork according to the presence of multiple "child" directories associated with product of multiple k-mer lengths
		if ITERATE:
			for PARENT in DIRECTORY_LISTING:
				for CHILD in range(int(MIN_KMER),int(MAX_KMER),int(KMER_STEP)):

					if LOG == 'ON':
						command.write(' '.join([
							"\nvelvetg",
							PARENT+"/"+key+"_"+str(CHILD),
							"-exp_cov auto",
							"-ins_length",
							INSERT_LENGTH,
							VELVETG_XTRA_PARAMS+"\n"
						]))
	
					os.system(' '.join([
						"velvetg",
						PARENT+"/"+key+"_"+str(CHILD),
						"-exp_cov auto",
						"-ins_length",
						INSERT_LENGTH,
						VELVETG_XTRA_PARAMS
					]))

		else:
			for PARENT in DIRECTORY_LISTING:
				if LOG == 'ON':
					command.write(' '.join([
						"\nvelvetg",
						PARENT,
						"-exp_cov auto",
						"-ins_length",
						INSERT_LENGTH,
						VELVETG_XTRA_PARAMS+"\n"
					]))

				os.system(' '.join([
					"velvetg",
					PARENT,
					"-exp_cov auto",
					"-ins_length",
					INSERT_LENGTH,
					VELVETG_XTRA_PARAMS
				]))
	

	print "\n\n--- Assembly Completed ---\n\nYou've Completed Assembly with Velvet.\n\n" 
	if LOG == 'ON':
		log.write("\n\n--- Assembly Completed ---\n\nYou've Completed Assembly with Velvet.\n\n")
		log.write("\nThe files you are most interested in are entitled \"contigs.fa\".")

#########
# RAY
#########

if re.search("3|RAY", ASSEMBLER):

	## Detail info for RAY use
	if LOG == 'ON':
		log.write("\nNOTE: RAY is a very straightfoward assembler, but has many additional perks for community analysis. It uses both paired and unpaired reads. Like the other two assemblers, you may want to carefully consider whether to use digital normalization, since read abundance is used in assmebly. RAY (or simply just RAY) does not automatically provide iteration for different k-mer assemblies, but this script will enable you to do so, if you chose.")

	## Set Variable Names
	DIRECTORY_LISTING = []

	print "\n\n-----Using RAY-----\n"

	############################
	## Create Main Assembly Directory	
	############################
	if os.path.exists('./' + OUTPUT + "/ASSEMBLY/RAY/"):
        	print "\n\n-- RAY Assembly Folder Exists - Caution: Files May Have Been Over-written --"

		if LOG == 'ON':
			log.write("\n-- RAY Assembly Folder Exists - Caution: Files May Have Been Over-written --\n")
	else:
       	        os.makedirs('./' + OUTPUT + "/ASSEMBLY/RAY/")

	if LOG == 'ON':
		log.write("\nNOTE: RAY uses unpaired reads and so we shall.\n\n")

	############################
	# Preparing Input Parameters
	############################
	
	## Enable User To Tailor Assembly Parameters
	if CUSTOM:
		print "\n\nThere are few \"vanilla\" parameter requirements for RAY, unless the user has specific needs.\n"
	
		os.system(' '.join([			
			"mpiexec Ray"
		]))
	
		print "\n\nAbove Are All Options Afforded by RAY. You will be given the option to provivde any additional arguments: \n"
	
		print "\n\nENTER the number of processors you have available now, hit enter if no preference [default: 12]: "
		THREADS = raw_input()

		print "\n\nENTER any additional parameters. NOTE: All of your upstream read files will be automatically specified. Do not specify additional reads, unless they are coming from outside of this immediate pipeline.\nENTER now, if NONE: hit ENTER: "
		RAY_XTRA_PARAMS = raw_input()
	
	## Set-up Custom Arguments (if none, specify defaults)
	try:
		THREADS
	except NameError:
		THREADS = '12'

	if len(THREADS) == 0:
		THREADS = '8'

	try:
		RAY_XTRA_PARAMS
	except NameError:
		RAY_XTRA_PARAMS = ''
	

	##################
	# Run the Analysis
	##################	

	for key, value in FILE_DICT.iteritems():
		# In the former case KMER has been set as a constant, the latter due to KMERGENIE, varies with each file
		if ITERATE == "Yes" or ITERATE == "YES" or ITERATE == "yes" or ITERATE == "y" or ITERATE == "Y" or MERGE:
			KMER = KMER
		else:
			KMER = KMER_DICT[key]

		if ITERATE:
			for ITER in range(int(MIN_KMER),int(MAX_KMER),int(KMER_STEP)):

				## Needs be string for join()
				ITER = str(ITER)

				## RAY will automatically create directories based on specified output (thank god)
				if np.size(value[0]) == 2:
					## Makes files are in correct format
					if re.search(".gz", value[0][0]):		
				                os.system(' '.join([
				                        "gunzip",
				                        value[0][0],
				                        value[0][1]
			        	        ]))

						PE_READS = Clean_GZ(value[0][0])
						SE_READS = Clean_GZ(value[0][1])

						if LOG == 'ON':
							log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
							command.write(' '.join(["gunzip",value[0][0],value[0][1]])+"\n")

					else:
						PE_READS = value[0][0]
						SE_READS = value[0][1]
	
				else:
					## Makes files are in correct format
					if re.search(".gz", value[0]):		
				                os.system(' '.join([
				                        "gunzip",
				                        value[0],
				                        value[1]
			        	        ]))

						PE_READS = Clean_GZ(value[0])
						SE_READS = Clean_GZ(value[1])

						if LOG == 'ON':
							log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
							command.write(' '.join(["gunzip",value[0],value[1]])+"\n")
	
					else:
						PE_READS = value[0]
						SE_READS = value[1]


				# Convert from .fastq to .fasta if necessary
				PE_READS = Convert_FQ(PE_READS)
				SE_READS = Convert_FQ(SE_READS)
								
				## Unfortunately, RAY auto-detects files based on extensions, and does not take kindly to .pe and .se produced from extract-paired-reads.py
				PE_READS = EXTENSION_STRIP(PE_READS,'.pe')
				SE_READS = EXTENSION_STRIP(SE_READS,'.se')

				if LOG == 'ON':
					log.write("\n\nNow Assembling @"+ITER+": "+str(PE_READS)+" and "+str(SE_READS)+" using various k-mer lengths\n")

				if LOG == 'ON':
					command.write(' '.join([
						"\nmpiexec",
						"-n",
						THREADS,
						"Ray",
						"-o",
						"./"+OUTPUT+"/ASSEMBLY/RAY/"+key+"_"+ITER,
						ITER,
						"-i",
						PE_READS,
						"-s",
						SE_READS,
						RAY_XTRA_PARAMS+"\n"
					]))
	
				os.system(' '.join([
					"mpiexec",
					"-n",
					THREADS,
					"Ray",
					"-o",
					"./"+OUTPUT+"/ASSEMBLY/RAY/"+key+"_"+ITER,
					ITER,
					"-i",
					PE_READS,
					"-s",
					SE_READS,
					RAY_XTRA_PARAMS
				]))
	

			## RAY does not need to have directories pre-emptively created.
		else:
			if np.size(value[0]) == 2:
				if re.search(".gz", value[0][0]):		
			                os.system(' '.join([
			                        "gunzip",
			                        value[0][0],
			                        value[0][1]
		        	        ]))

					PE_READS = Clean_GZ(value[0][0])
					SE_READS = Clean_GZ(value[0][1])

					if LOG == 'ON':
						log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
						command.write(' '.join(["gunzip",value[0][0],value[0][1]])+"\n")

				else:
					PE_READS = value[0][0]
					SE_READS = value[0][1]
	
			else:
				## Makes files are in correct format
				if re.search(".gz", value[0]):		
			                os.system(' '.join([
			                        "gunzip",
			                        value[0],
			                        value[1]
		        	        ]))

					PE_READS = Clean_GZ(value[0])
					SE_READS = Clean_GZ(value[1])

					if LOG == 'ON':
						log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
						command.write(' '.join(["gunzip",value[0],value[1]])+"\n")
	
				else:
					PE_READS = value[0]
					SE_READS = value[1]

			##If necessary convert from .fastq to .fasta
			PE_READS = Convert_FQ(PE_READS)
			SE_READS = Convert_FQ(SE_READS)

			## Prepare file extensions for RAY
			PE_READS = EXTENSION_STRIP(PE_READS,'.pe')
			SE_READS = EXTENSION_STRIP(SE_READS,'.se')

			if LOG == 'ON':
				log.write("\n\nNow Assembling @"+KMER+": "+str(PE_READS)+" and "+str(SE_READS)+" using k-mer length: "+KMER+"\n")
	
			if LOG == 'ON':
				command.write(' '.join([
					"\nmpiexec",
					"-n",
					THREADS,
					"Ray",
					"-o",
					"./"+OUTPUT+"/ASSEMBLY/RAY/"+key + "_" + KMER,
					KMER,
					"-i",
					PE_READS,
					"-s",
					SE_READS,
					RAY_XTRA_PARAMS+"\n"
			]))
	
			os.system(' '.join([
				"mpiexec",
				"-n",
				THREADS,
				"Ray",
				"-o",
				"./"+OUTPUT+"/ASSEMBLY/RAY/"+key + "_" + KMER,
				KMER,
				"-i",
				PE_READS,
				"-s",
				SE_READS,
				RAY_XTRA_PARAMS
			]))

	print "\n\n--- Assembly Completed ---\n\nYou've Completed Assembly with RAY.\n" 
	print "\nThe files you are most interested in are entitled \"Contigs.fasta\" and \"Scaffolds.fasta\"."
	print "\nThe assembly statistics can be found in the \"OutputNumbers.txt\"."

	if LOG == 'ON':
		log.write("\n\n--- Assembly Completed ---\n\nYou've Completed Assembly with RAY.\n\n")
		log.write("\nThe files you are most interested in are entitled \"Contigs.fasta\" and \"Scaffolds.fasta\".")
		log.write("\nThe assembly statistics can be found in the \"OutputNumbers.txt\".")

############
# ABySS
############

if re.search("4|ABySS|Abyss|abyss", ASSEMBLER):
	## Set Variable Names
	DIRECTORY_LISTING = []

	print "\n\n-----Using ABySS-----\n"

	############################
	## Create Main Assembly Directory	
	############################
	if os.path.exists('./' + OUTPUT + "/ASSEMBLY/ABYSS/"):
		print "\n\n-- ABYSS Assembly Folder Exists - Caution: Files May Have Been Over-written --"
	
		if LOG == 'ON':
			log.write("\n-- ABYSS Assembly Folder Exists - Caution: Files May Have Been Over-written --\n")
	else:
		os.makedirs('./' + OUTPUT + "/ASSEMBLY/ABYSS/")
	
	############################
	# Preparing Input Parameters
	############################
	
	## Enable User To Tailor Assembly Parameters
	if CUSTOM:
		print "\n\nThere are few parameter to select from for ABySS.\n"
	
		os.system(' '.join([			
			"abyss-pe --help"
		]))
	
		print "\n\nAbove Are All Options Afforded by ABySS. You will be given the option to provivde any additional arguments: \n"
	
		print "\n\nENTER any additional parameters. NOTE: All of your upstream read files will be automatically specified. Do not specify additional reads, unless they are coming from outside of this immediate pipeline.\nENTER now, if NONE: hit ENTER: "
		ABYSS_XTRA_PARAMS = raw_input()
	
	## Set-up Custom Arguments (if none, specify defaults)
	try:
		ABYSS_XTRA_PARAMS
	except NameError:
		ABYSS_XTRA_PARAMS = ''
	
	##################
	# Run the Analysis
	##################	

	for key, value in FILE_DICT.iteritems():
		# In the former case KMER has been set as a constant, the latter due to KMERGENIE, varies with each file
		if ITERATE == "Yes" or ITERATE == "YES" or ITERATE == "yes" or ITERATE == "y" or ITERATE == "Y" or MERGE:
			KMER = KMER
		else:
			KMER = KMER_DICT[key]

		if ITERATE:
			for ITER in range(int(MIN_KMER),int(MAX_KMER),int(KMER_STEP)):

				###################################
				## Create Sample and KMER  Specific Directory
				###################################
				if os.path.exists('./' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "/" + ITER):
					print "\n\n-- ABYSS Assembly Folder Exists - Caution: Files May Have Been Over-written --"

					if LOG == 'ON':
						log.write("\n-- ABYSS Assembly Folder Exists - Caution: Files May Have Been Over-written --\n")
				else:
					os.makedirs('./' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "/" + ITER)

				## Needs be string for join()
				ITER = str(ITER)

				## ABYSS will automatically create directories based on specified output (thank god)
				if np.size(value[0]) == 2:
					## Makes files are in correct format
					if re.search(".gz", value[0][0]):		
				                os.system(' '.join([
				                        "gunzip",
				                        value[0][0],
				                        value[0][1]
			        	        ]))

						PE_READS = Clean_GZ(value[0][0])
						SE_READS = Clean_GZ(value[0][1])

						if LOG == 'ON':
							log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
							command.write(' '.join(["gunzip",value[0][0],value[0][1]])+"\n")

					else:
						PE_READS = value[0][0]
						SE_READS = value[0][1]
	
				else:
					## Makes files are in correct format
					if re.search(".gz", value[0]):		
				                os.system(' '.join([
				                        "gunzip",
				                        value[0],
				                        value[1]
			        	        ]))

						PE_READS = Clean_GZ(value[0])
						SE_READS = Clean_GZ(value[1])

						if LOG == 'ON':
							log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
							command.write(' '.join(["gunzip",value[0],value[1]])+"\n")
	
					else:
						PE_READS = value[0]
						SE_READS = value[1]


				# Convert from .fastq to .fasta if necessary
				PE_READS = Convert_FQ(PE_READS)
				SE_READS = Convert_FQ(SE_READS)
								
				## Unfortunately, ABYSS auto-detects files based on extensions, and does not take kindly to .pe and .se produced from extract-paired-reads.py
				PE_READS = EXTENSION_STRIP(PE_READS,'.pe')
				SE_READS = EXTENSION_STRIP(SE_READS,'.se')

				if LOG == 'ON':
					log.write("\n\nNow Assembling @"+ITER+": "+str(PE_READS)+" and "+str(SE_READS)+" using various k-mer lengths\n")

				if LOG == 'ON':
				        command.write(' '.join([
				                "\nabyss-pe",
				                "name=" + key,
				                "k=" + ITER,
				                "in=" + "'" + PE_READS + "'",
			        	        "se=" + "'" + SE_READS + "'"
			        	]))

			        os.system(' '.join([
			                "abyss-pe",
			                "name=" + key,
			                "k=" + ITER,
			                "in=" + "'" + PE_READS + "'",
			                "se=" + "'" + SE_READS + "'"
			        ]))

			## Move Files to Appropriate Folder
		        os.system(' '.join([
		                "mv",
		                key+"-*",
				'./' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "/" + ITER + "/"	
		        ]))

		        os.system(' '.join([
		                "mv coverage.hist",
				'./' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "/" + ITER + "/"	
		        ]))
	
		else:

			###################################
			## Create Sample Specific Directory
			###################################
			if os.path.exists('./' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "_" + KMER):
				print "\n\n-- ABYSS Assembly Folder Exists - Caution: Files May Have Been Over-written --"

				if LOG == 'ON':
					log.write("\n-- ABYSS Assembly Folder Exists - Caution: Files May Have Been Over-written --\n")
			else:
				os.makedirs('./' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "_" + KMER)

			if np.size(value[0]) == 2:
				if re.search(".gz", value[0][0]):		
			                os.system(' '.join([
			                        "gunzip",
			                        value[0][0],
			                        value[0][1]
		        	        ]))

					PE_READS = Clean_GZ(value[0][0])
					SE_READS = Clean_GZ(value[0][1])

					if LOG == 'ON':
						log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
						command.write(' '.join(["gunzip",value[0][0],value[0][1]])+"\n")

				else:
					PE_READS = value[0][0]
					SE_READS = value[0][1]
	
			else:
				## Makes files are in correct format
				if re.search(".gz", value[0]):		
			                os.system(' '.join([
			                        "gunzip",
			                        value[0],
			                        value[1]
		        	        ]))

					PE_READS = Clean_GZ(value[0])
					SE_READS = Clean_GZ(value[1])

					if LOG == 'ON':
						log.write("\n\nWe unzipped your files: "+str(PE_READS)+" and "+str(SE_READS)+"\n")
						command.write(' '.join(["gunzip",value[0],value[1]])+"\n")
	
				else:
					PE_READS = value[0]
					SE_READS = value[1]

			##If necessary convert from .fastq to .fasta
			PE_READS = Convert_FQ(PE_READS)
			SE_READS = Convert_FQ(SE_READS)

			## Prepare file extensions for ABYSS
			PE_READS = EXTENSION_STRIP(PE_READS,'.pe')
			SE_READS = EXTENSION_STRIP(SE_READS,'.se')

			if LOG == 'ON':
				log.write("\n\nNow Assembling @"+KMER+": "+str(PE_READS)+" and "+str(SE_READS)+" using k-mer length: "+KMER+"\n")

			if LOG == 'ON':
			        command.write(' '.join([
			                "\nabyss-pe",
			                "name=" + key,
			                "k=" + KMER,
			                "in=" + "'" + PE_READS + "'",
		        	        "se=" + "'" + SE_READS + "'"
		        	]))

		        os.system(' '.join([
		                "abyss-pe",
		                "name=" + key,
		                "k=" + KMER,
		                "in=" + "'" + PE_READS + "'",
		                "se=" + "'" + SE_READS + "'"
		        ]))

			## Move Files to Appropriate Folder
		        os.system(' '.join([
		                "mv",
		                key+"-*",
				'./' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "_" + KMER + "/"	
		        ]))

		        os.system(' '.join([
		                "mv coverage.hist",
				'./' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "_" + KMER + "/"	
		        ]))

	print "\n\n--- Assembly Completed ---\n\nYou've Completed Assembly with ABYSS.\n" 
	print "\nThe files you are most interested in are entitled \"*SAMPLEID*-scaffolds.fa\" which links to a particular *SAMPLEID*.fa file."
	print "\nThe assembly statistics can be found in the \"OutputNumbers.txt\"."

	if LOG == 'ON':
		log.write("\n\n--- Assembly Completed ---\n\nYou've Completed Assembly with ABYSS.\n\n")
		log.write("\nThe files you are most interested in are entitled \"Contigs.fasta\" and \"Scaffolds.fasta\".")
		log.write("\nThe assembly statistics can be found in the \"OutputNumbers.txt\".")

################
# Save BESTK.log to ASSEMBLY Folder
################
if GENIE != "No" or GENIE != "NO" or GENIE != "no" or GENIE != "N" or GENIE != "n":
	os.system(' '.join([
		"mv ./bestk.log",
		'./' + OUTPUT + "/ASSEMBLY/"
	]))

################################
# BLAST PHYLOGENETIC INFORMATION
################################

#Initialize container for all of the blast outputs
compiled = []
compiled_stats = []

#Create a list for ID and Fasta separately
for key in SEQUENCE_DICT.iteritems():

	## Get Assembled Data
	## NOTE: If multiple assemblies were performed, the script will use the last one. It won't matter, since if the genome is contaminated, it will be detected in all assemblies.

	if re.search("1|IDBA", ASSEMBLER):
		CONTIGS = './' + OUTPUT + "/ASSEMBLY/IDBA/" + key + "/" + "scaffold.fa"
		CON_DIR = './' + OUTPUT + "/ASSEMBLY/IDBA/" + key + "/"

	if re.search("2|Velvet", ASSEMBLER):
#		if ITERATE:
#			CONTIGS = ''
#			CON_DIR = 

#		else:
		CONTIGS = './' + OUTPUT + "/ASSEMBLY/VELVET/" + key + "/" + "contigs.fa"
		CON_DIR = './' + OUTPUT + "/ASSEMBLY/VELVET/" + key + "/"

	if re.search("3|RAY", ASSEMBLER):
#		if ITERATE:
#			CONTIGS = ''
#			CON_DIR = 

#		else:
		CONTIGS = './' + OUTPUT + "/ASSEMBLY/RAY/" + key + "_" + KMER + "/" + "Scaffolds.fasta"
		CON_DIR = './' + OUTPUT + "/ASSEMBLY/RAY/" + key + "_" + KMER + "/"

	if re.search("4|ABySS|Abyss|abyss", ASSEMBLER):
#		if ITERATE:
#			CONTIGS = ''
#			CON_DIR = 

#		else:
		CONTIGS = './' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "_" + KMER + "/" + key + "-scaffolds.fa"
		CON_DIR = './' + OUTPUT + "/ASSEMBLY/ABYSS/" + key + "_" + KMER + "/"
	#Get Validated Sequence
	VALID_SEQ = SEQUENCE_DICT[key]

        #Construct Blast.db for Each Genome
        os.system(' '.join([
                "makeblastdb",
                "-in",
                CONTIGS,
                "-dbtype nucl -out",
                CONTIGS + ".db"
        ]))

        #Blastn VALID_SEQ 
        os.system(' '.join([
                "blastn",
                "-query",
                VALID_SEQ,
                "-db",
		CONTIGS + ".db",
                "-out",
                CON_DIR +  key + "_16S_confirm.txt",
                "-outfmt 7",
                "-max_target_seqs 10"                   #Change the number of hits returned by blast here.
        ]))

        #Step iv) is to compile one file with all the important data in it.
        blast_out = open(CON_DIR +  key + "_16S_confirm.txt", 'r')
        lines = blast_out.readlines()
        compiled.append(lines[5:15])

	#A Relic of the original script which used only ABYSS
	#This will compile the assembly statistics 
	#Very useful, but did not make an effort to do this for all assemblers
	if re.search("4|ABySS|Abyss|abyss", ASSEMBLER):
		#Files requires for priting compiled assembly stats
	        statfile = "./" + key + "/" + key + "-stats"
        	stats = open(statfile, 'r')
	        stats_temp = stats.readlines()
        	compiled_stats.append(stats_temp[1:4])

#Output all of the compiled Blast hits as a .tsv
blast_comp = open('./' + OUTPUT + '/ASSEMBLY/' + 'compiled_16S_verification.txt', 'w')
for Line in compiled:
        length = range(0,len(Line))
        count2 = 0

        for count2 in length:
                blast_comp.write(Line[count2])

blast_comp.close()


if re.search("4|ABySS|Abyss|abyss", ASSEMBLER):
	#Output all of the compiled Blast hits as a .tsv
	assembly_stats_comp = open('assembly_stats_compiled.txt', 'a')
	for Line in compiled_stats:
        	length = range(0,len(Line))
	        count2 = 0

        	for count2 in length:
                	assembly_stats_comp.write(Line[count2])

	assembly_stats_comp.close()

log.close()
command.close()
end = timeit.default_timer()

print end - now
