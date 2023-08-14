#!/usr/local/bin/python3

import shutil
import os
import subprocess
import argparse
import glob
import zipfile
import math
import sys
import json
import re

PARSER = argparse.ArgumentParser(
prog='CURED.py',
description='Classification Using Restriction Enzyme Diagnostics'
)

PARSER.add_argument(
"--species",
type=str,
nargs= '+',
help='Genus or species of interest'
)

PARSER.add_argument(
"-c",
"--case_control",
type=str,
help='text file with list of accession numbers to be downloaded from NCBI in case vs. control format'
)

PARSER.add_argument(
"-l",
"--list_of_case_accessions",
type=str,
help='text file with list of accession numbers that are your cases'
)

PARSER.add_argument(
"-g",
"--genomes",
type=str,
help='path to genomes'
)

PARSER.add_argument(
"--mlst",
action="store_true",
help='run mlst on genomes of interest to determine sequence type'
)

PARSER.add_argument(
"-st",
"--sequence_type",
type=int,
help='sequence type of interest'
)

PARSER.add_argument(
"-s",
"--minsupp",
type=int,
help='fsm-lite -s option. Sets sensitivity (default 2)'
)

PARSER.add_argument(
"-S",
"--maxsupp",
type=int,
help='fsm-lite -S option. Sets specificity (default inf)'
)

PARSER.add_argument(
"-m",
"--min",
type=int,
help='Minimum k-mer size to find. Default is 20.'
)

PARSER.add_argument(
"-M",
"--max",
type=int,
help='Maximum k-mer size to find. Default is 20.'
)

PARSER.add_argument(
'--summary',
action="store_true",
help='Run ncbi-datasets summary option to see how many genomes are to be downloaded prior to downloading. '
)

PARSER.add_argument(
'--sensitivity',
type=int,
help='Sensitivity threshold.'
)

PARSER.add_argument(
'--specificity',
type=int,
help='Specificity threshold.'
)

PARSER.add_argument(
'--keep_tmp',
action='store_true',
help='Option to keep the tmp files. Default is to remove all tmp files.'
)

PARSER.add_argument(
'-v',
'--version',
help='Print version and exit',
action='version',
version='%(prog)s 1.0'
)
ARGS = PARSER.parse_args()

if len(sys.argv) <= 1:
    print("No arguments selected. Please provide command line arguments. Exiting...")
    sys.exit(1)

if bool(vars(ARGS)["list_of_case_accessions"]) and not bool(vars(ARGS)["species"]):
    PARSER.exit(status=0, message="Error: You have to use -l with -s\n")
if bool(vars(ARGS)["mlst"]) and not bool(vars(ARGS)["sequence_type"]):
    PARSER.exit(status=0, message="Error: You have to use -st with -m\n")
if bool(vars(ARGS)["sequence_type"]) and not bool(vars(ARGS)["mlst"]):
    PARSER.exit(status=0, message="Error: You have to use -m with -st\n")
if bool(vars(ARGS)["summary"]) and not bool(vars(ARGS)["species"]):
    PARSER.exit(status=0, message="Error: You have to use --summary with --species\n")


#### Checking to see if fsm-lite, blast, and samtools are installed #####
try:
    fsmLite_check = subprocess.run(['fsm-lite', '--help'], capture_output=True, text=True)
    print("Found fsm-lite.")
except FileNotFoundError:
    PARSER.exit(status=0, message='ERROR: fsm-lite cannot be found. Please install.\n')

try:
    blast_check = subprocess.run(['blastn', '-h'], capture_output=True, text=True)
    print(f"Found blast.")
except FileNotFoundError:
    PARSER.exit(status=0, message='ERROR: blast cannot be found. Please install.\n')

try:
    samtools_check = subprocess.run(['samtools', '-h'], capture_output=True, text=True)
    print(f"Found samtools.")
except FileNotFoundError:
    PARSER.exit(status=0, message='ERROR: samtools cannot be found. Please install.\n')

#### If the user is using ncbi-datasets, check to make sure that the dependency is installed. #####
# TODO: Make sure version is 14.6.0 or higher
if not ARGS.genomes:
    try:
        datasets_version = subprocess.run(['datasets', '--version'], capture_output=True, text=True)
        print(f"Found ncbi-datasets: {datasets_version.stdout.strip()}")
    except FileNotFoundError:
        PARSER.exit(status=0, message='ERROR: ncbi-datasets cannot be found.\n')

##### If the user is using mlst, check to make sure that the dependency is installed. #########
if ARGS.mlst:
    try:
        mlst_version = subprocess.run(['mlst', '-version'], capture_output=True, text=True)
        print(f"Found mlst: {mlst_version.stdout.strip()}")
    except FileNotFoundError:
        PARSER.exit(status=0, message='ERROR: mlst cannot be found\n')

#### option to see the count of the number of genomes prior to downloading the genomes ####
# TODO: This should be an option by iteself (i.e. nothing else on the command line should be passed except for species)
if ARGS.summary:
    def extract_element_from_json(file_path, element_key):
        with open(file_path, 'r') as file:
            data = json.load(file)
            element = data.get(element_key)
            return element

    species = ARGS.species
    summary_cmd = subprocess.run([f'datasets summary genome taxon {species} > file.json'], capture_output=True, text=True, shell=True)

    json_file = 'file.json'
    element_key = 'total_count'

    extracted_element = extract_element_from_json(json_file, element_key)
    print('Number of genomes to be downloaded: ' + str(extracted_element))
    os.remove('file.json')

#### option to download genomes based on species of interest #####
# TODO: Make sure that this does not run if summary argument is passed
if ARGS.species:
        genus_species = ' '.join(ARGS.species)
        genus_species_command = '"' + genus_species + '"'
        run_datasets = subprocess.run([f"datasets download genome taxon {genus_species_command}"], shell=True, capture_output=True, text=True)
        print(run_datasets.stderr)

##### option to download genomes using a case vs control text file ######
CASES = []
CONTROLS = []
if ARGS.case_control:
    ids = []
    accession_numbers = ARGS.case_control
    with open(accession_numbers, 'r', encoding='utf-8-sig') as list:
        for line in list:
            line = line.rstrip()
            line_list = line.split(',')
            if line_list[-1] == 'case':
                CASES.append(line_list[0])
            else:
                if line_list[-1] == 'control':
                    CONTROLS.append(line_list[0])
            line = line.split(',')[0]
            ids.append(line)
    argument = ' '. join(ids)

    run_datasets = subprocess.run([f'datasets download genome accession {argument} --filename multiple_datasets.zip'], shell=True, capture_output=True, text=True)
    print(run_datasets.stderr)

##### option to download case genomes based on given accession numbers, as well as download genomes from a specific species #######
if ARGS.list_of_case_accessions:
    CASES = []
    case_accessions = ARGS.list_of_case_accessions
    with open(case_accessions, 'r') as cases:
        for case in cases:
            case = case.rstrip()
            CASES.append(case)

    if ARGS.species:
        genus_species = ' '.join(ARGS.species)
        genus_species_command = '"' + genus_species + '"'
        run_datasets = subprocess.run([f"datasets download genome taxon {genus_species_command}"], shell=True, capture_output=True, text=True)
        print(run_datasets.stderr)

##### extract fna files from all downloaded files ######
# TODO: Combine this block of code with above block of code that starts with if ARGS.case_control
if ARGS.case_control:

    # Function to extract files from a zip archive
    def extract_files_from_zip(zip_path, file_extension, destination_dir):
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            for member in zip_ref.infolist():
                if not member.is_dir() and member.filename.endswith(file_extension):
                    filename = os.path.basename(member.filename)
                    extracted_path = zip_ref.extract(member, path=destination_dir)
                    new_path = os.path.join(destination_dir, filename)
                    shutil.move(extracted_path, new_path)
                    print(f"Moved: {new_path}")

    # Specify the zip file
    zip_file = 'multiple_datasets.zip'

    # Specify the file extension
    file_extension = '.fna'

    # Specify the destination directory for the files to be moved to
    extractedFiles_directory = 'genomes'
    if not os.path.exists(extractedFiles_directory):
        os.mkdir(extractedFiles_directory)

    # Call the function to move files
    extract_files_from_zip(zip_file, file_extension, extractedFiles_directory)
    shutil.rmtree('genomes/ncbi_dataset/')

# if ARGS.species:
#     extract_command = subprocess.run(['unzip ncbi_datasets.zip'],shell=True)
#     print(extract_command.stdout)
#
#     destination_folder = 'genomes'
#     if not os.path.exists(destination_folder):
#         os.makedirs(destination_folder)
#
#     for dirpath, dirnames, filenames in os.walk('.'):
#         for filename in filenames:
#             if filename.endswith('.fna'):
#             # Build the source and destination paths
#                 source_path = os.path.join(dirpath, filename)
#                 destination_path = os.path.join(destination_folder, filename)
#             # Copy the file to the 'assemblies' folder
#                 shutil.copy2(source_path, destination_path)

############## run mlst on files in genomes directory ################
if ARGS.mlst:
    mlst_command = subprocess.run(['mlst genomes/*.fna > mlst_report.tsv'], shell=True)

    target_ST = []

    with open('mlst_report.tsv', 'r') as report:
        for line in report:
            line = line.rstrip()
            line_list = line.split()
            print(line_list)
            if line_list[2] == ARGS.sequence_type:
                target_ST.append(line_list[0])
    #print(target_ST)

# generate input file list required to run fsm-lite
if ARGS.genomes:
    directory = ARGS.genomes.rstrip('/')
else:
    directory = 'genomes'

# check to ensure that the directory exists before running the command
check_path = os.path.exists(directory)
if check_path == False:
    print('Directory does not exist. Exiting...')
    sys.exit(1)
else:
    pass

command = f'for f in {directory}/*.fna; do id=$(basename "$f" .fna); echo $id $f; done > input.list'
print('Preparing files for fsm-lite.')
inputList = subprocess.run(command, shell=True)

#TODO: Check that the names in the case_control file match those in the folder of genomes provided by the user.
# Block of code still needs to be corrected
if ARGS.genomes:
    genome_names = []
    for line in open('input.list', 'r'):
        line = line.rstrip()
        name = line.split(',')[0]
        genome_names.append(name)
    case_control_genome_names = CASES + CONTROLS
    sorted_genome_names = sorted(genome_names)
    sorted_case_control_genome_names = sorted(case_control_genome_names)
    for itemA, itemB in zip(sorted_genome_names, sorted_case_control_genome_names):
        if itemA != itemB:
            print('Names in case_control file do not match the names in the provided directory. Please correct. Exiting...')
            sys.exit(1)
        else:
            continue

print('Completed preparation. Running fsm-lite.')

# create command based on user-defined parameters for -s, -S, -m, -M. If no values are defined, default values are used. Default for kmer size is 20. Other default values come directly from fsm-lite.
if hasattr(ARGS, 'minsupp') and ARGS.minsupp:
    minsupp = int(ARGS.minsupp)
else:
    minsupp = 2

if hasattr(ARGS, 'maxsupp') and ARGS.maxsupp:
    maxsupp = int(ARGS.maxsupp)
else:
    maxsupp = math.inf

if hasattr(ARGS, 'min') and ARGS.min:
    min_val = int(ARGS.min)
else:
    min_val = 20

if hasattr(ARGS, 'max') and ARGS.max:
    max_val = int(ARGS.max)
else:
    max_val = 20

# Check to make sure that --min is smaller than or equal to --max (or else fsm-lite will raise error)
if min_val > max_val:
    print("-m,--min must be smaller than or equal to -M,--max. Exiting...")
    sys.exit(1)

# Check to make sure that --minsupp is smaller than or equal to --maxsupp (or else fsm-lite will raise error)
if minsupp > maxsupp:
    print("ERROR: -s,--minsupp must be smaller than or equal to -S,--maxsupp. Exiting...")
    sys.exit(1)

# Command is written based on whether or not a maxsupp value is defined by the user
if ARGS.maxsupp:
    fsmLite_cmd = 'fsm-lite -l input.list -m {} -M {} -S {} -s {} -t tmp | gzip - > output.txt.gz'.format(min_val, max_val, maxsupp, minsupp)
else:
    fsmLite_cmd = 'fsm-lite -l input.list -m {} -M {} -s {} -t tmp | gzip - > output.txt.gz'.format(min_val, max_val, minsupp)

# run fsm-lite
fsmLite = subprocess.run(fsmLite_cmd, shell=True)
print('fsm-lite is complete.')

# once fsm-lite is finished, parse the output file
print('Parsing output. Finding the unique kmers.')

# Set the thresholds for sensitivity and specificity. If the user does not define threshold values, use 100 sensitivity and specificity as default.
if ARGS.sensitivity:
    sensitivity = ARGS.sensitivity
else:
    sensitivity = 100
if ARGS.specificity:
    specificity = ARGS.specificity
else:
    specificity = 100

sensitivity = sensitivity / 100
sensitivity = round(sensitivity, 2)
specificity = specificity / 100
specificity = round(specificity, 2)

# Get the lengths of control list and case list. Necessary for parsing the output file.
cases = len(CASES)
controls = len(CONTROLS)

# Find range of allowed list lengths (effectively finding max false negative threshold for sensitivity)
sensitivity_step1 = cases / sensitivity
sensitivity_step2 = int(sensitivity_step1 - cases)
lower_limit = cases - sensitivity_step2
#print(lower_limit)

# Find maximum nuber of controls that can be found (effectively finding max false positive threshold for specificity)
specificity_step1 = controls / specificity
specificity_step2 = int(specificity_step1 - controls)
upper_limit = cases + specificity_step2

# fsm-lite output is gzipped so need to gunzip it first
unzip_command = subprocess.run(['gunzip output.txt.gz'], shell=True)

#TODO: Add in lines of code to check that there is information in the output file (sometimes it will just output an empty file and not return any errors)
check_file = os.stat('output.txt').st_size
if(check_file == 0):
    print('Fsm-lite output is empty. Exiting...')
    sys.exit(1)
else:
    pass

# Find lists with a length equal to or more than lower limit but not greater than upper limit. For those lists passing threshold, check that the number of false positives is not greater than specified limit (specificity_step2)
passing_sensitivity_threshold = []
for line in open('output.txt', 'r'):
    line_list = line.rstrip().split()
    kmer = line_list[0]
    genome_list = line_list[2:]
    if (len(genome_list) >= lower_limit) and (len(genome_list) <= upper_limit):
        passing_sensitivity_threshold.append(line_list)

passing_specificity_threshold = []
failing_specificity_threshold = []
for data_list in passing_sensitivity_threshold:
    false_positive_counter = 0
    kmer = data_list[0]
    separator = data_list[1]
    genome_list = data_list[2:]

    for genome in genome_list:
        genome = genome.split(':')[0]
        genomes_list = genome.split('_')
        genome = '_'.join(genomes_list[:2])
        for CONTROL in CONTROLS:
            if genome == CONTROL:
                false_positive_counter += 1
    if false_positive_counter > specificity_step2:
        failing_specificity_threshold.append(data_list)
    else:
        data_list.append(false_positive_counter) # Append counter as last element in list to calculate specificity in final step
        passing_specificity_threshold.append(data_list)

# Enzyme dictionary. Taken from https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities and modified accordingly. For positions where multiple bases could be matches, modified the sequences so that they are in a format compatible with regex.
ENZYME_DICT = {'BsaAI': '[CT]ACGT[AG]', 'EaeI': '[CT]GGCC[AG]', 'BsaWI': '[AT]CCGG[AT]', 'PspXI': '[ACG]CTCGAG[CGT]', 'DraI': 'TTTAAA', 'PacI': 'TTAATTAA', 'PsiI-v2': 'TTATAA', 'BstBI': 'TTCGAA', 'PI-PspI': 'TGGCAAACAGCTATTATGGGTATTATGGGT', 'MscI': 'TGGCCA', 'FspI': 'TGCGCA', 'HpyCH4V': 'TGCA', 'Hpy188I': 'TC[ACTG]GA', 'NruI-HF®': 'TCGCGA', 'MmeI': 'TCC[AG]AC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'Hpy188III': 'TC[ACTG][ACTG]GA', 'I-SceI': 'TAGGGATAACAGGGTAAT', 'SnaBI': 'TACGTA', 'I-CeuI': 'TAACTATAACGGTCCTAAGGTAGCGAA', 'MseI': 'TTAA', 'BsrGI-HF®': 'TGTACA', 'BclI\xa0BclI-HF': 'TGATCA', 'XbaI': 'TCTAGA', 'TaqI-v2': 'TCGA', 'BspEI': 'TCCGGA', 'BspHI': 'TCATGA', 'HaeII': '[AG]GCGC[CT]', 'PpuMI': '[AG]GG[AT]CC[CT]', 'EcoO109I': '[AG]GG[ACTG]CC[CT]', 'CviKI-1': '[AG]GC[CT]', 'NspI': '[AG]CATG[CT]', 'BstYI': '[AG]GATC[CT]', 'BsrFI-v2': '[AG]CCGG[CT]', 'ApoI-HF': '[AG]AATT[CT]', 'TspRI': '[ACTG][ACTG]CA[CG]TG[ACTG][ACTG]', 'BsiHKAI': 'G[AT]GC[AT]C', 'HincII': 'GT[CT][AG]AC', 'PmeI': 'GTTTAAAC', 'HpaI': 'GTTAAC', 'Hpy166II': 'GT[ACTG][ACTG]AC', 'BsgI': 'GTGCAG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'Nt.BsmAI': 'GTCTC[ACTG]', 'BcoDI\xa0BsmAI': 'GTCTC[ACTG]', 'BciVI': 'GTATCC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'BstZ17I-HF®': 'GTATAC', 'AccI': 'GT[AC][GT]AC', 'RsaI': 'GTAC', 'BanII': 'G[AG]GC[CT]C', 'BsaHI': 'G[AG]CG[CT]C', 'BaeGI': 'G[GT]GC[AC]C', 'HphI': 'GGTGA[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'BsaI-HF®v2': 'GGTCTC[ACTG]', 'KpnI-HF®': 'GGTACC', 'NlaIV': 'GG[ACTG][ACTG]CC', 'ApaI': 'GGGCCC', 'BsmFI': 'GGGAC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'EciI': 'GGCGGA[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'PluTI': 'GGCGCC', 'SfiI': 'GGCC[ACTG][ACTG][ACTG][ACTG][ACTG]GGCC', 'FseI': 'GGCCGGCC', 'SfoI': 'GGCGCC', 'FokI': 'GGATG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'BtsCI': 'GGATG[ACTG][ACTG]', 'Nt.AlwI': 'GGATC[ACTG][ACTG][ACTG][ACTG]', 'AlwI': 'GGATC[ACTG][ACTG][ACTG][ACTG]', 'AscI': 'GGCGCGCC', 'NarI': 'GGCGCC', 'HaeIII': 'GGCC', 'Bsp1286I': 'G[AGT]GC[ACT]C', 'Nt.BspQI': 'GCTCTTC[ACTG]', 'BspQI\xa0SapI': 'GCTCTTC[ACTG]', 'BmtI-HF®': 'GCTAGC', 'MwoI': 'GC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]GC', 'Cac8I': 'GC[ACTG][ACTG]GC', 'BtgZI': 'GCGATG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'AsiSI': 'GCGATCGC', 'HhaI': 'GCGC', 'BglI': 'GCC[ACTG][ACTG][ACTG][ACTG][ACTG]GGC', 'NmeAIII': 'GCCGAG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'SrfI': 'GCCCGGGC', 'NaeI': 'GCCGGC', 'SphI-HF®\xa0SphI': 'GCATGC', 'SfaNI': 'GCATC[ACTG][ACTG][ACTG][ACTG][ACTG]', 'BstAPI': 'GCA[ACTG][ACTG][ACTG][ACTG][ACTG]TGC', 'Nb.BtsI': 'GCAGTG', 'BtsI-v2': 'GCAGTG[ACTG][ACTG]', 'BbvI': 'GCAGC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'Nb.BsrDI': 'GCAATG', 'BsrDI': 'GCAATG[ACTG][ACTG]', 'BlpI': 'GCT[ACTG]AGC', 'Fnu4HI': 'GC[ACTG]GC', 'NotI-HF®\xa0NotI': 'GCGGCCGC', 'BsaBI': 'GAT[ACTG][ACTG][ACTG][ACTG]ATC', 'EcoRV-HF®\xa0EcoRV': 'GATATC', 'MlyI': 'GAGTC[ACTG][ACTG][ACTG][ACTG][ACTG]', 'WarmStart®\xa0Nt.BstNBI\xa0Nt.BstNBI': 'GAGTC[ACTG][ACTG][ACTG][ACTG]', 'PleI': 'GAGTC[ACTG][ACTG][ACTG][ACTG]', 'BseRI': 'GAGGAG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'SacI-HF®': 'GAGCTC', 'Eco53kI': 'GAGCTC', 'DrdI': 'GAC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]GTC', 'AhdI': 'GAC[ACTG][ACTG][ACTG][ACTG][ACTG]GTC', 'PshAI': 'GAC[ACTG][ACTG][ACTG][ACTG]GTC', 'Tth111I\xa0PflFI': 'GAC[ACTG][ACTG][ACTG]GTC', 'AatII': 'GACGTC', 'HgaI': 'GACGC[ACTG][ACTG][ACTG][ACTG][ACTG]', 'ZraI': 'GACGTC', 'Nb.BsmI': 'GAATGC', 'BsmI': 'GAATGC[ACTG]', 'XmnI': 'GAA[ACTG][ACTG][ACTG][ACTG]TTC', 'BbsI-HF®\xa0BbsI': 'GAAGAC[ACTG][ACTG]', 'MboII': 'GAAGA[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'DpnI': 'GATC', 'ApaLI': 'GTGCAC', 'SalI\xa0SalI-HF®': 'GTCGAC', 'CviQI': 'GTAC', 'BanI': 'GG[CT][AG]CC', 'AvaII': 'GG[AT]CC', 'BstEII-HF®': 'GGT[ACTG]ACC', 'Acc65I': 'GGTACC', 'Sau96I': 'GG[ACTG]CC', 'PspOMI': 'GGGCCC', 'KasI': 'GGCGCC', 'BamHI-HF®\xa0BamHI': 'GGATCC', 'ApeKI\xa0TseI': 'GC[AT]GC', 'NheI-HF®': 'GCTAGC', 'BssHII': 'GCGCGC', 'HinP1I': 'GCGC', 'NgoMIV': 'GCCGGC', 'TfiI': 'GA[AT]TC', 'HinfI': 'GA[ACTG]TC', 'EcoRI-HF®\xa0EcoRI': 'GAATTC', 'BpuEI': 'CTTGAG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'BpmI': 'CTGGAG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'PstI-HF®\xa0PstI': 'CTGCAG', 'AcuI': 'CTGAAG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'EarI': 'CTCTTC[ACTG]', 'BspCNI': 'CTCAG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'SgrAI': 'C[AG]CCGG[CT]G', 'MspJI': 'C[ACTG][ACTG][AG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'AbaSI': 'C[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]G', 'MspA1I': 'C[AC]GC[GT]G', 'Hpy99I': 'CG[AT]CG', 'Esp3I': 'CGTCTC[ACTG]', 'BsmBI-v2': 'CGTCTC', 'BsiEI': 'CG[AG][CT]CG', 'PvuI-HF®': 'CGATCG', 'RsrII': 'CGG[AT]CCG', 'BstUI': 'CGCG', 'HpyAV': 'CCTTC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'EcoNI': 'CCT[ACTG][ACTG][ACTG][ACTG][ACTG]AGG', 'Bpu10I': 'CCT[ACTG]AGC', 'SbfI-HF®': 'CCTGCAGG', 'Nb.BbvCI': 'CCTCAGC', 'Nt.BbvCI': 'CCTCAGC', 'BbvCI': 'CCTCAGC', 'MnlI': 'CCTC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'BslI': 'CC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]GG', 'BsrBI': 'CCGCTC', 'SacII': 'CCGCGG', 'AciI': 'CCGC', 'LpnPI': 'CC[AGT]G[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'Nt.CviPII': 'CC[AGT]', 'FauI': 'CCCGC[ACTG][ACTG][ACTG][ACTG]', 'BseYI': 'CCCAGC', 'SmaI': 'CCCGGG', 'BccI': 'CCATC[ACTG][ACTG][ACTG][ACTG]', 'BstXI': 'CCA[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]TGG', 'XcmI': 'CCA[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]TGG', 'PflMI': 'CCA[ACTG][ACTG][ACTG][ACTG][ACTG]TGG', 'BstNI': 'CC[AT]GG', 'Bsu36I': 'CCT[ACTG]AGG', 'NciI': 'CC[CG]GG', 'ScrFI': 'CC[ACTG]GG', 'FspEI': 'CC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'MslI': 'CA[CT][ACTG][ACTG][ACTG][ACTG][AG]TG', 'NlaIII': 'CATG', 'BtsIMutI': 'CAGTG[ACTG][ACTG]', 'AlwNI': 'CAG[ACTG][ACTG][ACTG]CTG', 'EcoP15I': 'CAGCAG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'PvuII-HF®\xa0PvuII': 'CAGCTG', 'DraIII-HF®': 'CAC[ACTG][ACTG][ACTG]GTG', 'AleI-v2': 'CAC[ACTG][ACTG][ACTG][ACTG]GTG', 'BmgBI': 'CACGTC', 'Nb.BssSI': 'CACGAG', 'BssSI-v2': 'CACGAG', 'PaqCI': 'CACCTGC[ACTG][ACTG][ACTG][ACTG]', 'PmlI': 'CACGTG', 'NdeI': 'CATATG', 'AvaI\xa0BsoBI': 'C[CT]CG[AG]G', 'SmlI': 'CT[CT][AG]AG', 'AflII': 'CTTAAG', 'SfcI': 'CT[AG][CT]AG', 'DdeI': 'CT[ACTG]AG', 'XhoI\xa0PaeR7I': 'CTCGAG', 'BfaI': 'CTAG', 'BsiWI\xa0BsiWI-HF®': 'CGTACG', 'EagI-HF®': 'CGGCCG', 'StyI-HF®': 'CC[AT][AT]GG', 'AvrII': 'CCTAGG', 'BtgI': 'CC[AG][CT]GG', 'BsaJI': 'CC[ACTG][ACTG]GG', 'HpaII\xa0MspI': 'CCGG', 'TspMI\xa0XmaI': 'CCCGGG', 'NcoI-HF®\xa0NcoI': 'CCATGG', 'CviAII': 'CATG', 'MfeI-HF®': 'CAATTG', 'SwaI': 'ATTTAAAT', 'NsiI-HF®\xa0NsiI': 'ATGCAT', 'PI-SceI': 'ATCTATGTCGGGTGCGGAGAAAGAGGTAAT', 'AseI': 'ATTAAT', 'ClaI\xa0BspDI': 'ATCGAT', 'ScaI-HF®': 'AGTACT', 'StuI': 'AGGCCT', 'AfeI': 'AGCGCT', 'AluI': 'AGCT', 'BmrI': 'ACTGGG[ACTG][ACTG][ACTG][ACTG][ACTG]', 'BsrI': 'ACTGG[ACTG]', 'HpyCH4III': 'AC[ACTG]GT', 'BceAI': 'ACGGC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'BspMI\xa0BfuAI': 'ACCTGC[ACTG][ACTG][ACTG][ACTG]', 'SspI-HF®': 'AATATT', 'AclI': 'AACGTT', 'BglII': 'AGATCT', 'SpeI-HF®': 'ACTAGT', 'AflIII': 'AC[AG][CT]GT', 'HpyCH4IV': 'ACGT', 'MluI-HF®': 'ACGCGT', 'SexAI': 'ACC[AT]GGT', 'AgeI-HF®': 'ACCGGT', 'PciI': 'ACATGT', 'HindIII\xa0HindIII-HF®': 'AAGCTT', 'Tsp45I': 'GT[CG]AC', 'DpnII\xa0Sau3AI\xa0MboI': 'GATC', 'PspGI': 'CC[AT]GG', 'StyD4I': 'CC[ACTG]GG', 'FatI': 'CATG', 'MluCI': 'AATT', 'BsaXI': '[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]AC[ACTG][ACTG][ACTG][ACTG][ACTG]CTCC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'CspCI': '[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]CAA[ACTG][ACTG][ACTG][ACTG][ACTG]GTGG[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'BaeI': '[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]AC[ACTG][ACTG][ACTG][ACTG]GTA[CT]C[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]', 'BcgI': '[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]CGA[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]TGC[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]'}


# Write to a report the kmers that are passing both sensitivity and specificity thresholds.
##TODO: Need to make the CURED output report a unique name (i.e. add a time stamp, etc.)
cured_report = 'CURED_KmerReport.tsv'
unique_kmers = []
with open(cured_report, 'w') as out:
    out.write('Kmer\tGenomes With Kmer: Frequency of Kmer Occurrence\tSpecificity\tSensitivity\n')
    for LIST in passing_specificity_threshold:
        kmer = LIST[0]
        if kmer not in unique_kmers:
            unique_kmers.append(kmer)
        else:
            continue
        number_of_controls_found = LIST[-1] # NEED THIS VALUE TO CALCULATE FINAL SPECIFICITY
        GENOMES_LIST = LIST[2:-1]
        GENOMES = ' '.join(GENOMES_LIST)
        sensitivity = len(GENOMES_LIST) - number_of_controls_found # Need to subtract number of controls found in the genome list so to not get >100% sensitivity
        final_sensitivity = sensitivity / cases
        specificity = number_of_controls_found / len(GENOMES_LIST)
        unrounded_specificity = 1 - specificity
        final_specificity = round(unrounded_specificity, 2)

        out.write(f'{kmer}\t{GENOMES}\t{final_specificity}\t{final_sensitivity}\n')

                # match_found = False
                # for enzyme_name, enzyme_sequence in ENZYME_DICT.items():
                #     matches = re.findall(enzyme_sequence, kmer)
                #     for match in matches:
                #         out.write(f'{kmer}\t{GENOMES}\t{final_specificity}\t{final_sensitivity}\t{enzyme_name}\t{match}\n')
                #         match_found = True
                #
                # if not match_found:
                #     out.write(f'{kmer}\t{GENOMES}\t{final_specificity}\t{final_sensitivity}\tNo match\t-\n')
unique_kmer_count = str(len(unique_kmers))
print(f'{unique_kmer_count} unique kmers were found in your case genomes.')
# Function needed to find the unique restriction sites
def find_enzymes_with_different_positions(sequences_dict, enzymes_dict):
    enzymes_with_different_positions = []

    for enzyme, enzyme_sequence in enzymes_dict.items():
        enzyme_occurrences = {sequence_name: [] for sequence_name in sequences_dict}

        for sequence_name, sequence in sequences_dict.items():
            matches = re.finditer(enzyme_sequence, sequence)
            for match in matches:
                match_seq = match.group()
                start_position = match.start()
                end_position = match.end()
                enzyme_occurrences[sequence_name].append((enzyme, start_position, end_position))
        enzymes_with_different_positions.append(enzyme_occurrences)
    return enzymes_with_different_positions

# Function to remove temporary files that are created when finding unique restriction sites.
def remove_files_with_pattern(directory, pattern):
    matching_files = glob.glob(os.path.join(directory, pattern))

    for file_path in matching_files:
        try:
            os.remove(file_path)
            # print(f"File {file_path} removed successfully.")
        except OSError as e:
            print(f"Error while removing {file_path}: {e}")

# Function to check if a file exists and is not empty.
def check_file_exists_and_not_empty(file_path):
    if os.path.exists(file_path):
        if os.path.getsize(file_path) > 0:
            return True
        else:
            return False
    else:
        return False

# Get paths of controls to be used during the for loop used for finding unique restriction sites.
control_filepaths = []
for control in CONTROLS:
    for file in os.listdir(extractedFiles_directory):
        if control in file:
            filepath = os.path.join(extractedFiles_directory, file)
            control_filepaths.append(filepath)

# Get case that will be used to compare all of the controls. Eventually will need to check that the kmer of interest is in ths case.
case_genome_to_blast = CASES[0]

# Iterates over each unique kmer, blasts kmer in the chosen case. Then extracts +/-200 bps. Then, blasts that extracted region (containing the unique kmer) in each control. Then iterates over each extracted region, searching for unique restriction sites.
for kmer in unique_kmers:
    #print(kmer)
    sequence_dictionary = {}

# Writes each kmer to a file. This is necessary for blast.
    with open('uniqueKmer.fa', 'w') as out:
        out.write('>Kmer\n')
        out.write(kmer)

# Find the representative case in the genomes/ directory. Blast the kmer. Parameters needed to be tweaked because kmer size [20] is so small.
    for file in os.listdir(extractedFiles_directory):
        if case_genome_to_blast in file:
            filepath = os.path.join(extractedFiles_directory, file)
            genome_name = file.split('.fna')[0]
            blast_cmd = f'blastn -query uniqueKmer.fa -subject {filepath} -word_size 10 -evalue 1e-2 -outfmt "6 qseqid sseqid sstart send sstrand" > {genome_name}_blast.out'
            #print(blast_cmd)
            subprocess.run(blast_cmd, shell=True)
        else:
            continue

        check_case_exists = check_file_exists_and_not_empty(f'{genome_name}_blast.out')
        #print('Blast output for case exist:' + str(check_case_exists))
        if not check_case_exists:
            print(f'Blast could not find {kmer} in case')
            break
# Using the blast output, extract the region (+/- 200 bps) using samtools. Command is built using blast output variables.
        if check_case_exists:
            output_file = f'{genome_name}_samtools.fa'
            with open(f'{genome_name}_blast.out', 'r') as blast_out:
                header_blast_case = blast_out.readline()
                #print(header_blast_case)
                line_list = header_blast_case.rstrip('\n').split('\t')
                query_strand = line_list[0]
                # subj_strand = line_list[1]
                contig_loc = line_list[1]
                start = line_list[2]
                end = line_list[3]
                strand = line_list[4]

                if strand == 'plus':
                    start_adjusted = int(start) - 200
                    end_adjusted = int(end) + 200
                    samtools_extract_cmd = f'samtools faidx {filepath} {contig_loc}:{start_adjusted}-{end_adjusted} > {output_file}'
                    #print(samtools_extract_cmd)
                    subprocess.run(samtools_extract_cmd, shell=True)
                else:
                    if strand == 'minus':
                        start_adjusted = int(end) - 200
                        end_adjusted = int(start) + 200
                        rev_samtools_extract_cmd = f'samtools faidx {filepath} {contig_loc}:{start_adjusted}-{end_adjusted} -i > {output_file}'
                        #print(rev_samtools_extract_cmd)
                        subprocess.run(rev_samtools_extract_cmd, shell=True)
# Add extracted region for case to sequence dictionary. This will be used later on to find unique restriction sites.
            check_samtools_exists = check_file_exists_and_not_empty(output_file)
            #print('Samtools output for case exist: ' + str(check_samtools_exists))
            if check_samtools_exists:
                with open(output_file, 'r') as out:
                    header_line = output_file.split('_')
                    header = '_'.join(header_line[0:2])
                    contents = out.readlines()
                    seqs = []
                    for line in contents[1:]:
                        line = line.rstrip('\n')
                        seqs.append(line)
                        seq = ''.join(seqs)
                        sequence_dictionary[header] = seq

    if not check_case_exists:
        continue
# Iterate over controls, blast the extracted case region in each of the controls, if that region exists, extract +/- 200 bps using samtools. Add this to the sequence dictionary.
    for control in control_filepaths:
        control_basename = control.split('/')[1]
        control_basename =  control_basename.split('.fna')[0]
        blast_cmd = f'blastn -query {output_file} -subject {control} -outfmt "6 qseqid sseqid sstart send sstrand" > {control_basename}_blast.out'
        subprocess.run(blast_cmd, shell=True)
        #print(blast_cmd)
        check_control_exists = check_file_exists_and_not_empty(f'{control_basename}_blast.out')
        #print(check_control_exists)
        if check_control_exists:
            with open(f'{control_basename}_blast.out', 'r') as control_blast_out:
                header_of_blast_control = control_blast_out.readline()
                #print(header_of_blast_control)
                control_line = header_of_blast_control.rstrip('\n')
                line_list = control_line.split('\t')
                query_strand = line_list[0]
                contig_loc = line_list[1]
                start = line_list[2]
                end = line_list[3]
                strand = line_list[4]
                samtools_output_file = f'{control_basename}_samtools.fa'

                if strand == 'plus':
                    start_adjusted = int(start)
                    end_adjusted = int(end)
                    samtools_extract_cmd = f'samtools faidx {control} {contig_loc}:{start_adjusted}-{end_adjusted} > {samtools_output_file}'
                    subprocess.run(samtools_extract_cmd, shell=True)
                else:
                    if strand == 'minus':
                        start_adjusted = int(end)
                        end_adjusted = int(start)
                        samtools_extract_cmd = f'samtools faidx {control} {contig_loc}:{start_adjusted}-{end_adjusted} -i > {samtools_output_file}'
                        subprocess.run(samtools_extract_cmd, shell=True)

                with open(f'{samtools_output_file}', 'r') as samtools_out:
                    header_line = samtools_output_file.split('_')
                    header = '_'.join(header_line[0:2])
                    contents = samtools_out.readlines()
                    seq = []
                    for line in contents[1:]:
                        line = line.rstrip()
                        seq.append(line)
                        full_sequence = ''.join(seq)
                        sequence_dictionary[header] = full_sequence
        else:
            if not check_control_exists:
                print(f'No region found in {control_basename}. Checking next control...')
                continue
    #print(sequence_dictionary)
# Search extracted region for enzymes in enzyme dictionary.
    list_of_enzymes = find_enzymes_with_different_positions(sequence_dictionary, ENZYME_DICT)
    #print(list_of_enzymes)
# Find unique restriction enzyme sites.
    case_basename_list = genome_name.split('_')
    case_basename = '_'.join(case_basename_list[0:2])
    for enzyme in list_of_enzymes:
        if enzyme[case_basename]:
            first_value_second_element = enzyme[case_basename][0][1]
            first_value_third_element = enzyme[case_basename][0][2]
            list_length = len(enzyme)
            target_length = list_length - 1
            i = 0
            j = 0
            for genome, match in enzyme.items():
                if match:
                    second_element = match[0][1] # Starting position
                    third_element = match[0][2]
                    enzyme_match = match[0][0]
                    if second_element == first_value_second_element:
                        pass
                    else:
                        i += 1
                        if i == target_length: # For positions, need to correct for regex 0-based indexing
                            print(f'Unique restriction enzyme site for {enzyme_match} in {case_basename} from positions: {first_value_second_element}-{first_value_third_element}')
                else:
                    if not match:
                        j += 1
                        if j == target_length: # For positions, need to correct for regex 0-based indexing
                            print(f'{enzyme_match} is only found in {case_basename} from positions: {first_value_second_element}-{first_value_third_element}')

        else:
            continue

    # Remove temporary files before proceeding to the next unique kmer
    directory = os.getcwd()
    blast_tmp_files = '*_blast.out'
    samtools_tmp_files = '*_samtools.fa'
    samtools_indices = '*.fai'

    remove_files_with_pattern(directory, blast_tmp_files)
    remove_files_with_pattern(directory, samtools_tmp_files)
    remove_files_with_pattern(extractedFiles_directory, samtools_indices)
    os.remove('uniqueKmer.fa')

# Remove all temporary files, unless otherwise specified.
if ARGS.keep_tmp:
    pass
else:
    os.remove('multiple_datasets.zip')
    os.remove('tmp.meta')
    os.remove('tmp.tmp')
    os.remove('input.list')
    os.remove('output.txt')
    shutil.rmtree('genomes/')


# Print final message. Tell the user what the name of output file.
print(f'Pipeline finished. Results printed to {cured_report}.')
