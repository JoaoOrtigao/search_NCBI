#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Created on Wed Nov 27 2019

@author: joao ortigao

"""

from datetime import datetime
from os.path import join as pathjoin
from os import getcwd as osgetcwd
from  Bio import Entrez

##############################################################################

import argparse

parser = argparse.ArgumentParser(
    description='''Seacrh for pacient data at NCBI databases''',
    epilog='''Search_NCBI_v2.0 27-11-2019''')

parser.add_argument('--o', required=True,
                    help='OUTDIR: output directory')

parser.add_argument('--l', required=True,
                    help='DISEASE_LIST.txt: list containing 1 disease per line')

parser.add_argument('--e', required=True,
                    help='user@email.com: email adress to use as Entrez identification')

parser.add_argument('--version', action='version', version='%(prog)s 2.0')


args=parser.parse_args()

##############################################################################

print("\n\nSTART TIME:\t",datetime.now())

try:
    OUTDIR = args.o
except:
    OUTDIR = 'OUTDIR'

print("\n\noutput directory is:",pathjoin(osgetcwd(),OUTDIR))
    
try:
    DISEASE_LIST = args.l
    print("\n\nDISEASE_LIST is:",pathjoin(osgetcwd(),DISEASE_LIST))
except:
    print("\n\nNo DISEASE_LIST found. Exiting...\n\n")
    exit()
    
try:
    Entrez.email = args.e
    print("\n\nEntrez.email is: ",Entrez.email,"\n\n")
except:
    print("\n\nNo valid Entrez.email was declared. Exiting")
    exit()

##############################################################################
    # MAKE SEARCH AND DOWNLOAD

from Search_NCBI_v02_functions import efetch_found_bioprojects 
from Search_NCBI_v02_functions import esearch_disease 
efetch_found_bioprojects(esearch_disease(DISEASE_LIST,OUTDIR))

##############################################################################
    # CREATE db_para_curagem

from Search_NCBI_v02_functions import collect_XML_file
from Search_NCBI_v02_functions import classify_disease
from Search_NCBI_v02_functions import writer
df=collect_XML_file(OUTDIR)
df2=classify_disease(df,OUTDIR,DISEASE_LIST)
writer(df2,OUTDIR)

##############################################################################

print("END TIME:\t",datetime.now())


