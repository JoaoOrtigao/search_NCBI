#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Nov 21 2019

@author: joao ortigao
"""

from pandas import DataFrame as pd_DataFrame
from pandas import ExcelWriter as pd_ExcelWriter
from os.path import join as pathjoin
from os.path import exists as pathexists
from os.path import isfile 
from os import mkdir as osmkdir
from os import getcwd as osgetcwd
from  Bio import Entrez
from re import sub
from glob import glob as gb
import xml.etree.ElementTree as ET
parser = ET.XMLParser(encoding="utf-8")

##############################################################################

def CREATE_DIR(OUTDIR):
    if not pathexists(pathjoin(OUTDIR)):
        osmkdir(pathjoin(OUTDIR))
        
    for DIR in ["IdListDIR","IdListDIR/disease","IdListDIR/query"]:
        if not pathexists(pathjoin(OUTDIR,DIR)):
            osmkdir(pathjoin(OUTDIR,DIR))

##############################################################################
    
def MAKE_DICIONARY(DISEASE_LIST):
    DISEASES=[]
    DISEASES = [line.rstrip('\n') for line in open(DISEASE_LIST)]
    CODES = [s.replace(' ', '_') for s in DISEASES]
    CODES = [s.replace('\'', '') for s in CODES]
    DIC = dict(zip(DISEASES,CODES))
    return(DIC)
    
##############################################################################
    
def esearch_disease(DISEASE_LIST,OUTDIR):

    CREATE_DIR(OUTDIR)
    
    DISEASE_DIC = MAKE_DICIONARY(DISEASE_LIST)
    
    # data frame to store all Counts
    # +2 for one extra line for "COUNTS" and "TOTAL1"
    df=pd_DataFrame(index=range(0,len(DISEASE_DIC)+2),columns=range(0,8)) 
    df.columns=["disease","COD","QUERY1","QUERY2","QUERY3","QUERY4",\
                "QUERY5","TOTAL2"]
    COL1=list(DISEASE_DIC); COL1.append('COUNTS'); COL1.append('TOTAL1')
    df['disease']=COL1
    
    # data frame to store all the commands used for each search
    COMMAND=pd_DataFrame(index=range(0,len(DISEASE_DIC)),columns=range(0,8))
    COMMAND.columns=["disease","COD","QUERY1","QUERY2","QUERY3","QUERY4",\
                     "QUERY5","END"]
    COMMAND["disease"]=COL1[0:len(DISEASE_DIC)]
    COMMAND["END"]='.'
    
    # data frameto store the queries' explanations
    QUERY_description=pd_DataFrame(index=range(0,5),columns=range(0,1))
    QUERY_description.columns=["DESCRIPTION"]
    QUERY_description.index=["QUERY1","QUERY2","QUERY3","QUERY4","QUERY5"]
    QUERY1_desc='Procura o nome da doença em todos os campos e filtra por'\
                ' experimentos de expressão gênica feitos com amostras '\
                'humanas. Essa é a QUERY mais abrangente.'
    QUERY2_desc='Igual a QUERY1 só que também procura por "patient" OU '\
                '"patients" em todos os campos'
    QUERY3_desc='Igual a QUERY2 só que também filtra por bioprojects '\
                'presentes na base de dados SRA'
    QUERY4_desc='Procura o nome da doença somente no título do bioproject, '\
                'procura por "patient" OU "patients" em todos os campos e '\
                'filtra por experimentos de expressão gênica feitos com '\
                'amostras humanas'
    QUERY5_desc='Igual a QUERY4 só que também filtra por bioprojects '\
                'presentes na base de dados SRA'
    QUERY_description["DESCRIPTION"]=[QUERY1_desc,QUERY2_desc,QUERY3_desc,\
                                      QUERY4_desc,QUERY5_desc]
    
    IdList_QUERY1=[]
    IdList_QUERY2=[]
    IdList_QUERY3=[]
    IdList_QUERY4=[]
    IdList_QUERY5=[]
    IdList_total=[]

    N=0
    for DISEASE in list(DISEASE_DIC):
        
        print(str(N)+'\t'+DISEASE)
            
        COD=DISEASE_DIC[DISEASE]
        df["COD"][N]=COD
        COMMAND["COD"][N]=COD
         
        QUERY_DIC={'1':'("'+DISEASE+'"[All Fields])AND'\
                        '("transcriptome gene expression"[Filter]AND"org '\
                        'human"[Filter])',
                   '2':'("'+DISEASE+'"[All Fields]AND'\
                        '("patient"[All Fields]OR"patients"[All Fields])AND'\
                        '("transcriptome gene expression"[Filter]AND"org '\
                        'human"[Filter])',
                   '3':'("'+DISEASE+'"[All Fields]AND'\
                        '("patient"[All Fields]OR"patients"[All Fields])AND'\
                        '("transcriptome gene expression"[Filter]AND"org '\
                        'human"[Filter]AND"bioproject sra"[Filter])',
                   '4':'("'+DISEASE+'"[Title]AND'\
                        '("patient"[All Fields]OR"patients"[All Fields])AND'\
                        '("transcriptome gene expression"[Filter]AND"org '\
                        'human"[Filter])',
                   '5':'("'+DISEASE+'"[Title]AND'\
                        '("patient"[All Fields]OR"patients"[All Fields])AND'\
                        '("transcriptome gene expression"[Filter]AND"org '\
                        'human"[Filter])AND"bioproject sra"[Filter])'}

        Idlist_disease=[]

        ROUND=['1','2','3','4','5']
        for R in ROUND:
            QUERY='QUERY'+R
            TERM=QUERY_DIC[R]
#            COMMAND[locals[QUERY]][N]=TERM
            
            handle = Entrez.esearch(db="bioproject", retmax=1000, 
                                        term=TERM)
            record = Entrez.read(handle)
            handle.close()   
            if int(record["Count"]) > 1000:
                print('\nATTENTION!\nn'+record["Count"]+' bioprojects are '\
                      'related to this esearch and only 1000 will be written '\
                      'to the Idlist for the further analysis.\n\n'+QUERY+\
                      'for '+DISEASE+'\n\n'+QUERY_DIC[R]+'\n')
                exit

            # MONTAR LISTA POR DOENÇA
            Idlist_disease+=list(record["IdList"])
            IdList_total+=list(record["IdList"])

            # ADD IDS TO QUERY AND TOTAL LISTS
#            IdList_total+=record["IdList"]
            if R == '1':
                IdList_QUERY1+=list(record["IdList"])
                COMMAND['QUERY1'][N]=TERM
                df['QUERY1'][N]=int(record["Count"])
            elif R == '2':
                IdList_QUERY2+=list(record["IdList"])
                COMMAND['QUERY2'][N]=TERM
                df['QUERY2'][N]=int(record["Count"])
            elif R == '3':
                IdList_QUERY3+=list(record["IdList"])
                COMMAND['QUERY3'][N]=TERM
                df['QUERY3'][N]=int(record["Count"])
            elif R == '4':
                IdList_QUERY4+=list(record["IdList"])
                COMMAND['QUERY4'][N]=TERM
                df['QUERY4'][N]=int(record["Count"])
            elif R == '5':
                IdList_QUERY5+=list(record["IdList"])
                COMMAND['QUERY5'][N]=TERM
                df['QUERY5'][N]=int(record["Count"])

        #remove replicates from the list      
        Idlist_disease=list(set(Idlist_disease)) 
        
        df['TOTAL2'][N]=len(Idlist_disease)
        
        outfile=pathjoin(OUTDIR,"IdListDIR/disease",COD+".txt")
        with open(outfile, 'w') as f:
            print( "\n".join(Idlist_disease), file = f)
            f.close()

        N+=1

    #preencher a linha com totais
    for COL in list(df)[2:len(df)]: #COL da terceira coluna até a última
        df[COL][len(DISEASE_DIC)]=df[COL][0:len(DISEASE_DIC)].sum(axis=0)
    
    # ESCREVER DEMAIS LISTAS PARA ARQUIVOS TXT
    IdList_total=list(set(IdList_total))
    outfile=pathjoin(OUTDIR,"IdListDIR/IdList_total.txt")
    with open(outfile, 'w') as f:
        print( "\n".join(IdList_total), file = f)
        f.close()

    IdList_QUERY1=list(set(IdList_QUERY1))
    df.loc[len(DISEASE_DIC)+1,"QUERY1"] = len(IdList_QUERY1)
    outfile=pathjoin(OUTDIR,"IdListDIR/query","IdList_QUERY1.txt")
    with open(outfile, 'w') as f:
        print( "\n".join(IdList_QUERY1), file = f)
        f.close()
    
    
    IdList_QUERY2=list(set(IdList_QUERY2))
    df.loc[len(DISEASE_DIC)+1,"QUERY2"] = len(IdList_QUERY2)
    outfile=pathjoin(OUTDIR,"IdListDIR/query","IdList_QUERY2.txt")
    with open(outfile, 'w') as f:
        print( "\n".join(IdList_QUERY2), file = f)
        f.close()
    
    
    IdList_QUERY3=list(set(IdList_QUERY3))
    df.loc[len(DISEASE_DIC)+1,"QUERY3"] = len(IdList_QUERY3)
    outfile=pathjoin(OUTDIR,"IdListDIR/query","IdList_QUERY3.txt")
    with open(outfile, 'w') as f:
        print( "\n".join(IdList_QUERY3), file = f)
        f.close()
        
        
    IdList_QUERY4=list(set(IdList_QUERY4))
    df.loc[len(DISEASE_DIC)+1,"QUERY4"] = len(IdList_QUERY4)
    outfile=pathjoin(OUTDIR,"IdListDIR/query","IdList_QUERY4.txt")
    with open(outfile, 'w') as f:
        print( "\n".join(IdList_QUERY4), file = f)
        f.close()
    
    IdList_QUERY5=list(set(IdList_QUERY5))
    df.loc[len(DISEASE_DIC)+1,"QUERY5"] = len(IdList_QUERY5)
    outfile=pathjoin(OUTDIR,"IdListDIR/query","IdList_QUERY5.txt")
    with open(outfile, 'w') as f:
        print( "\n".join(IdList_QUERY5), file = f)
        f.close()

    #ESCREVER TODOS OS RESULTADOS PARA UM ARQUIVO EXCEL
    writer = pd_ExcelWriter(pathjoin(OUTDIR,'search_NCBI_RESULT.xlsx'), 
                            engine='xlsxwriter')
    df.to_excel(writer, sheet_name='counts')
    COMMAND.to_excel(writer, sheet_name='command_lines')
    QUERY_description.to_excel(writer, sheet_name='query_description')
    writer.save()

    return(pathjoin(osgetcwd(),OUTDIR))
    
##############################################################################    

def efetch_found_bioprojects(OUTDIR):

    
    def printProgressBar (iteration, total, prefix = '', suffix = '', \
                      decimals = 1, length = 100, fill = '█'):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent \
                                      complete (Int)
                                      length      - Optional  : character length of bar (Int)
                                      fill        - Optional  : bar fill character (Str)
        """
        percent = ("{0:." + str(decimals) + "f}")\
                                    .format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
        # Print New Line on Complete
        if iteration == total: 
            print()
    
    
    """
    COLETAR INFORMAÇOES SOBRE BIOPROJECTS ECONTRADOS
    """
    
    if pathexists(OUTDIR):
        
        for DIR in ['Bioprojects','Bioprojects/xml']:
            if not pathexists(pathjoin(OUTDIR,DIR)):
                osmkdir(pathjoin(OUTDIR,DIR))
        
        path_to_list=pathjoin(OUTDIR,'IdListDIR/IdList_total.txt')
        if isfile(path_to_list):
            with open(path_to_list,'r') as f:
                IdList_total=list(filter(None, f.read().splitlines()))
        else:
            print('File '+f+' was not found. Run esearch_disease(OUTDIR) '\
                  'for making it.')
            exit()
    else:
        print('Directory '+pathjoin(OUTDIR)+' is not accessible. Did you run'\
              'esearch_disease() previously? If not, do it and try again.')
        exit()
    
    df2=pd_DataFrame(index=range(0,len(IdList_total)),columns=range(0,7)) 
    df2.columns=["ID","accession","GEO","title","abstract","disease","COD"]
    df2["ID"]=IdList_total
    
    print("\n\n") # ESSE PRINT SERVE PARA DISTANCIAR A BARRA DE PROCESSAMENTO 
                  # QUE VEM LOGO ABAIXO DENTRO DO LOOPING

    # prepare bar progress                  
    l = len(IdList_total)
    i=0
    printProgressBar(0, l, prefix = 'Download:', suffix = 'Complete', 
                     length = 50)    

    RECALL=[] # if download fails, the ID is stored in RECALL
    DIC_ID={}
    for ID in IdList_total:

        try:
            handle = Entrez.efetch(db="bioproject", id=ID)
        except:
            RECALL+=[ID]
            print('handle = Entrez.efetch(db="bioproject", id='+ID+')\tFAILED')
            continue # avoid catastrophic event in case NCBI fails to give 
                     # the informatio for one ID
        try:    
            record = handle.read()
            root = ET.fromstring(record)
            DIC = root.find(".//ProjectID/ArchiveID").attrib
            DIC_ID[DIC['accession']] = DIC_ID.get(DIC['accession'],DIC['id'])
    
            outfile=pathjoin(OUTDIR,'Bioprojects/xml',DIC['accession']+\
                                 '_'+DIC['id']+'.xml')
    
            #print(outfile)
            with open(outfile, "w", encoding="utf-8") as f:
                print(record, file = f)

        except:
            RECALL+=[ID]
            print('FAILED to process '+ID+' during the first trial')
            continue    

        printProgressBar(i+1, l, prefix = 'Download:', suffix = 'Complete', 
                         length = 50)
        i+=1

    # RECALL for failure IDs
    if len(RECALL) > 0:
        print("\n\nFailure to download IDs. STARTING RECALL.")

        l = len(RECALL)
        i=0
        printProgressBar(0, l, prefix = 'Download:', suffix = 'Complete', 
                         length = 50)
        RECALL2=[]
        for ID in RECALL:
    
            try:
                handle = Entrez.efetch(db="bioproject", id=ID)
            except:
                RECALL2+=[ID]
                print('handle = Entrez.efetch(db="bioproject", id='+ID+')'\
                      '\tFAILED in RECALL')
                continue 
            try:                             
                record = handle.read()
                root = ET.fromstring(record)
                DIC = root.find(".//ProjectID/ArchiveID").attrib
                DIC_ID[DIC['accession']] = DIC_ID.get(DIC['accession'],DIC['id'])
        
                outfile=pathjoin(OUTDIR,'Bioprojects/xml',DIC['accession']+\
                                     '_'+DIC['id']+'.xml')
        
                #print(outfile)
                with open(outfile, "w", encoding="utf-8") as f:
                    print(record, file = f)
            except:
                RECALL2+=[ID]
                print('FAILED to process '+ID+' during the RECALL')
                continue
    
            printProgressBar(i+1, l, prefix = 'RECALL:', suffix = 'Complete', 
                             length = 50)
            i+=1
    
        if len(RECALL2) > 0:
            outfile=pathjoin(OUTDIR,'Bioprojects/','RECALL_failure.txt')    
            open(outfile,'w').write(str(RECALL2))
            print("It was not possible to get ID even during the RECALL\nYou"\
                  "can find the problematic IDs on file:\n"+outfile)


    outfile=pathjoin(OUTDIR,'Bioprojects/','dict_ID_ACC.txt')    
    open(outfile,'w').write(str(DIC_ID))
    
##############################################################################
    
def collect_XML_file(OUTDIR):
    # aqui são importados os xml com a descrição de cada bioproject
    files = gb(pathjoin(OUTDIR,'Bioprojects/xml/*.xml'))
    
    df=pd_DataFrame(index=range(0,len(files)),columns=range(0,13)) 
    df.columns=["ID","accession","GEO","title","abstract","disease1",\
                "COD1","disease2","COD2","disease3","COD3","disease4","COD4"]
    
    DIC_ID={}
    N=0
    for file in files:
        #with open(file, "r", encoding="utf-8") as f:
            #contents = f.read()
            #tree = ET.fromstring(contents)

        tree = ET.parse(file,parser=ET.XMLParser(encoding="utf-8"))
        root = tree.getroot()
        
        try:
            GEO = root.find(".//ExternalLink/dbXREF/ID").text
        except:    
            GEO = None # declare empty variable
        title = root.find(".//ProjectDescr/Title").text
        abstract = root.find(".//ProjectDescr/Description").text
        
        DIC = root.find(".//ProjectID/ArchiveID").attrib
        DIC_ID[DIC['accession']] = DIC_ID.get(DIC['accession'],DIC['id'])
        accession=DIC['accession']
        ID=DIC['id']
    
        for COL in ['ID','accession','GEO','title','abstract']:
            df[COL][N]=locals()[COL]
            #print(N)
        
        N+=1
        
    return df

##############################################################################

def classify_disease(df2,OUTDIR,DISEASE_LIST):
    
    DATADIR=OUTDIR
    
    COD_DIC = MAKE_DICIONARY(DISEASE_LIST)
    COD_DIC = {v: k for k, v in COD_DIC.items()} # invert the dictionary map
    
    files2 = gb(pathjoin(OUTDIR,'IdListDIR/disease/*.txt'))
    
    DISEASE1=[]
    DISEASE2=[]
    DISEASE3=[]
    DISEASE4=[]
    
    for file2 in files2:
        
        COD = sub('.txt','',sub('.*IdListDIR/disease\\\\','',file2))
        DISEASE = COD_DIC[COD]
    
        with open(file2,'r') as f:
            IDs = filter(None, f.read().splitlines())
            f.close()
    
        ROUND=['1','2','3','4']
        for ID in IDs:        
            #print(ID)
            for R in ROUND:
                if ID not in locals()["DISEASE"+R]:
                    POS=df2[df2["ID"] == ID].index[0]
                    df2.loc[[POS],'disease'+R] = DISEASE
                    df2.loc[[POS],'COD'+R] = COD
                    locals()["DISEASE"+R].append(ID)
                    break
    return df2

##############################################################################

def writer(df2, OUTDIR):
    writer = pd_ExcelWriter(pathjoin(OUTDIR,'db_para_curagem.xlsx'), 
                            engine='xlsxwriter')
    df2.to_excel(writer, sheet_name='db_completo_nov18')
    writer.save()

    df2.to_csv(pathjoin(OUTDIR,'db_para_curagem.tsv'),sep='\t')
