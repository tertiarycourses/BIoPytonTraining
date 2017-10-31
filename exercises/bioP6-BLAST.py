############################### BLAST sequencing ###################################

# if you know the nucelotide sequence in database number

from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")

# blasting a fasta file

from Bio.Blast import NCBIWWW
fasta_string = open("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\AK292608.fasta").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)


# read in the FASTA file as a SeqRecord

from Bio.Blast import NCBIWWW
from Bio import SeqIO
record = SeqIO.read("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\AK292608.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

#if we wish to include the existing identifier
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))


### Parsing BLAST output

from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)
len(blast_record.alignments)


# or, if you have lots of results (i.e., multiple query sequences):
from Bio.Blast import NCBIXML
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)

# Hit objects represent all query results from a single database entry
# HSP (high-scoring pair) represents region(s) in the hit sequence that contains significant # alignment(s) to the query sequence.
# HSPFragment represents a single, contiguous match between the query and hit sequences.



###################### NCBI blast

## we are given an unknown sequence > find out which species it is from

unknown_gene="CATGCTACGGTGCTAAAAGCATTACGCCCTATAGTGATTTTCGAGACATACTGTGTTTTTAAATATAGTATTGCC"


resultHandle=NCBIWWW.qblast("blastn","nt",unknown_gene)       # program, database, sequence
blastRecord=NCBIXML.read(resultHandle)                        # returns XML object
len(blastRecord.alignments)

E_VALUE_THRESH = 0.01
for alignment in blastRecord.alignments:
    for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                    print('****Alignment****')
                    print('sequence:', alignment.title)
                    print('length:', alignment.length)
                    print('e value:', hsp.expect)
                    print(hsp.query)
                    print(hsp.match)
                    print(hsp.sbjct)

import numpy
import pandas as pd

pd.DataFrame({'a':1, 'b':2}, index=[0])   # simple pandas data frame > index[0] means one row




### create a file

dataArray=[]

E_VALUE_THRESH = 0.01
for alignment in blastRecord.alignments:
    for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                   df= pd.DataFrame({'sequence':alignment.title,
                                 'length':alignment.length,
                                  'e_value':hsp.expect}, index=[0])
                   dataArray.append(df)


dataArray[0]  # see the first value

pd.concat(dataArray, ignore_index=True)

dataArray.head()
dataArray.shape()



