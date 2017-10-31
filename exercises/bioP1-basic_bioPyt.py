seq=input('Enter DNA or RNA sequence:')

## code
abases=seq.count('a')
cbases=seq.count('c')
gbases=seq.count('g')
tbases=seq.count('t')
ubases=seq.count('u')
seq_length=len(seq)
seq_gc_content=(cbases + gbases) *100 / seq_length
seq_gc_content=round(seq_gc_content,2)
if 'u' in seq :
      print("sequence is a RNA molecule")
else:
      print("sequence is DNA molecule")
print("the length of sequence is",seq_length,"and the GC content is",seq_gc_content,"%")




##### functions


def has_stop_codon(seq):
    stop_codon_found=False
    stop_codons=["tga","tag","taa"]
    for i in range(0,len(seq),3):
        codon=seq[i:i+3].lower()
        if codon in stop_codons:
            stop_codon_found=True
            break
    return stop_codon_found



def complement(seq):
    basecomplement={'A':'T','G':'C','C':'T','T':'A','N':'N','n':'n','a':'t','g':'c','c':'t','t':'a'}
    letters=list(seq)
    letter=[basecomplement[base] for base in letters]
    return ''.join(letters)



################ Install Biopython ##########################

pip install reportlab

testing
>>> from reportlab.graphics import renderPDF

pip install biopython

#### check installed

import Bio
print(Bio.__version__)



################ IUPAC

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


dna=Seq("ATGCTAGCTAGCTAGCTAGCTACAT")

dna.complement()
dna.reverse_complement()


rna=dna.transcribe()
rna

protein=rna.translate()
protein

protein2=dna.translate()   # go from DNA straight to RNA
protein2

################# Codons

from Bio.Data import CodonTable

sorted(CodonTable.unambiguous_dna_by_name["Standard"])


################ SeqRecords Object

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


seq=SeqIO.read('C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\p53.gb','genbank')
print(seq)



seqs=[]
for s in SeqIO.parse('C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\hydra.fasta','fasta'):
   seqs.append(s)

len(seqs)
seqs


### pasrsing through objects one at a time 

seqFH=SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\hydra.fasta","fasta")

seqObj=next(seqFH)
seqObj.seq
seqObj.id
seqObj.description


for seq in seqFH:
    print(seq.id)
