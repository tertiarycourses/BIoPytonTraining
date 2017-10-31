

########################### working with sequences ########################

from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
my_seq

print(my_seq)

my_seq.alphabet

my_seq.complement()

my_seq.reverse_complement()


### FASTA parsing example

from Bio import SeqIO
for seq_record in SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

### GeneBank parsing example

from Bio import SeqIO
for seq_record in SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.gbk", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))




# where possible you should specify the alphabet explicitly when creating 
# your sequence objects - in this case an unambiguous DNA alphabet object:

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
my_seq
my_seq.alphabet



### counting sequence

print(my_seq[0]) #first letter
print(my_seq[2]) #third letter
print(my_seq[-1]) #last letter

"AAAA".count("AA")
Seq("AAAA").count("AA")


#  you may actually want an overlapping count

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC', IUPAC.unambiguous_dna)
len(my_seq)

my_seq.count("G")

100 * float(my_seq.count("G") + my_seq.count("C")) / len(my_seq)


### GC content

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
my_seq = Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC', IUPAC.unambiguous_dna)
GC(my_seq)

### slicing sequence

my_seq[4:12]  # slice of the sequence:

my_seq[1::3]   # first, second and third codon positions of this DNA sequence:

my_seq[::-1]   # reverse the sequence

### adding sequnces (concatanating)

# You may often have many sequences to add together, which can be done with a for loop like this:

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
list_of_seqs = [Seq("ACGT", generic_dna), Seq("AACC", generic_dna), Seq("GGTT", generic_dna)]
concatenated = Seq("", generic_dna)
for s in list_of_seqs:
	     concatenated += s

concatenated

### changing case

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
dna_seq = Seq("acgtACGT", generic_dna)
dna_seq

dna_seq.upper()
dna_seq.lower()


### complement


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
my_seq

my_seq.complement()
my_seq.reverse_complement()


#### Transciptions & translation of sequences

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
coding_dna

messenger_rna = coding_dna.transcribe()
messenger_rna

messenger_rna.translate()

# You can also translate directly from the coding strand DNA sequence:
coding_dna.translate()


### codons

# Suppose we are dealing with a mitochondrial sequence. We need to tell the # translation function to use the relevant genetic code instead:
coding_dna.translate(table="Vertebrate Mitochondrial")


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA" + \
  "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT" + \
  "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT" + \
  "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT" + \
  "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA",
   generic_dna)

gene.translate(table="Bacterial")
gene.translate(table="Bacterial", to_stop=True)

### Translation table

from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

print(standard_table)
print(mito_table)


### comparing seq objects

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
seq1 = Seq("ACGT", IUPAC.unambiguous_dna)
seq2 = Seq("ACGT", IUPAC.ambiguous_dna)
str(seq1) == str(seq2)
str(seq1) == str(seq1)


###  Editing Sequences

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)


# convert it into a mutable sequence (a MutableSeq object) and do pretty much # anything you want with it:

mutable_seq = my_seq.tomutable()
mutable_seq

mutable_seq[5] = "C"
mutable_seq.remove("T")
mutable_seq.reverse()

#  Once you have finished editing your a MutableSeq object, it�s easy to get #  back to a read-only Seq object should you need to:
new_seq = mutable_seq.toseq()
new_seq


