
########################### Sequence Input/Output ##############

### Reading Sequence Files

from Bio import SeqIO
for seq_record in SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.fasta", "fasta"):   # fasta file
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


from Bio import SeqIO
for seq_record in SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.gbk", "genbank"):   # geneBank file
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


### Iterating over the records in a sequence file


# Instead of using a for loop, can also use the next() function on an 
# iterator to step through the entries, like this:

from Bio import SeqIO
record_iterator = SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.fasta", "fasta")

first_record = next(record_iterator)
print(first_record.id)
print(first_record.description)

second_record = next(record_iterator)
print(second_record.id)
print(second_record.description)

### Getting a list of the records in a sequence file

from Bio import SeqIO
records = list(SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.gbk", "genbank"))

print("Found %i records" % len(records))

print("The last record")
last_record = records[-1] #using Python's list tricks
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))

print("The first record")
first_record = records[0] #remember, Python counts from zero
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))

###  Extracting data

from Bio import SeqIO
record_iterator = SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.gbk", "genbank")
first_record = next(record_iterator)
print(first_record)


print(first_record.annotations)
print(first_record.annotations.keys())
print(first_record.annotations.values())
print(first_record.annotations["source"])
print(first_record.annotations["organism"])


# list of the species each orchid sequence is from:

from Bio import SeqIO
all_species = []
for seq_record in SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.gbk", "genbank"):
    all_species.append(seq_record.annotations["organism"])
print(all_species)


### Parsing sequences from compressed files

from Bio import SeqIO
print(sum(len(r) for r in SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.gbk", "gb")))


from Bio import SeqIO
with open("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ls_orchid.gbk") as handle:
        print(sum(len(r) for r in SeqIO.parse(handle, "gb")))

# dealing with gzip compressed file 

import gzip               OR  import bz2 (for bzip2 compressed files)
from Bio import SeqIO
with gzip.open("ls_orchid.gbk.gz", "rt") as handle:
     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))


#### Parsing sequences over Intenet > online databases

### Entrez

from Bio import Entrez
from Bio import SeqIO
Entrez.email = "A.N.Other@example.com"
with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="6273291") as handle:
    seq_record = SeqIO.read(handle, "fasta")
print("%s with %i features" % (seq_record.id, len(seq_record.features)))



from Bio import Entrez
from Bio import SeqIO
Entrez.email = "A.N.Other@example.com"
with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="6273291") as handle:
    seq_record = SeqIO.read(handle, "gb") #using "gb" as an alias for "genbank"
print("%s with %i features" % (seq_record.id, len(seq_record.features)))


# several records. This time the handle contains multiple records, so we # must use the Bio.SeqIO.parse() function:

from Bio import Entrez
from Bio import SeqIO
Entrez.email = "A.N.Other@example.com"
with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="6273291,6273290,6273289") as handle:
    for seq_record in SeqIO.parse(handle, "gb"):
        print("%s %s..." % (seq_record.id, seq_record.description[:50]))
        print("Sequence length %i, %i features, from: %s"
              % (len(seq_record), len(seq_record.features), seq_record.annotations["source"]))

### SwissProt

from Bio import ExPASy
from Bio import SeqIO
with ExPASy.get_sprot_raw("O23729") as handle:
    seq_record = SeqIO.read(handle, "swiss")
print(seq_record.id)
print(seq_record.name)
print(seq_record.description)
print(repr(seq_record.seq))
print("Length %i" % len(seq_record))
print(seq_record.annotations["keywords"])


#### Writing Sequence Files

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

rec1 = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD" \
                    +"GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK" \
                    +"NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM" \
                    +"SSAC", generic_protein),
                 id="gi|14150838|gb|AAK54648.1|AF376133_1",
                 description="chalcone synthase [Cucumis sativus]")

rec2 = SeqRecord(Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ" \
                    +"DMVVVEIPKLGKEAAVKAIKEWGQ", generic_protein),
                 id="gi|13919613|gb|AAK33142.1|",
                 description="chalcone synthase [Fragaria vesca subsp. bracteata]")

rec3 = SeqRecord(Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC" \
                    +"EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP" \
                    +"KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN" \
                    +"NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV" \
                    +"SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW" \
                    +"IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT" \
                    +"TGEGLEWGVLFGFGPGLTVETVVLHSVAT", generic_protein),
                 id="gi|13925890|gb|AAK49457.1|",
                 description="chalcone synthase [Nicotiana tabacum]")

my_records = [rec1, rec2, rec3]

from Bio import SeqIO
SeqIO.write(my_records, "my_example.faa", "fasta")   # write into fasta format


# converting formats

from Bio import SeqIO
count = SeqIO.convert("ls_orchid.gbk", "genbank", "my_example.fasta", "fasta")


### Converting a file of sequences to their reverse complements

from Bio import SeqIO
for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
     print(record.id)
     print(record.seq.reverse_complement())






