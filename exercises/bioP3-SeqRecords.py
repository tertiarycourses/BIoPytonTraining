##########################  Creating a SeqRecord  #########################


### SeqRecord objects from scratch

# To create a SeqRecord at a minimum you just need a Seq object:

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

simple_seq = Seq("GATC")
simple_seq_r.id = "AC12345"
simple_seq_r.description = "Made up sequence I wish I could write a paper about"

simple_seq_r = SeqRecord(simple_seq, id="AC12345")


### SeqRecord objects from FASTA files

from Bio import SeqIO
record = SeqIO.read("NC_005816.fna", "fasta")
record.seq
record.id
record.name
record.description

### SeqRecord objects from GenBank files

from Bio import SeqIO
record = SeqIO.read("NC_005816.gb", "genbank")
record.seq
record.id
record.name
record.description

### Sequence described by a feature or location

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
example_parent = Seq("ACCGAGACGGCAAAGGCTAGCATAGGTATGAGACTTCCTTCCTGCCAGTGCTGAGGAACTGGGAGCCTAC")
example_feature = SeqFeature(FeatureLocation(5, 18), type="gene", strand=-1)


feature_seq = example_parent[example_feature.location.start:example_feature.location.end].reverse_complement()
print(feature_seq)

### comapre sequences

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
record1 = SeqRecord(Seq("ACGT"), id="test")
record2 = SeqRecord(Seq("ACGT"), id="test")

record1 == record2


### Slicing a SeqRecord

from Bio import SeqIO
record = SeqIO.read("NC_005816.gb", "genbank")
len(record)
len(record.features)

print(record.features[20])  # focus in on the pim gene, YP_pPCP05.


sub_record = record[4300:4800]  # slice this parent record from 4300 to 4800
sub_record

print(sub_record.features[0])
print(sub_record.features[1])

sub_record.id
sub_record.name
sub_record.description


#### Adding SeqRecord objects

from Bio import SeqIO
record = next(SeqIO.parse("example.fastq", "fastq"))
len(record)
25
print(record.seq)

left = record[:20]
print(left.seq)

right = record[21:]
print(right.seq)

edited = left + right  # add the two parts together

edited = record[:20] + record[21:]  #You can make this shorter with just:









