#################### Multiple Sequence Alignment objects #########################

### Parsing or Reading Sequence Alignments

### using ClustalW

from Bio.Align.Applications import ClustalwCommandline
cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
print(cline)

import os
from Bio.Align.Applications import ClustalwCommandline
clustalw_exe = r"C:\Program Files\new clustal\clustalw2.exe"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile="opuntia.fasta")
assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
stdout, stderr = clustalw_cline()

import os
from Bio.Align.Applications import ClustalwCommandline
clustalw_exe = r"C:\Program Files\new clustal\clustalw2.exe"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile="opuntia.fasta")
assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
stdout, stderr = clustalw_cline()

from Bio import Phylo
tree = Phylo.read("opuntia.dnd", "newick")
Phylo.draw_ascii(tree)


### using  MUSCLE

from Bio.Align.Applications import MuscleCommandline

from Bio.Align.Applications import MuscleCommandline
cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.txt")

muscle -in opuntia.fasta -out opuntia.txt


from Bio.Align.Applications import MuscleCommandline
muscle_cline = MuscleCommandline(input="opuntia.fasta")
stdout, stderr = muscle_cline()
from StringIO import StringIO
from Bio import AlignIO
align = AlignIO.read(StringIO(stdout), "fasta")
print(align)


import subprocess
from Bio.Align.Applications import MuscleCommandline
muscle_cline = MuscleCommandline(input="opuntia.fasta")
child = subprocess.Popen(str(muscle_cline),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True,
                             shell=(sys.platform!="win32"))
from Bio import AlignIO
align = AlignIO.read(child.stdout, "fasta")
print(align)




#### Biopython�s pairwise2

# global pairwise alignment between the same two hemoglobin sequences from above 
# (HBA_HUMAN, HBB_HUMAN) stored in alpha.faa and beta.faa:

from Bio import pairwise2
from Bio import SeqIO
seq1 = SeqIO.read("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\alpha.faa", "fasta")
seq2 = SeqIO.read("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\beta.faa", "fasta")
alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)


len(alignments)

print(pairwise2.format_alignment(*alignments[0]))



# Better alignments are usually obtained by penalizing gaps: higher costs for opening a gap and 
# lower costs for extending an existing gap.

from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
vseq1 = SeqIO.read("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\alpha.faa", "fasta")
seq2 = SeqIO.read("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\beta.faa", "fasta")
alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
len(alignments)

print(pairwise2.format_alignment(*alignments[0]))

