############################# Graphics including GenomeDiagram ##########################

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
record = SeqIO.read("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\NC_005816.gb", "genbank")

gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

for feature in record.features:
    if feature.type != "gene":
        #Exclude this feature
        continue
    if len(gd_feature_set) % 2 == 0:
        color = colors.blue
    else:
        color = colors.lightblue
    gd_feature_set.add_feature(feature, color=color, label=True)

gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4',
                fragments=4, start=0, end=len(record))
gd_diagram.write("plasmid_linear.pdf", "PDF")
gd_diagram.write("plasmid_linear.eps", "EPS")
gd_diagram.write("plasmid_linear.svg", "SVG")
gd_diagram.write("plasmid_linear.png", "PNG")


# If you want to do a circular figure, then try this:

gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
                start=0, end=len(record), circle_core=0.7)
gd_diagram.write("plasmid_circular.pdf", "PDF")


#### arrow sigils

#Full height shafts, giving pointed boxes:
gd_feature_set.add_feature(feature, sigil="ARROW", color="brown",
                           arrowshaft_height=1.0)


### example > including RE sites

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

record = SeqIO.read("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\NC_005816.gb", "genbank")

gd_diagram = GenomeDiagram.Diagram(record.id)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

for feature in record.features:
    if feature.type != "gene":
        #Exclude this feature
        continue
    if len(gd_feature_set) % 2 == 0:
        color = colors.blue
    else:
        color = colors.lightblue
    gd_feature_set.add_feature(feature, sigil="ARROW",
                               color=color, label=True,
                               label_size = 14, label_angle=0)

#I want to include some strandless features, so for an example
#will use EcoRI recognition sites etc.
for site, name, color in [("GAATTC","EcoRI",colors.green),
                          ("CCCGGG","SmaI",colors.orange),
                          ("AAGCTT","HindIII",colors.red),
                          ("GGATCC","BamHI",colors.purple)]:
    index = 0
    while True:
        index  = record.seq.find(site, start=index)
        if index == -1 : break
        feature = SeqFeature(FeatureLocation(index, index+len(site)))
        gd_feature_set.add_feature(feature, color=color, name=name,
                                   label=True, label_size = 10,
                                   label_color=color)
        index += len(site)

gd_diagram.draw(format="linear", pagesize='A4', fragments=4,
                start=0, end=len(record))
gd_diagram.write("plasmid_linear_nice.pdf", "PDF")
gd_diagram.write("plasmid_linear_nice.eps", "EPS")
gd_diagram.write("plasmid_linear_nice.svg", "SVG")

gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
                start=0, end=len(record), circle_core = 0.5)
gd_diagram.write("plasmid_circular_nice.pdf", "PDF")
gd_diagram.write("plasmid_circular_nice.eps", "EPS")
gd_diagram.write("plasmid_circular_nice.svg", "SVG")



################# Histograms of sequence lengths




from Bio import SeqIO
sizes = [len(rec) for rec in SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\hydra.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)
sizes

import pylab
pylab.hist(sizes, bins=20)

#pylab.title("%i orchid sequences\nLengths %i to %i" \
#            % (len(sizes),min(sizes),max(sizes)))

pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.show()


####################### Plots of sequence by GC

from Bio import SeqIO
from Bio.SeqUtils import GC

gc_values = sorted(GC(rec.seq) 
for rec in SeqIO.parse("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\hydra.fasta", "fasta"))

import pylab
pylab.plot(gc_values)

#pylab.title("%i orchid sequences\nGC%% %0.1f to %0.1f" \
#            % (len(gc_values),min(gc_values),max(gc_values)))

pylab.xlabel("Genes")
pylab.ylabel("GC%")
#pylab.show()
