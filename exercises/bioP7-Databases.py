############################### ENTREZ #############################

from Bio import Entrez
Entrez.email = "A.N.Other@example.com"    # supply your own NCBI email id


handle = Entrez.einfo()
record = Entrez.read(handle)

record.keys()
record["DbList"]

### Searching the Entrez databases using ESearch

# search in PubMed about BioPython

handle = Entrez.esearch(db="pubmed", term="biopython")
record = Entrez.read(handle)
print(record["IdList"])


# search GenBank

handle = Entrez.esearch(db="nucleotide", term="Cypripedioideae[Orgn] AND matK[Gene]", idtype="acc")
record = Entrez.read(handle)
record["Count"]
record["IdList"]

# list of computational journal titles

handle = Entrez.esearch(db="nlmcatalog", term="computational[Journal]", retmax='20')
record = Entrez.read(handle)
print("{} computational journals found".format(record["Count"]))


###  ESummary: Retrieving summaries from primary IDs

handle = Entrez.esummary(db="nlmcatalog", id="101660833")
record = Entrez.read(handle)
info = record[0]['TitleMainList'][0]
print("Journal info\nid: {}\nTitle: {}".format(record[0]["Id"], info["Title"]))


### EFetch: Downloading full records from Entrez

# downloading sequences in the FASTA or GenBank/GenPept plain text formats  downloading sequences 
# in the FASTA or GenBank/GenPept plain text formats 

handle = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")
print(handle.read())


### ELink: Searching for related items in NCBI Entrez

pmid = "19304878"
record = Entrez.read(Entrez.elink(dbfrom="pubmed", id=pmid))

record[0]["DbFrom"]
record[0]["IdList"]

len(record[0]["LinkSetDb"])


record[0]["LinkSetDb"][0]["Link"][0]
record[0]["LinkSetDb"][0]["Link"][1]

####### parsing Medline records 

handle = Entrez.esearch(db="pubmed", term="biopython")
record = Entrez.read(handle)
record["IdList"]
idlist = record["IdList"]
handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")

###### parsing GEO records

handle = Entrez.esearch(db="gds", term="GSE16")
record = Entrez.read(handle)
handle.close()
record["Count"]
record["IdList"]


### examples

from Bio import Entrez
Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
handle = Entrez.egquery(term="orchid")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
     if row["DbName"]=="pubmed":
         print(row["Count"])

handle = Entrez.esearch(db="pubmed", term="orchid", retmax=1678)   # assuming we have 463 orchird articles
record = Entrez.read(handle)
handle.close()
idlist = record["IdList"]

print(idlist)


from Bio import Medline
handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
records = Medline.parse(handle)


###### parsing Entrez Nucleotide records

handle = Entrez.egquery(term="Cypripedioideae")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
     if row["DbName"]=="nuccore":
        print(row["Count"])


from Bio import Entrez
handle = Entrez.esearch(db="nucleotide", term="Cypripedioideae", retmax=814, idtype="acc")
record = Entrez.read(handle)
handle.close()

print(record.keys())

print(record["Count"])

len(record["IdList"])
record["IdList"][:5]

#print(record[0].keys())
#print(record[0]["GBSeq_primary-accession"])
#print(record[0]["GBSeq_other-seqids"])
#print(record[0]["GBSeq_organism"])


###### parsing GenBank records

handle = Entrez.egquery(term="Opuntia AND rpl16")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
     if row["DbName"]=="nuccore":
         print(row["Count"])

handle = Entrez.esearch(db="nuccore", term="Opuntia AND rpl16")
record = Entrez.read(handle)
gi_list = record["IdList"]
gi_list

gi_str = ",".join(gi_list)
handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")

text = handle.read()


#To get the records in a more Python-friendly form, we can use Bio.SeqIO 

from Bio import SeqIO
handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
records = SeqIO.parse(handle, "gb")

#### lineage of an organism

handle = Entrez.esearch(db="Taxonomy", term="Cypripedioideae")
record = Entrez.read(handle)
record["IdList"]


handle = Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
records = Entrez.read(handle)

records[0].keys()
records[0]["Lineage"]


############################### SWISS-PROT #############################

from Bio import ExPASy
handle = ExPASy.get_sprot_raw(myaccessionnumber)        # from ExPASy database
handle = gzip.open("myswissprotfile.dat.gz", "rt")      # from gzipped file
handle = open("myswissprotfile.dat")                    # local file
handle = urllib.urlopen("http://www.somelocation.org/data/someswissprotfile.dat")    # online file



from Bio import SwissProt
record = SwissProt.read(handle)

print(record.description)


from Bio import SwissProt
handle = open("uniprot_sprot.dat")
descriptions = [record.description for record in SwissProt.parse(handle)]
len(descriptions)
descriptions[:5]


############################## PROSITE #################################

from Bio.ExPASy import Prosite
handle = open("myprositefile.dat")
records = Prosite.parse(handle)

## multiple records
records = Prosite.parse(handle)
record = next(records)
record.accession

record = next(records)
record.accession


# scanning the database

sequence = "MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFTCRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN"

from Bio.ExPASy import ScanProsite
handle = ScanProsite.scan(seq=sequence)

result = ScanProsite.read(handle)
type(result)


handle = ScanProsite.scan(seq=sequence, lowscore=1)
result = ScanProsite.read(handle)
result.n_match


############################ EXPASY (enzymes) ###########################


from Bio import ExPASy
from Bio import SwissProt

accessions = ["O23729", "O23730", "O23731"]
records = []

for accession in accessions:
     handle = ExPASy.get_sprot_raw(accession)
     record = SwissProt.read(handle)
     records.append(record)


#from Bio import ExPASy
#import re

#handle = ExPASy.sprot_search_de("Orchid Chalcone Synthase")
#html_results = handle.read()
#if "Number of sequences found" in html_results:
#     ids = re.findall(r'HREF="/uniprot/(\w+)"', html_results)
#     else:
#     ids = re.findall(r'href="/cgi-bin/niceprot\.pl\?(\w+)"', html_results)


# Retrieving Prosite and Prosite documentation records

#from Bio import ExPASy
#handle = ExPASy.get_prosite_raw('PS00001')
#text = handle.read()
#print(text)

#record = Prosite.read(handle)             # parse it into a Bio.Prosite.Record object

#record = Prodoc.read(handle)   # parse it into a Bio.ExPASy.Prodoc.Record object:


#################################### KEGG ##########################################

from Bio.KEGG import Enzyme
records = Enzyme.parse(open("C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\ec_5.4.2.2.txt"))
record = list(records)[0]
record.classname


##### KEGG API

from Bio.KEGG import REST
from Bio.KEGG import Enzyme
request = REST.kegg_get("ec:5.4.2.2")
open("ec_5.4.2.2.txt", 'w').write(request.read())
records = Enzyme.parse(open("ec_5.4.2.2.txt"))
record = list(records)[0]
record.classname


# combination of querying the KEGG API. This will demonstrate how to extract a unique set of all # human pathway gene symbols which relate to DNA repair. 

from Bio.KEGG import REST

human_pathways = REST.kegg_list("pathway", "hsa").read()

# Filter all human pathways for repair pathways
repair_pathways = []
for line in human_pathways.rstrip().split("\n"):
    entry, description = line.split("\t")
    if "repair" in description:
        repair_pathways.append(entry)

# Get the genes for pathways and add them to a list
repair_genes = [] 
for pathway in repair_pathways:
    pathway_file = REST.kegg_get(pathway).read()  # query and read each pathway

    # iterate through each KEGG pathway file, keeping track of which section
    # of the file we're in, only read the gene in each pathway
    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section
        
        if current_section == "GENE":
            gene_identifiers, gene_description = line[12:].split("; ")
            gene_id, gene_symbol = gene_identifiers.split()

            if not gene_symbol in repair_genes:
                repair_genes.append(gene_symbol)

print("There are %d repair pathways and %d repair genes. The genes are:" % \
      (len(repair_pathways), len(repair_genes)))
print(", ".join(repair_genes))


