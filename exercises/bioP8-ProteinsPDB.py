
######################## Reading and writing crystal structure files ##########################

##### determining secondary structure

from Bio.PDB.PDBParser import PDBParser
parser = PDBParser()
structure_id = "1AGR"
filename = "C:\\Users\\user\\Desktop\\tertiary\\PYTprogramming\\BioPython\\datasets\\protein_align\\1AGR_A.pdb"

# Iterating through all atom
p = PDBParser()
structure = p.get_structure(structure_id, filename)
for model in structure:
     for chain in model:
         for residue in chain:
             for atom in residue:
                 print(atom)



## writing files
io = PDBIO()
io.set_structure(s)
io.save('out.pdb')


# The Structure object is at the top of the hierarchy. Its id is a user given string
# The id of the Model object is an integer, which is derived from the position of the model in the # parsed file (they are automatically numbered starting from 0).
# The id of a Chain object is derived from the chain identifier in the PDB/mmCIF file, and is a # single character (typically a letter)
# A residue id is a tuple with three elements: The hetero-field (hetfield) ;
#                                             : The sequence identifier (resseq), an integer describing the position of the residue in the chain (e.g., 100);
#                                             :  insertion code (icode); a string



########################### Accessing the Protein Data Bank  ##################################

pdbl = PDBList()
pdbl.retrieve_pdb_file('1FAT')
