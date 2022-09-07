#Biopython parser pdb
import math
from Bio.PDB import *
import numpy
from scipy.spatial.distance import *

## 1) Parser pdb pour obtenir les atomes et les coordonnées

parser = PDBParser(PERMISSIVE=1)  #il ignore les erreurs
structure_id = "7YL9"
filename = "7YL9.pdb"
structure = parser.get_structure(structure_id, filename)

residue = []
atom_id = []
atom_co = []
for model in structure:
    for chain in model:
        for residue in chain:
            residue.append(residue)
            for atom in residue:
                atom_id.append(atom)
                atom_co.append(atom.get_coord())

#print(atom_id)
#print(atom_co)

## 2) Calcule de la matrice de distance entre chaque atome

atom_co_array=numpy.array(atom_co) 
distances_array = squareform(pdist(atom_co_array))

print(distances_array.shape)