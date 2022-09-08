#Biopython parser pdb
import math
from Bio.PDB import *
import numpy
from scipy.spatial.distance import pdist, squareform
import pandas

## 1) Parser pdb pour obtenir les atomes et les coordonnées

parser = PDBParser(PERMISSIVE=1)  #il ignore les erreurs
structure_id = "7YL9"
filename = "7YL9.pdb"
structure = parser.get_structure(structure_id, filename)

resi = []
atom_id = []
atom_co = []
for model in structure:
    for chain in model:
        for residue in chain:
            resi.append(residue)
            for atom in residue:
                atom_id.append(atom)
                atom_co.append(atom.get_coord())

#print(atom_id)
#print(atom_co)

atom_co_df = pandas.DataFrame(atom_co, columns=['x_co','y_co','z_co'], index=atom_id) #df des atomes et de leurs coordonées
#atom_co_df.insert(1, "residue", resi, allow_duplicates=True) 
#print(atom_co_df)

## 2) Calcule de la matrice de distance entre chaque atome en data frame

atom_co_array=numpy.array(atom_co) #matrice avec les coordonnées
distances_df = pandas.DataFrame(squareform(pdist(atom_co_array)), columns=atom_id, index=atom_id) #calcule des distances entre chaque atome et archivage dans un df

print(distances_df)

## 3) Recherche des atomes voisins à 10A

atom_voison=[]
for i in range(distances_df):
    numpy.where(distances_df[:][i]<10) #extrait les indexes ou la distance entre 2 atomes est inf a 10A
    











# 3) Dictionnaire des rayons de Vann der Waals

vdw_rad = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "P": 1.8, "S":1.8, "H20": 1.7}