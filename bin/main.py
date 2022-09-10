from Bio.PDB import *
import numpy as np
from scipy.spatial.distance import pdist, squareform
import function_prot_expo as fpe

## Dictionnaire des rayons de Vann der Waals
vdw_ray = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "P": 1.8, "S":1.8}

## Parser pdb pour obtenir les atomes et les coordonnées
parser = PDBParser(PERMISSIVE=1)  #permet d'ignorer les erreurs
structure_id = "3i40"
filename = "./data/3i40.pdb"
structure = parser.get_structure(structure_id, filename)

resi = []
atom_id = []
atom_co = []
for model in structure:
    for chain in model:
        for residue in chain:
            resi.append(residue)
            for atom in residue:
                if str(atom)[6] != "d": #supprime les d
                    atom_id.append(str(atom)[6])
                    atom_co.append(atom.get_coord())

atom_co = np.array(atom_co)

## Calcule de la matrice de distance entre chaque atome en data frame
distances_atom = squareform(pdist(atom_co)) 

## Recherche des atomes voisins à 10A
point_surface = []
indexe_voisin = []
distance = []
for id, coor in zip(range(len(atom_id)), atom_co):
    indexe_voisin = np.where((distances_atom[:][id] <10) & (distances_atom[:][id]  >0))[0] #extrait les indexes ou la distance entre 2 atomes est inf a 10A
    point_sphere = fpe.fibonacci_sphere(coor,vdw_ray[atom_id[id]])
    point_surface_atom = 0
    for point in point_sphere:
        flag_break=True
        for voisin in indexe_voisin:
            distance = fpe.distance_euclidienne(point, atom_co[voisin])
            if distance<vdw_ray[atom_id[voisin]]:
                flag_break=False
                break
        if flag_break:
            point_surface_atom +=1
    point_surface.append(point_surface_atom)
print(point_surface)