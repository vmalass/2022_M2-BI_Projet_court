from Bio.PDB import *
import numpy as np
from scipy.spatial.distance import pdist, squareform
import function_prot_expo as fpe
import sys

# Import pdb file.
# try:
#     file_name = sys.argv[1]
# except:
#     print('Enter the file.pdb')
#     sys.exit("Usage:python test.py file.pdb")

# Constant.
H20_RAY = 1.7
CST_SPHERE = 92

# Dictionary of Van der Waals rays.
vdw_ray = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "P": 1.8, "S": 1.8}

# Dictionary of residue.
residue_aera = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLU": 223,
 "GLN": 225, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,"MET": 224, 
 "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}

# Parser pdb for obtain atoms' coordonates, identifications and residues.
parser = PDBParser(PERMISSIVE=1)  # Allows errors to be ignored.
# structure_id = file_name.replace(".pdb", "")
# filename = file_name
structure_id = '3i40'
filename = '3i40.pdb'
structure = parser.get_structure(structure_id, filename)

resi = []
unique_residue = []
atom_id = []
atom_co = []

for chain in structure[0]:  # Choose first model in pdb
    for residue in chain:
        # Choose just ATOM.
        if is_aa(residue):       ## augmente aera avec perte de 35 atomes) ##
            unique_residue.append(str(residue).split(" ")[1])
            for atom in residue:
                # Removes the d (disordered atoms).
                if str(atom)[6] != "d":
                    # Obtain the atoms' identifications.
                    atom_id.append(str(atom)[6])
                    # Obtain the atoms' coordonates.
                    atom_co.append(atom.get_coord())
                    # Obtain just the residue.
                    resi.append(str(atom.get_parent()).split(" ")[1])

atom_co = np.array(atom_co)  # Created a array matrix.

# Calculation of the distance matrix between each atom.
distances_atom = squareform(pdist(atom_co))

# Search neighbours' atoms.
surface_point = []
surface_atom_expose = []
ratio = []
index_neighbour = []
distance = []
total_surface_prot = 0
for id, coor in zip(range(len(atom_id)), atom_co):
    # neighbours extraction.
    index_neighbour = np.where((distances_atom[:][id] < 10) &
                               (distances_atom[:][id] != 0))[0]
    # Creation of a sphere around the atom.
    point_sphere = fpe.fibonacci_sphere(coor, vdw_ray[atom_id[id]])
    surface_point_atom = 0  # Counter for the exposed points of the sphere.
    for point in point_sphere:
        flag_break = True
        for neighbour in index_neighbour:
            # Calculation of the distance between a point on the sphere and the
            # neighbouring atom.
            distance = fpe.distance_euclidienne(point, atom_co[neighbour])
            if distance < (vdw_ray[atom_id[neighbour]] + H20_RAY):
                flag_break = False
                break
        if flag_break:
            surface_point_atom += 1
        # Calculation the ratio exposure point and surface in angstroms per atoms
    total_surface_prot += (4 * np.pi * (vdw_ray[atom_id[id]] + H20_RAY)**2)
    ratio = (surface_point_atom / CST_SPHERE) * (4 * np.pi * (vdw_ray[atom_id[id]] + H20_RAY)**2)
    # Storage of the points of the free sphere in a list.
    surface_point.append(surface_point_atom)
    surface_atom_expose.append(ratio)

# Calculation the surface exposed per atoms.
area = sum(surface_atom_expose)

# Calculation the surface exposed per residues.
area_residue = {}
for key, value in zip(resi, surface_atom_expose):
    area_residue[key] = area_residue.get(key, 0) + value

relative_area=0
for res in unique_residue:
    relative_area += area_residue[res]/residue_aera[res]

access_100 = 100 * area / total_surface_prot

