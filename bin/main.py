from Bio.PDB import *
import numpy as np
from scipy.spatial.distance import pdist, squareform
import function_prot_expo as fpe

# Dictionary of Van der Waals rays.
vdw_ray = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "P": 1.8, "S": 1.8}

# Parser pdb for obtain atoms' coordonates.
parser = PDBParser(PERMISSIVE=1)  # Allows errors to be ignored.
structure_id = "3i40"
filename = "./data/3i40.pdb"
structure = parser.get_structure(structure_id, filename)

resi = []
atom_id = []
atom_co = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                # Removes the d (disordered atoms).
                if str(atom)[6] != "d":
                    # Obtain the atoms' identifications.
                    atom_id.append(str(atom)[6])
                    # Obtain the atoms' coordonates.
                    atom_co.append(atom.get_coord())
                    # Obtain the residue.
                    resi.append(atom.get_parent())

atom_co = np.array(atom_co)  # Created a array matrix.

# Calculation of the distance matrix between each atom.
distances_atom = squareform(pdist(atom_co))

# Search neighbours' atoms.
surface_point = []
index_neighbour = []
distance = []
for id, coor in zip(range(len(atom_id)), atom_co):
    # neighbours extraction.
    index_neighbour = np.where((distances_atom[:][id] < 10) &
                               (distances_atom[:][id] > 0))[0]
    # Creation of a sphere around the atom.
    point_sphere = fpe.fibonacci_sphere(coor, vdw_ray[atom_id[id]])
    surface_point_atom = 0  # Counter for the exposed points of the sphere.
    for point in point_sphere:
        flag_break = True
        for neighbour in index_neighbour:
            # Calculation of the distance between a point on the sphere and the
            # neighbouring atom.
            distance = fpe.distance_euclidienne(point, atom_co[neighbour])
            if distance < vdw_ray[atom_id[neighbour]]:
                flag_break = False
                break
        if flag_break:
            surface_point_atom += 1
    # Storage of the points of the free sphere in a list.
    surface_point.append(surface_point_atom)
