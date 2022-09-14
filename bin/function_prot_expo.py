from Bio.PDB import *
import numpy as np

# Constant
H2O_RAY = 1.7
CST_SPHERE = 92

def parser_filter(parser_object):
    """
    The parser_filter function takes a parser object (in Bio.PDB) and
    returns the residues, unique residues, atomic identification numbers
    and coordinates.

    :param parser_object: the parser object
    :return: residues, unique residues, atom identification numbers
             and coordinates.
    """
    resi = []
    unique_residue = []
    atom_id = []
    atom_co = []

    for chain in parser_object[0]:  # Choose first model in pdb
        for residue in chain:
            # Choose just ATOM.
            if is_aa(residue):
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
    return resi, unique_residue, atom_id, atom_co



def fibonacci_sphere(coordonnee, vdw_ray):
    """
    The fibonacci_sphere function generates a sphere of points using the
    Fibonacci sequence around an atom.
    The function takes two arguments:
        coordonnee : Position of the atom in the protein (x,y,z)
        vdw_ray : The Van der Waals radius of the atom in angstroms

    :param coordonnee : Set the centre of the sphere
    :param vdw_ray : Define the radius of the sphere which will be used to
    generate points
    :return : A list of points on the surface of a sphere
    """
    gold = np.pi * (3. - np.sqrt(5.))  # Golden number.
    sx = []
    sy = []
    sz = []
    xyz = []
    for i in range(CST_SPHERE):
        y = 1 - (i / float(CST_SPHERE - 1)) * 2
        radius = np.sqrt(1 - y * y)
        theta = gold * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        rayon = vdw_ray + (H2O_RAY)

        sx = (rayon * x + coordonnee[0])
        sy = (rayon * y + coordonnee[1])
        sz = (rayon * z + coordonnee[2])
        xyz.append([sx, sy, sz])

    points = np.array(xyz)
    return points


def distance_euclidienne(co_pts_sphere, co_voisin):
    """
    The function Euclidean_distance calculates the distance between
    two points in 3D space.
    It takes as input two lists :
        co_pts_sphere : coordinates of the points of the sphere.
        co_voisin : coordinates of the neighbouring atom to the sphere.

    :param co_pts_sphere : Coordinate of a point of the sphere
    :param co_neighbour : Coordinate of the neighbouring atom
    :return : The distance between two points
    """
    distance = np.sqrt((co_pts_sphere[0] - co_voisin[0])**2 +
                       (co_pts_sphere[1] - co_voisin[1])**2 +
                       (co_pts_sphere[2] - co_voisin[2])**2)
    return distance


def residue_area(indexe_residue, expose_surface_atom):
    """
    The residue_area function takes in a list of residue and the corresponding
    exposure surface area for each atom in residues. It returns a dictionary
    with the residues as keys and their corresponding exposure surface area
    as values.

    :param indexe_residue: The residue index.
    :param expose_surface_atom: The area of each atoms in residue.
    :return: A dictionary with the key being the residue index and value
    being the total surface area of that residue exposed to solvent
    """
    area_residue = {}
    for key, value in zip(indexe_residue, expose_surface_atom):
        area_residue[key] = area_residue.get(key, 0) + value
    return area_residue


def area_relative(chain_residue, area_residue):
    """
    The area_relative function takes one liste and dictionarie as input. The
    list contains each residue in a protein, and the dictionary contains the
    total area of each residue in a protein. The function then returns an
    integer representing how much more exposed relative surface exposition
    in protein.

    :param chain_residue: Specify which residues to include in the calculation
    of relative area
    :param area_residue: Calculate the relative area of each residue in a chain
    :return: The relative area of the chain to the total surface area
    """
    relative_area = 0
    for res in chain_residue:
        relative_area += area_residue[res]/area_residue[res]
    return relative_area
