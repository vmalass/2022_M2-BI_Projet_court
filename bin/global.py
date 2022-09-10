import math
from readline import append_history_file
from Bio.PDB import *
import numpy
from scipy.spatial.distance import pdist, squareform

## creation sphere centre sur l'atomes avec ajout du rayon de van der waals de l'h2o
def fibonacci_sphere(coordonnee, vdw_ray):
    """
    The fibonacci_sphere function generates a sphere of points on the surface of a
    fibonacci sphere. The function takes two arguments, coordonnee and vdw_ray.
    coordonnee is an array containing the x, y and z coordinates for the center of 
    the sphere. vdw_ray is an integer representing the radius in angstroms to be used 
    for each point on the fibonacci sphere.
    
    :param coordonnee: Define the center of the sphere
    :param vdw_ray: Set the radius of the sphere
    :return: A list of points
    """
    gold = numpy.pi * (3. - numpy.sqrt(5.))  # nombre d'or
    sx =[]
    sy =[]
    sz =[]
    xyz = []
    for i in range(92):
        y = 1 - (i / float(92 - 1)) * 2
        radius = numpy.sqrt(1 - y * y)  # radius at y
        theta = gold * i  # golden angle increment
        x = numpy.cos(theta) * radius
        z = numpy.sin(theta) * radius
        rayon=vdw_ray+(1.7*2)
        sx = (rayon * x + coordonnee[0])
        sy = (rayon * y + coordonnee[1])
        sz = (rayon * z + coordonnee[2])
        xyz.append([sx, sy, sz])
    points = numpy.array(xyz) 
    return points

def distance_euclidienne(co_pts_sphere, co_voisin):
    """
    The distance_euclidienne function computes the Euclidean distance between two points.
    
    
    
    :param co_pts_sphere: Get the coordinates of the point in a sphere
    :param co_voisin: Store the coordinates of a point in the sphere
    :return: The distance between the two points
    :doc-author: Trelent
    """

    distance = math.sqrt((co_pts_sphere[0]-co_voisin[0])**2+(co_pts_sphere[1]-co_voisin[1])**2+(co_pts_sphere[2]-co_voisin[2])**2)
    return distance

## Dictionnaire des rayons de Vann der Waals
vdw_ray = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "P": 1.8, "S":1.8}

## Parser pdb pour obtenir les atomes et les coordonnées
parser = PDBParser(PERMISSIVE=1)  #il ignore les erreurs
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

atom_co = numpy.array(atom_co)

## Calcule de la matrice de distance entre chaque atome en data frame

distances_atom = squareform(pdist(atom_co)) #calcule des distances entre chaque atome et archivage dans un df

## Recherche des atomes voisins à 10A
point_surface = []
indexe_voisin = []
distance = []

for id, coor in zip(range(len(atom_id)), atom_co):
    indexe_voisin = numpy.where((distances_atom[:][id] <10) & (distances_atom[:][id]  >0))[0] #extrait les indexes ou la distance entre 2 atomes est inf a 10A
    point_sphere = fibonacci_sphere(coor,vdw_ray[atom_id[id]])
    point_surface_atom = 0
    for point in point_sphere:
        flag_break=True
        for voisin in indexe_voisin:
            distance = distance_euclidienne(point, atom_co[voisin])
            if distance<vdw_ray[atom_id[voisin]]:
                flag_break=False
                break
        if flag_break:
            point_surface_atom +=1
    point_surface.append(point_surface_atom)
