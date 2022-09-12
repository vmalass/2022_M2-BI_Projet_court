from Bio.PDB import *
import numpy as np
from scipy.spatial.distance import pdist, squareform

# Constant
H2O_RAY = 1.7
CST_SPHERE = 92


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
