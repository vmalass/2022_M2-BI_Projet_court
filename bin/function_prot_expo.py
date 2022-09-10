import math
from readline import append_history_file
from Bio.PDB import *
import numpy
from scipy.spatial.distance import pdist, squareform

def fibonacci_sphere(coordonnee, vdw_ray):
    """
    The fibonacci_sphere function generates a sphere of points using the Fibonacci sequence.
    The function takes two arguments:
        coordonnee : The center of the sphere (x,y,z)
        vdw_ray : The radius of the sphere in angstroms
    
    :param coordonnee: Define the center of the sphere
    :param vdw_ray: Define the radius of the sphere that will be used to generate points
    :return: A list of points on the surface of a sphere
    :doc-author: Trelent
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
    The distance_euclidienne function calculates the distance between two points in 3D space.
    It takes as input a list of coordinates for one point and a list of coordinates for another point.
    The function returns the distance between those two points.
    
    :param co_pts_sphere: Calculate the distance between a point and its neighbor
    :param co_voisin: Store the coordinates of a point in the sphere
    :return: The distance between two points
    :doc-author: Trelent
    """
    distance = math.sqrt((co_pts_sphere[0]-co_voisin[0])**2+(co_pts_sphere[1]-co_voisin[1])**2+(co_pts_sphere[2]-co_voisin[2])**2)
    return distance