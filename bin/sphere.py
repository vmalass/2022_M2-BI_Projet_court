import math
import plotly.express as px
import pandas
import numpy

#Creation de la sphere avec 92 points
def fibonacci_sphere():
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
    points = []
    gold = numpy.pi * (3. - numpy.sqrt(5.))  # nombre d'or
    sx =[]
    sy =[]
    sz =[]
    for i in range(92):
        y = 1 - (i / float(92 - 1)) * 2

        radius = numpy.sqrt(1 - y * y)  # radius at y
        theta = gold * i  # golden angle increment
        x = numpy.cos(theta) * radius
        z = numpy.sin(theta) * radius
        rayon=1.7
        sx.append(rayon * x + 32.465)
        sy.append(rayon * y + 55.196)
        sz.append(rayon * z + 24.877)
        points.append((sx, sy, sz,)) 
    return sx, sy, sz



df=fibonacci_sphere()
fig = px.scatter_3d(df, x=df[0], y=df[1], z=df[2])
fig.show()