import math
import plotly.express as px
import pandas
import numpy


def fibonacci_sphere(n=92):

    points = []
    phi = numpy.pi * (3. - numpy.sqrt(5.))  # nombre d'or
    sx =[]
    sy =[]
    sz =[]
    for i in range(n):
        y = 1 - (i / float(n - 1)) * 2

        radius = numpy.sqrt(1 - y * y)  # radius at y
        theta = phi * i  # golden angle increment
        x = numpy.cos(theta) * radius
        z = numpy.sin(theta) * radius
        rayon=1.7
        sx.append(rayon * x + 32.465)
        sy.append(rayon * y + 55.196)
        sz.append(rayon * z + 24.877)
        points.append((sx, sy, sz,)) 
    return sx,sy,sz



df=fibonacci_sphere()
fig = px.scatter_3d(df, x=df[0], y=df[1], z=df[2])
fig.show()