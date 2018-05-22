from numpy import array
from basegeometry import Triangle
coords = array([(2, 1), (2, 8), (7, 1)])
triangle = Triangle(coords)
triangle.packCircles(want2plot=True)